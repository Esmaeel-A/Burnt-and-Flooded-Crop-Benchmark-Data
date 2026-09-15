#!/usr/bin/env python3
"""
Stage A:  build the burn EVENT LAYER.

Finds candidate agricultural fires(WHERE and WHEN a fire
happened) and confirms each against three independent fire
products. 

    MCD64A1 monthly burned area years 2020 2021 2022 2023 2024 2025
      -connected-component clustering into discrete fire patches
      -vectorise, filter by area
      -cropland (ESA WorldCover)
      -3-way verification
      -spatiotemporal de-duplication
      -dataset/manifests/events.csv 

 all three must agree, or the event is discarded:

    MODIS MCD64A1   burned area   >= 20% of the perimeter   (locates the event)
    VIIRS VNP64A1   burned area   >= 20% of the perimeter   (independent sensor)
    VIIRS VNP14A1   active fire   >= 3 thermal detections   (independent signal)

    MCD64A1 is the source, so it is not independent of itself; the two VIIRS products
    are the independent confirmations (`n_indep = 2`). Events that fail any test are
    dropped.

    Cells are visited (month, year, continent) year-major, so each pass touches every
    year. To prevent high burn region from dominating the quota before the 
    the later years are visited. 
    
    Re-running appends to an existing events.csv, skipping anything already present by
    id and by spatiotemporal proximity.

"""
from __future__ import annotations
import os
import sys, csv, json, math, argparse, datetime as dt
from pathlib import Path

import ee

if not os.environ.get("EE_PROJECT"):
    sys.exit("Set EE_PROJECT in the environment (export EE_PROJECT=your-project-id).")
ee.Initialize(project=os.environ["EE_PROJECT"])

PROJ = Path.cwd()
OUT = PROJ / "dataset" / "manifests"

COUNTRIES = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017")

# --- gates -----------------------------------------------------------------
# MIN_AREA is larger than one chip (512 px @ 10 m = 5.12 km, ~26 km2) 
MIN_AREA, MAX_AREA = 50.0, 20000.0      # km2
CROP_MIN = 0.15                         # cropland share of the perimeter
VERIFY_MIN, AF_MIN = 0.20, 3            # burned-area fraction / active-fire pixel count
MIN_PATCH_PX = 20                       # >=20 MODIS px @500 m before the km2 filter
VERIFY_PAD_DAYS = 45                    # widen the window when confirming
DEDUP_KM, DEDUP_DAYS = 25.0, 30         # same fire if this close in space AND time

CONTINENTS = [
    ("NorthAmerica", (-170, 12, -50, 84)), ("CentralAmericaCaribbean", (-95, 5, -58, 27)),
    ("SouthAmerica", (-93, -57, -32, 13)), ("Europe", (-25, 34, 45, 72)),
    ("Africa", (-20, -37, 52, 38)), ("MiddleEast", (34, 12, 63, 42)),
    ("SouthAsia", (60, 5, 98, 38)), ("EastAsia", (98, 18, 150, 55)),
    ("SoutheastAsia", (92, -11, 142, 21)), ("CentralAsiaRussia", (45, 40, 180, 78)),
    ("Oceania", (110, -50, 180, -10)),
]

COLS = ["event_id", "source", "tier", "bbox_w", "bbox_s", "bbox_e", "bbox_n",
        "date_start", "date_end", "event_date_start", "event_date_end", "year",
        "continent", "country", "place", "area_km2_reported",
        "match_modis", "match_viirs", "match_af", "n_indep", "aoi_tightened",
        "modis_frac", "viirs_frac", "af_count", "burn_km2", "crop_frac", "n_chips", "notes"]


def continent_of(lon, lat):
    for name, (w, s, e, n) in CONTINENTS:
        if w <= lon <= e and s <= lat <= n:
            return name
    return "Other"


def cropland():
    """ESA WorldCover v200 class 40 = cropland."""
    return ee.ImageCollection("ESA/WorldCover/v200").first().select("Map").eq(40).rename("crop")


def iso(ms, shift_days=0):
    return (dt.datetime.utcfromtimestamp(ms / 1000)
            + dt.timedelta(days=shift_days)).strftime("%Y-%m-%d")


def doy_to_ms(year, doy):
    return int(dt.datetime(int(year), 1, 1).timestamp() * 1000) + int((int(doy) - 1) * 86400000)


def bbox_of(coords):
  
    xs, ys = [], []

    def rec(x):
        if x and isinstance(x[0], (int, float)):
            xs.append(x[0]); ys.append(x[1])
        else:
            for y in x:
                rec(y)
    rec(coords)
    return round(min(xs), 5), round(min(ys), 5), round(max(xs), 5), round(max(ys), 5)


# --------------------------------------------------------------------------- #
def fire_patches(region, d0, d1):

    burn_date = ee.ImageCollection("MODIS/061/MCD64A1").filterDate(d0, d1).select("BurnDate").max()
    burned = burn_date.gt(0).selfMask()
    sizes = burned.connectedPixelCount(256, True)
    patches = burned.connectedComponents(ee.Kernel.plus(1), 256).updateMask(sizes.gte(MIN_PATCH_PX))
    vec = patches.select("labels").reduceToVectors(
        geometry=region, scale=500, geometryType="polygon", eightConnected=True,
        maxPixels=1e9, bestEffort=True, labelProperty="lab")

    crop = cropland()

    def enrich(f):
        g = f.geometry()
        centre = g.centroid(1000).coordinates()
        span = burn_date.reduceRegion(ee.Reducer.minMax(), g, scale=500,
                                      maxPixels=1e7, bestEffort=True)
        return f.set("km2", g.area(1000).divide(1e6),
                     "lon", centre.get(0), "lat", centre.get(1),
                     "doy_min", span.get("BurnDate_min"), "doy_max", span.get("BurnDate_max"),
                     "crop_frac", crop.reduceRegion(ee.Reducer.mean(), g, scale=200,
                                                    maxPixels=1e8, bestEffort=True).get("crop"))

    return (vec.map(enrich)
            .filter(ee.Filter.gte("km2", MIN_AREA))
            .filter(ee.Filter.lte("km2", MAX_AREA))
            .filter(ee.Filter.gte("crop_frac", CROP_MIN)))


def verify(ig_ms, fd_ms, geom):
    """Confirm one patch against MODIS burned area, VIIRS burned area and VIIRS active fire.
    """
    w0 = ee.Date(ig_ms).advance(-VERIFY_PAD_DAYS, "day")
    w1 = ee.Date(fd_ms).advance(VERIFY_PAD_DAYS, "day")
    modis = ee.ImageCollection("MODIS/061/MCD64A1").filterDate(w0, w1).select("BurnDate").max().gt(0).unmask(0)
    viirs = ee.ImageCollection("NASA/VIIRS/002/VNP64A1").filterDate(w0, w1).select("Burn_Date").max().gt(0).unmask(0)
    active = ee.ImageCollection("NASA/VIIRS/002/VNP14A1").filterDate(w0, w1).select("FireMask").max().gte(7).unmask(0)
    country = COUNTRIES.filterBounds(geom.centroid(1000)).first()
    return ee.Dictionary({
        "modis_frac": modis.reduceRegion(ee.Reducer.mean(), geom, scale=250, maxPixels=1e8, bestEffort=True).get("BurnDate"),
        "viirs_frac": viirs.reduceRegion(ee.Reducer.mean(), geom, scale=250, maxPixels=1e8, bestEffort=True).get("Burn_Date"),
        "af_count": active.reduceRegion(ee.Reducer.sum(), geom, scale=1000, maxPixels=1e8, bestEffort=True).get("FireMask"),
        "country": ee.Algorithms.If(country, country.get("country_na"), ""),
    }).getInfo()


# --------------------------------------------------------------------------- #
class Deduper:
    """remove duplicates: a patch that is the same physical fire as one already accepted.

    Since the same fire can surface in adjacent months, and two
    clusterings of overlapping windows produce different ids for one burn. Two events
    are the same fire when their centroids are within DEDUP_KM *and* their start dates
    within DEDUP_DAYS. Not just Ids. 
    """

    def __init__(self, existing):
        self.pts = []
        for r in existing:
            try:
                lon = (float(r["bbox_w"]) + float(r["bbox_e"])) / 2.0
                lat = (float(r["bbox_s"]) + float(r["bbox_n"])) / 2.0
                self.pts.append((lon, lat, self.day(r["event_date_start"])))
            except (ValueError, KeyError):
                continue

    @staticmethod
    def day(s):
        return dt.date.fromisoformat(str(s)[:10]).toordinal()

    def seen(self, lon, lat, day):
        for xlon, xlat, xday in self.pts:
            if abs(day - xday) > DEDUP_DAYS:
                continue
            dy = (lat - xlat) * 111.0
            dx = (lon - xlon) * 111.0 * math.cos(math.radians(lat))
            if math.hypot(dx, dy) <= DEDUP_KM:
                return True
        return False

    def add(self, lon, lat, day):
        self.pts.append((lon, lat, day))


def event_id(lon, lat, year):

    return f"MC{year}{abs(hash((round(lon, 3), round(lat, 3), year))) % 10 ** 8:08d}"


def build(n, years, per_call, per_cell, fresh):
    OUT.mkdir(parents=True, exist_ok=True)
    events_csv = OUT / "events.csv"

    existing = []
    if events_csv.exists() and not fresh:
        existing = list(csv.DictReader(open(events_csv)))
        print(f"resuming — {len(existing)} events already in the layer")
    have = {r["event_id"] for r in existing}
    dedup = Deduper(existing)
    kept = []

    conts = [(nm, bx) for nm, bx in CONTINENTS]
    cells = [(y, mi, nm, bx) for mi in range(12) for y in years for nm, bx in conts]
    print(f"scanning {len(cells)} (month x year x continent) cells for up to {n} events ...")

    for (y, mi, nm, bx) in cells:
        if len(kept) >= n:
            break
        d0 = dt.date(y, mi + 1, 1).isoformat()
        d1 = dt.date(y + (mi == 11), ((mi + 1) % 12) + 1, 1).isoformat()
        region = ee.Geometry.Rectangle(list(bx))
        try:
            feats = fire_patches(region, d0, d1).limit(per_call).getInfo()["features"]
        except Exception:
            continue                       # a cell that times out is skipped, not fatal
        if not feats:
            continue

        cell_n = 0
        for f in feats:
            if len(kept) >= n or cell_n >= per_cell:
                break
            p = f["properties"]
            if p.get("crop_frac") is None or p.get("doy_min") is None:
                continue
            ig = doy_to_ms(y, p["doy_min"])
            fd = doy_to_ms(y, p["doy_max"])
            eid = event_id(p["lon"], p["lat"], y)
            if eid in have or dedup.seen(p["lon"], p["lat"], Deduper.day(iso(ig))):
                continue
            try:
                v = verify(ig, fd, ee.Geometry(f["geometry"]))
            except Exception:
                continue
            mf = v["modis_frac"] or 0
            vf = v["viirs_frac"] or 0
            afc = int(v["af_count"] or 0)
            if not (mf >= VERIFY_MIN and vf >= VERIFY_MIN and afc >= AF_MIN):
                continue                   # all three must confirm

            have.add(eid)
            dedup.add(p["lon"], p["lat"], Deduper.day(iso(ig)))
            bw, bs, be, bn = bbox_of(f["geometry"]["coordinates"])
            cont = continent_of(p["lon"], p["lat"])
            kept.append({
                "event_id": eid, "source": "mcd64-cluster", "tier": "A",
                "bbox_w": bw, "bbox_s": bs, "bbox_e": be, "bbox_n": bn,
                "date_start": iso(ig, -VERIFY_PAD_DAYS), "date_end": iso(fd, VERIFY_PAD_DAYS),
                "event_date_start": iso(ig), "event_date_end": iso(fd), "year": y,
                "continent": cont, "country": v["country"], "place": v["country"],
                "area_km2_reported": round(p["km2"], 2),
                "match_modis": True, "match_viirs": True, "match_af": True,
                "n_indep": 2, "aoi_tightened": True,
                "modis_frac": round(mf, 3), "viirs_frac": round(vf, 3), "af_count": afc,
                "burn_km2": round(p["km2"] * mf, 2), "crop_frac": round(p["crop_frac"], 3),
                "n_chips": "", "notes": ""})
            cell_n += 1
            print(f"  [{len(kept):4d}/{n}] {eid} {y} {cont[:12]:12s} {str(v['country'])[:14]:14s} "
                  f"crop={100*p['crop_frac']:3.0f}% V={100*vf:3.0f}% AF={afc:3d} area={p['km2']:.0f}km2")

    rows = existing + kept
    if not rows:
        print("no events found — nothing written"); return 1

    with open(events_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in COLS})

    feats = [{"type": "Feature",
              "geometry": {"type": "Polygon", "coordinates": [[
                  [r["bbox_w"], r["bbox_s"]], [r["bbox_e"], r["bbox_s"]],
                  [r["bbox_e"], r["bbox_n"]], [r["bbox_w"], r["bbox_n"]],
                  [r["bbox_w"], r["bbox_s"]]]]},
              "properties": {k: r[k] for k in ("event_id", "continent", "country", "crop_frac", "year")}}
             for r in rows]
    (OUT / "events.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats}))

    (OUT / "run_report_events.json").write_text(json.dumps({
        "stage": "event_layer", "method": "MCD64A1 connected-component clustering (GlobFire method)",
        "n_new": len(kept), "n_total": len(rows), "years": sorted(years),
        "gates": {"min_area_km2": MIN_AREA, "max_area_km2": MAX_AREA, "crop_min": CROP_MIN,
                  "verify_min_frac": VERIFY_MIN, "active_fire_min_px": AF_MIN,
                  "verify": "MODIS MCD64A1 AND VIIRS VNP64A1 AND VIIRS VNP14A1"},
        "by_continent": {c: sum(1 for r in rows if r["continent"] == c)
                         for c in sorted({r["continent"] for r in rows})},
        "by_year": {str(y): sum(1 for r in rows if str(r["year"]) == str(y))
                    for y in sorted({r["year"] for r in rows}, key=str)},
    }, indent=2))

    print(f"\nevent layer -> {events_csv}  ({len(kept)} new, {len(rows)} total)")
    print("  by year:     ", {y: sum(1 for r in kept if r["year"] == y)
                              for y in sorted({r["year"] for r in kept})})
    print("  by continent:", {c: sum(1 for r in kept if r["continent"] == c)
                              for c in sorted({r["continent"] for r in kept})})
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description="Stage A — build the burn event layer")
    ap.add_argument("--n", type=int, default=1600, help="events to add this run")
    ap.add_argument("--years", nargs="+", type=int, default=[2020, 2021, 2022, 2023, 2024, 2025])
    ap.add_argument("--per-call", type=int, default=40, help="patches fetched per cell")
    ap.add_argument("--per-cell", type=int, default=8,
                    help="cap per (year, month, continent) — stops one region monopolising the quota")
    ap.add_argument("--fresh", action="store_true", help="start a new layer instead of appending")
    a = ap.parse_args(argv)
    return build(a.n, a.years, a.per_call, a.per_cell, a.fresh)


if __name__ == "__main__":
    sys.exit(main())
