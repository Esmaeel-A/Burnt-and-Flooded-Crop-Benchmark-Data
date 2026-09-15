#!/usr/bin/env python3
"""
Stage B: sample burned-cropland chips and write the full chip package.

 For each event in the layer, finds the windows containing the most burned cropland and
 downloads the label, both Sentinel-2 , and both Sentinel-1. 

 PER CHIP (512 px @ 10 m = 5.12 km, in the event's UTM zone)

  <chip>_label.tif     5-band uint8
  <chip>_s2_pre.tif    12-band float32   Sentinel-2 before the burn
  <chip>_s2_post.tif   12-band float32   Sentinel-2 after
  <chip>_s1_pre.tif     2-band float32   Sentinel-1 VV+VH before   
  <chip>_s1_post.tif    2-band float32   Sentinel-1 VV+VH after
  <chip>.json          sidecar metadata

  Band layout

      S2:  B2 B3 B4 B5 B6 B7 B8 B8A B9 B11 B12 NDVI        NDVI = ND(B8A, B4)
      S1:  VV VH  (dB)

 THE LABEL is  Sentinel-2 dNBR, clipped to the event perimeter.

    NBR  = (B8 - B12) / (B8 + B12)
    dNBR = (NBR_pre - NBR_post) x 1000

 Since one scene rarely covers a whole chip, and cloud gaps leave holes — both produce black
  chips. Candidate scenes are masked, sorted, and mosaicked with the best on top: it
  supplies most pixels and the rest is gap-filled by the next best. S2 sorts by clear
  fraction, S1 by temporal distance. The top scene's date is reported as the nominal
  acquisition.

  the MCD64A1 patch here is
  used only to bound where a burn may be claimed. It never masks the label. 

  A cropland pixel is burned when dNBR >= 100 and it falls inside the perimeter.
  Severity follows USGS / Key & Benson (2006).

 

"""
from __future__ import annotations
import os
import sys, csv, json, argparse, datetime as dt, urllib.request
from pathlib import Path

import numpy as np
import rasterio
import rasterio.warp as rwarp
from rasterio.transform import xy as rxy

import ee

if not os.environ.get("EE_PROJECT"):
    sys.exit("Set EE_PROJECT in the environment (export EE_PROJECT=your-project-id).")
ee.Initialize(project=os.environ["EE_PROJECT"])

PROJ = Path.cwd()
MAN = PROJ / "dataset" / "manifests"
CHIPS = PROJ / "dataset" / "processed" / "chips"

SIZE, RES = 512, 10                       # 512 px @ 10 m = 5.12 km
HALF = SIZE * RES / 2
CHIPS_PER_EVENT = 5


MIN_CROP_FRAC = 0.20
MIN_BURN_FRAC = 0.02

SEV_T = {"low": 100, "modlow": 270, "modhigh": 440, "high": 660}   # dNBR x1000
CS_THRESH = 0.60                          # Cloud Score+ clear threshold
CONT_CAP = 700                            # max events chipped per continent

# flood's exact stack and order — NDVI is derived and always last
S2_BANDS = ["B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B9", "B11", "B12", "NDVI"]
REFL = [b for b in S2_BANDS if b != "NDVI"]

S2C, S1C = "COPERNICUS/S2_SR_HARMONIZED", "COPERNICUS/S1_GRD"
CSP = "GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED"

# search windows, relative to the event's observed burn dates
PRE_W = (-110, -8)                        # pre-fire S2
POST_W = (10, 75)                         # post-fire S2
S1_PRE_BACK, S1_POST_FWD = 60, 75


def utm_epsg(lon, lat):
    return (32600 if lat >= 0 else 32700) + int((lon + 180) // 6) + 1


def iso(d, shift=0):
    return (dt.date.fromisoformat(str(d)[:10]) + dt.timedelta(days=shift)).isoformat()


def download(img, region, scale, epsg, path, tries=3, timeout=420):
    for i in range(tries):
        try:
            url = img.getDownloadURL({"region": region, "scale": scale, "crs": f"EPSG:{epsg}",
                                      "format": "GEO_TIFF", "filePerBand": False})
            with urllib.request.urlopen(url, timeout=timeout) as r:
                path.write_bytes(r.read())
            return True
        except Exception:
            if i == tries - 1:
                return False
    return False


# --------------------------------------------------------------------------- #
# imagery selection
# --------------------------------------------------------------------------- #
def _snow(im):
    """NDSI snow test. Snow flips NBR  and would otherwise read as severe burn"""
    return im.normalizedDifference(["B3", "B11"]).gt(0.42).And(im.select("B8").gt(1100))


def s2_mosaic(aoi, d0, d1):
    """Clearest-first Sentinel-2 mosaic over [d0, d1], cloud- and snow-masked.

    Returns (mosaic_all_bands, nominal_scene). The mosaic keeps every band so the same
    image serves both the NBR used for the label and the 12-band export stack.
    """
    coll = (ee.ImageCollection(S2C).filterBounds(aoi).filterDate(d0, d1)
            .sort("CLOUDY_PIXEL_PERCENTAGE").limit(25)
            .linkCollection(ee.ImageCollection(CSP), ["cs"]))

    def score(im):
        clear = im.select("cs").gte(CS_THRESH).And(_snow(im).Not())
        frac = clear.reduceRegion(ee.Reducer.mean(), aoi, 120,
                                  maxPixels=1e8, bestEffort=True).get("cs")
        return im.updateMask(clear).set("clear_frac", frac)

    scored = coll.map(score).filter(ee.Filter.notNull(["clear_frac"])).sort("clear_frac")
    nominal = ee.Image(scored.sort("clear_frac", False).first())
    mosaic = scored.mosaic()               # ascending sort -> clearest last -> on top
    return mosaic, nominal.set("n_scenes", scored.size())


def s2_export(mosaic):
    """The 12-band flood stack, float32, in flood's fixed order."""
    ndvi = mosaic.normalizedDifference(["B8A", "B4"]).rename("NDVI")
    return mosaic.select(REFL).addBands(ndvi).select(S2_BANDS).toFloat()


def s1_mosaic(aoi, target, back, fwd):
    """Nearest-first Sentinel-1 GRD (IW, VV+VH) mosaic around `target`.

    Returns (mosaic, nominal_scene, n_scenes). n_scenes is 0 where the AOI has no SAR
    acquisition at all — a coverage gap, not a failure.
    """
    t = ee.Date(target)
    coll = (ee.ImageCollection(S1C).filterBounds(aoi)
            .filterDate(t.advance(-back, "day"), t.advance(fwd, "day"))
            .filter(ee.Filter.eq("instrumentMode", "IW"))
            .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VV"))
            .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VH")))
    coll = coll.map(lambda im: im.set("dt", ee.Number(im.date().difference(t, "day")).abs()))
    nominal = ee.Image(coll.sort("dt").first())
    mosaic = coll.sort("dt", False).mosaic()      # descending dt -> nearest on top
    return mosaic, nominal, coll.size()


def event_layers(ev):
    """Every per-event Earth Engine image, built once and reused for all its chips."""
    lon = (float(ev["bbox_w"]) + float(ev["bbox_e"])) / 2
    lat = (float(ev["bbox_s"]) + float(ev["bbox_n"])) / 2
    aoi = ee.Geometry.Rectangle([float(ev["bbox_w"]), float(ev["bbox_s"]),
                                 float(ev["bbox_e"]), float(ev["bbox_n"])], None, False)
    ig, fd = ev["event_date_start"], ev["event_date_end"]
    pre_w = (iso(ig, PRE_W[0]), iso(ig, PRE_W[1]))
    post_w = (iso(fd, POST_W[0]), iso(fd, POST_W[1]))

    s2pre, s2pre_nom = s2_mosaic(aoi, *pre_w)
    s2post, s2post_nom = s2_mosaic(aoi, *post_w)

    nbr_pre = s2pre.normalizedDifference(["B8", "B12"])
    nbr_post = s2post.normalizedDifference(["B8", "B12"])
    dnbr = nbr_pre.subtract(nbr_post).multiply(1000).rename("dnbr")
    sev = (ee.Image(0).where(dnbr.gte(SEV_T["low"]), 1).where(dnbr.gte(SEV_T["modlow"]), 2)
           .where(dnbr.gte(SEV_T["modhigh"]), 3).where(dnbr.gte(SEV_T["high"]), 4)
           .rename("sev").toInt())
    obs = nbr_pre.mask().And(nbr_post.mask()).rename("obs").toInt()

    # the MCD64A1 patch that defined this event — bounds where a burn may be claimed
    perim = (ee.ImageCollection("MODIS/061/MCD64A1").filterDate(iso(ig, -20), iso(fd, 20))
             .select("BurnDate").max().gt(0).unmask(0).focalMax(1).rename("inperim").toInt())
    crop = ee.ImageCollection("ESA/WorldCover/v200").first().select("Map").eq(40).rename("crop")

    s1pre, s1pre_nom, s1pre_n = s1_mosaic(aoi, iso(ig, -5), S1_PRE_BACK, 0)
    s1post, s1post_nom, s1post_n = s1_mosaic(aoi, iso(fd, 8), 0, S1_POST_FWD)

    return dict(
        epsg=utm_epsg(lon, lat), aoi=aoi, pre_w=pre_w, post_w=post_w,
        s2_pre=s2_export(s2pre), s2_post=s2_export(s2post),
        s2_pre_nom=s2pre_nom, s2_post_nom=s2post_nom,
        s1_pre=s1pre.select(["VV", "VH"]).toFloat(),
        s1_post=s1post.select(["VV", "VH"]).toFloat(),
        s1_pre_nom=s1pre_nom, s1_post_nom=s1post_nom,
        s1_pre_n=s1pre_n, s1_post_n=s1post_n,
        label_stack=dnbr.addBands(sev).addBands(crop).addBands(obs).addBands(perim).toFloat(),
        scout=obs.And(sev.gte(1)).And(crop).And(perim).unmask(0).rename("bc").toFloat())


def pick_centres(path, n):
    """Non-overlapping windows with the most burned cropland, best first.

    A box filter at chip size turns 'burned cropland density' into a per-pixel score;
    taking the argmax and then suppressing its neighbourhood gives disjoint chips.
    """
    from scipy.ndimage import uniform_filter
    with rasterio.open(path) as s:
        a = np.nan_to_num(s.read(1)); tr = s.transform; res = abs(tr.a)
    win = max(int(round(SIZE * RES / res)), 3)
    d = uniform_filter(a.astype("float32"), size=win, mode="constant")
    H, W = a.shape; m = win // 2
    d[:m, :] = d[-m:, :] = d[:, :m] = d[:, -m:] = -1      # no chip may run off the AOI
    out = []
    for _ in range(n):
        r, c = np.unravel_index(int(np.argmax(d)), d.shape)
        if d[r, c] <= 0:
            break
        out.append((*rxy(tr, r, c), float(d[r, c]), int(r), int(c)))
        r0, r1 = max(0, r - win), min(H, r + win)
        c0, c1 = max(0, c - win), min(W, c + win)
        d[r0:r1, c0:c1] = -1                              # suppress to disjoint chips
    return out


def scene_metadata(L, s1_pre_ok, s1_post_ok):
    """One round trip for every scene property the sidecars need."""
    meta = ee.Dictionary({
        "s2pre_date": L["s2_pre_nom"].date().format("YYYY-MM-dd"),
        "s2pre_clear": L["s2_pre_nom"].get("clear_frac"),
        "s2post_date": L["s2_post_nom"].date().format("YYYY-MM-dd"),
        "s2post_clear": L["s2_post_nom"].get("clear_frac"),
        "s2pre_n": L["s2_pre_nom"].get("n_scenes"),
        "s2post_n": L["s2_post_nom"].get("n_scenes"),
    }).getInfo()
    if s1_pre_ok:
        meta.update(ee.Dictionary({
            "s1pre_date": L["s1_pre_nom"].date().format("YYYY-MM-dd"),
            "s1pre_pass": L["s1_pre_nom"].get("orbitProperties_pass"),
            "s1pre_orbit": L["s1_pre_nom"].get("relativeOrbitNumber_start")}).getInfo())
    if s1_post_ok:
        meta.update(ee.Dictionary({
            "s1post_date": L["s1_post_nom"].date().format("YYYY-MM-dd"),
            "s1post_pass": L["s1_post_nom"].get("orbitProperties_pass"),
            "s1post_orbit": L["s1_post_nom"].get("relativeOrbitNumber_start")}).getInfo())
    return meta


def write_label(path, arr):
    """Collapse the raw stack into the 5-band uint8 label and overwrite in place.

    Returns the derived masks and counts the sidecar needs.
    """
    dnbr, sev, crop, obs, inperim = [np.nan_to_num(b) for b in arr]
    observed = obs >= 0.5
    cropm = crop >= 0.5
    burned = observed & (sev >= 1) & (inperim >= 0.5)
    excluded = cropm & ~observed              # cropland where cloud made dNBR unusable

    label = np.full(observed.shape, 4, np.uint8)          # 4 = non-cropland
    label[cropm] = 2                                      # 2 = unburned cropland
    label[cropm & burned] = 1                             # 1 = burned cropland
    label[cropm & excluded] = 3                           # 3 = excluded cropland

    sev_u8 = np.clip(sev, 0, 4).astype(np.uint8); sev_u8[~observed] = 255
    dn_u8 = np.clip((dnbr + 500) / 8, 0, 254).astype(np.uint8); dn_u8[~observed] = 255
    stack = np.stack([label, dn_u8, excluded.astype(np.uint8),
                      cropm.astype(np.uint8), sev_u8])

    with rasterio.open(path) as s:
        prof = s.profile
    prof.update(count=5, dtype="uint8", nodata=255, compress="deflate",
                height=SIZE, width=SIZE, photometric="MINISBLACK")
    with rasterio.open(path, "w", **prof) as d:
        d.write(stack)
    return label, sev_u8, observed, cropm, excluded


def build_event(ev):
    eid = ev["event_id"]
    odir = CHIPS / eid
    odir.mkdir(parents=True, exist_ok=True)

    # `_done` marks an event whose sampling already ran, including zero-yield ones, so a
    # re-run over a larger event pool skips it without re-downloading its scout raster.
    marker = odir / "_done"
    if marker.exists():
        return len(list(odir.glob("*.json"))), "cached"

    L = event_layers(ev)
    proj = ee.Projection(f"EPSG:{L['epsg']}")

    scout = odir / "_scout.tif"
    if not download(L["scout"], L["aoi"], 100, L["epsg"], scout):
        return 0, "scout download failed"
    centres = pick_centres(scout, CHIPS_PER_EVENT)
    scout.unlink(missing_ok=True)
    if not centres:
        marker.write_text("0")
        return 0, "no burned-cropland candidates"

    try:
        have = ee.Dictionary({"pre": L["s1_pre_n"], "post": L["s1_post_n"]}).getInfo()
        s1_pre_ok, s1_post_ok = int(have["pre"]) > 0, int(have["post"]) > 0
        meta = scene_metadata(L, s1_pre_ok, s1_post_ok)
    except Exception as e:
        return 0, f"scene metadata failed: {type(e).__name__}"

    made = 0
    for k, (x, y, _score, rr, cc) in enumerate(centres):
        chip_id = f"{eid}_r{rr:04d}_c{cc:04d}"
        side = odir / f"{chip_id}.json"
        if side.exists():
            made += 1
            continue

        rect = ee.Geometry.Rectangle([x - HALF, y - HALF, x + HALF, y + HALF], proj, False, True)
        lab_p = odir / f"{chip_id}_label.tif"
        if not download(L["label_stack"], rect, RES, L["epsg"], lab_p):
            continue
        with rasterio.open(lab_p) as s:
            arr = s.read()[:, :SIZE, :SIZE]

        label, sev_u8, observed, cropm, excluded = write_label(lab_p, arr)

        crop_px = int(cropm.sum())
        burn_px = int((label == 1).sum())
        crop_frac = crop_px / label.size
        burn_frac_of_crop = burn_px / max(crop_px, 1)
        if crop_frac < MIN_CROP_FRAC or burn_frac_of_crop < MIN_BURN_FRAC:
            lab_p.unlink(missing_ok=True)
            continue

        ok = {}
        jobs = [("s2_pre", L["s2_pre"]), ("s2_post", L["s2_post"])]
        if s1_pre_ok:
            jobs.append(("s1_pre", L["s1_pre"]))
        if s1_post_ok:
            jobs.append(("s1_post", L["s1_post"]))
        for nm, img in jobs:
            ok[nm] = download(img, rect, RES, L["epsg"], odir / f"{chip_id}_{nm}.tif")

        with rasterio.open(lab_p) as s:
            bounds, crs = s.bounds, s.crs
        w_, s_, e_, n_ = rwarp.transform_bounds(crs, "EPSG:4326", *bounds)
        sv = sev_u8[(label == 1) & (sev_u8 != 255)]

        # Paths are relative to the chips root so a sidecar means the same thing on any
        # machine. NOTE ON DENOMINATORS: `cropland_frac` and `burned_crop_frac` are
        # shares of the WHOLE CHIP, while `excluded_crop_frac` is a share OF CROPLAND.
    
        rel = lambda name: f"{eid}/{chip_id}_{name}.tif"
        side.write_text(json.dumps({
            "chip_id": chip_id, "event_id": eid, "rank": k + 1,
            "source": ev["source"], "tier": ev["tier"], "continent": ev["continent"],
            "country": ev["country"], "place": ev["place"],
            "provenance": {
                "origin": "MCD64A1 burned-area clustering (GlobFire method)",
                "confirmed_by": ["MODIS MCD64A1", "VIIRS VNP64A1", "VIIRS VNP14A1 active fire"],
                "n_independent": int(ev["n_indep"]), "tier": ev["tier"],
                "burn_km2": float(ev["burn_km2"]),
                "modis_frac": float(ev["modis_frac"]), "viirs_frac": float(ev["viirs_frac"]),
                "active_fire_px": int(ev["af_count"])},
            "event_date_start": ev["event_date_start"], "event_date_end": ev["event_date_end"],
            "crs": f"EPSG:{L['epsg']}", "res_m": RES, "size": SIZE,
            "label_path": rel("label"),
            "burn_severity_mean": round(float(sv.mean()), 2) if sv.size else 0.0,
            "cropland_frac": round(crop_frac, 4),
            "burned_crop_frac": round(burn_px / label.size, 4),
            "burned_crop_px": burn_px,
            "excluded_crop_frac": round(float(excluded.sum()) / max(crop_px, 1), 4),
            "min_lon": round(w_, 6), "min_lat": round(s_, 6),
            "max_lon": round(e_, 6), "max_lat": round(n_, 6),
            "imagery": {
                "s2_pre": {"date": meta["s2pre_date"],
                           "clear_frac": round(meta["s2pre_clear"] or 0, 4),
                           "n_scenes": meta["s2pre_n"],
                           "path": rel("s2_pre"), "bands": S2_BANDS},
                "s2_post": {"date": meta["s2post_date"],
                            "clear_frac": round(meta["s2post_clear"] or 0, 4),
                            "n_scenes": meta["s2post_n"],
                            "path": rel("s2_post"), "bands": S2_BANDS},
                "s1_pre": ({"date": meta["s1pre_date"], "path": rel("s1_pre"),
                            "pass": meta["s1pre_pass"], "rel_orbit": meta["s1pre_orbit"],
                            "bands": ["VV", "VH"]} if s1_pre_ok else None),
                "s1_post": ({"date": meta["s1post_date"], "path": rel("s1_post"),
                             "pass": meta["s1post_pass"], "rel_orbit": meta["s1post_orbit"],
                             "bands": ["VV", "VH"]} if s1_post_ok else None),
                "s1_note": (None if (s1_pre_ok and s1_post_ok) else
                            "no Sentinel-1 IW VV+VH acquisition over this AOI in the "
                            "pre and/or post window (genuine SAR coverage gap)"),
                "s2_pre_window": list(L["pre_w"]), "s2_post_window": list(L["post_w"])},
            "qc": {"all_rasters_ok": all(ok.values()),
                   "rasters": {k2: ("ok" if v else "missing") for k2, v in ok.items()},
                   "observed_frac": round(float(observed.mean()), 4),
                   "verdict": "clean" if all(ok.values()) and observed.mean() > 0.5 else "check",
                   "source": "chip QC",
                   "method": "download success per raster + S2 observability "
                             "(clear pre AND post) fraction"},
        }, indent=2))
        made += 1

    marker.write_text(str(made))
    return made, "ok"


def main(argv=None):
    ap = argparse.ArgumentParser(description="Stage B — sample and download chips")
    ap.add_argument("--limit", type=int, default=None, help="only build the first N events")
    ap.add_argument("--all", action="store_true", help="build every event in the layer")
    ap.add_argument("--start", type=int, default=0, help="skip the first N selected events")
    ap.add_argument("--cont-cap", type=int, default=CONT_CAP,
                    help="max events per continent — limits the Africa/South Asia skew")
    a = ap.parse_args(argv)

    events_csv = MAN / "events.csv"
    if not events_csv.exists():
        sys.exit(f"no event layer at {events_csv} — run 01_build_events.py first")
    rows = list(csv.DictReader(open(events_csv)))

    seen, sel = {}, []
    for r in rows:
        c = r["continent"]
        if seen.get(c, 0) >= a.cont_cap:
            continue
        seen[c] = seen.get(c, 0) + 1
        sel.append(r)
    sel = sel[a.start:]
    if a.limit:
        sel = sel[:a.limit]

    CHIPS.mkdir(parents=True, exist_ok=True)
    print(f"chipping {len(sel)} events (continent cap {a.cont_cap}) -> {CHIPS}")
    total = 0
    for i, ev in enumerate(sel, 1):
        try:
            n, msg = build_event(ev)
        except Exception as e:
            n, msg = 0, f"{type(e).__name__}: {str(e)[:70]}"
        total += n
        print(f"  [{i}/{len(sel)}] {ev['event_id']} {ev['continent'][:12]:12s} "
              f"{str(ev['country'])[:12]:12s} -> {n} chips ({msg}) | total {total}")
    print(f"\nDONE: {total} chips from {len(sel)} events")
    return 0


if __name__ == "__main__":
    sys.exit(main())
