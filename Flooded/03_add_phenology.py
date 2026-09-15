#!/usr/bin/env python3
"""
Stage C1: add phenology information to every chip sidecar.


    VNP22Q2 runs 2013-2024. Chips from later years get `phenology: null` plus a
    `phenology_note` saying why 

    For each pixel, read its OWN transition dates from VIIRS VNP22Q2 across four
    candidate growth cycles (granules Y-1 and Y, cycles 1 and 2). Keep a cycle only if
    it is fill-free and monotonic — greenness increase < maximum < decrease < minimum.
    Locate the burn date within the first surviving cycle to get that pixel's stage,
    then take a majority vote over the chip.

      dormant        burn falls outside every coherent cycle
      green-up       between greenness increase and maximum
      peak/maturity  between maximum and decrease
      senescence     between decrease and minimum

"""
from __future__ import annotations
import argparse, datetime as dt, json, os, sys, time
from collections import defaultdict
from pathlib import Path

REPO = Path.cwd()
ROOT = REPO / "release" / "chips"
COLL, EPOCH, FILL = "NASA/VIIRS/002/VNP22Q2", dt.date(2000, 1, 1), 32000
PHEN_MAX_YEAR = 2024             # coverage ends here; later floods -> null + note
STAGE = {0: "dormant", 1: "green-up", 2: "peak/maturity", 3: "senescence"}
DB = ["Onset_Greenness_Increase", "Onset_Greenness_Maximum",
      "Onset_Greenness_Decrease", "Onset_Greenness_Minimum"]
SHORT = ["green_up_onset", "maximum", "senescence_onset", "dormancy_onset"]


def init_ee():
    """Initialise Earth Engine with the Cloud project in EE_PROJECT."""
    import ee
    project = os.environ.get("EE_PROJECT", "").strip()
    if not project:
        raise RuntimeError(
            "Set EE_PROJECT in the environment (export EE_PROJECT=your-project-id).")
    try:
        ee.Initialize(project=project)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(
            f"ee.Initialize(project={project!r}) failed: {exc}. Confirm the project is registered "
            "for Earth Engine and that `earthengine authenticate` has been run.") from exc
    return ee


def d2k(s):
    return (dt.date.fromisoformat(s[:10]) - EPOCH).days


def k2d(v):
    return (EPOCH + dt.timedelta(days=int(round(v)))).isoformat()


def build(ee, yr, F):
    """stage/cycle images + per-cycle coherence-masked date bands."""
    anchor = None; parts = []
    stage = has = anyok = cid = None
    for gi, gy in enumerate((yr - 1, yr)):
        img = ee.ImageCollection(COLL).filterDate(f"{gy}-01-01", f"{gy+1}-01-01").first()
        b = [img.select(f"{x}_1").toFloat() for x in DB]
        if anchor is None:
            anchor = b[0].multiply(0).unmask(0)
            stage = has = anyok = cid = anchor
        for c in (1, 2):
            b = [img.select(f"{x}_{c}").toFloat() for x in DB]
            ok = (b[0].lt(FILL).And(b[1].lt(FILL)).And(b[2].lt(FILL)).And(b[3].lt(FILL))
                  .And(b[0].lt(b[1])).And(b[1].lt(b[2])).And(b[2].lt(b[3]))).unmask(0)
            cont = ok.And(b[0].lte(F)).And(b[3].gte(F)).unmask(0)
            s = (anchor.where(cont.And(b[1].gt(F)), 1)
                       .where(cont.And(b[1].lte(F)).And(b[2].gt(F)), 2)
                       .where(cont.And(b[2].lte(F)), 3))
            fresh = cont.And(has.Not())
            stage = stage.where(fresh, s); cid = cid.where(fresh, gi * 10 + c)
            has = has.Or(cont); anyok = anyok.Or(ok)
            n = gi * 10 + c
            for nm, bb in zip(("inc", "max", "dec", "min"), b):
                parts.append(bb.updateMask(ok.selfMask()).rename(f"c{n}_{nm}"))
    return (stage.updateMask(anyok.gt(0)).rename("stage"),
            cid.updateMask(anyok.gt(0)).rename("cyc"), ee.Image.cat(parts))


def collect(force, limit):
    """Sidecars needing phenology, grouped by event (chips of an event share flood date + year)."""
    by_ev = defaultdict(list)
    skipped = 0
    for jf in ROOT.glob("*/*.json"):
        try:
            d = json.loads(jf.read_text())
        except Exception:
            continue
        if not force and "phenology" in d:
            skipped += 1; continue
        if not d.get("event_date_start") or d.get("min_lon") is None:
            continue
        by_ev[d["event_id"]].append((jf, d))
    if limit:
        by_ev = dict(sorted(by_ev.items())[:limit])
    return by_ev, skipped


def label_event(ee, items, yr, F):
    """(hist, meds) for one event's chips, or (None, exc) after 3 failed attempts."""
    feats = [ee.Feature(ee.Geometry.Rectangle([float(d["min_lon"]), float(d["min_lat"]),
                                               float(d["max_lon"]), float(d["max_lat"])]),
                        {"cid": d["chip_id"]}) for _, d in items]
    fc = ee.FeatureCollection(feats)
    for attempt in (1, 2, 3):                      # GEE rate-limit tolerance
        try:
            st, cy, dates = build(ee, yr, F)
            hist = st.addBands(cy).reduceRegions(fc, ee.Reducer.frequencyHistogram(), 500).getInfo()
            meds = dates.reduceRegions(fc, ee.Reducer.median(), 500).getInfo()
            time.sleep(0.25)                       # pace to stay under the rate limit
            return hist, meds
        except Exception as exc:
            if attempt == 3:
                return None, exc
            time.sleep(2.0 * attempt)


def write_phenology(items, hist, meds):
    """Write the `phenology` block into each chip GEE returned. Returns (labelled, null)."""
    done = nulls = 0
    medmap = {f["properties"].get("cid"): f["properties"] for f in meds["features"]}
    for f in hist["features"]:
        p = f["properties"]; cidn = p.get("cid")
        jf, d = next(((a, b) for a, b in items if b["chip_id"] == cidn), (None, None))
        if jf is None:
            continue
        sh = {int(float(k)): v for k, v in (p.get("stage") or {}).items()}
        if not sh:
            d["phenology"] = None
            d["phenology_note"] = "no fill-free, monotonic VNP22Q2 cycle over this chip"
            nulls += 1
        else:
            n = sum(sh.values()); dom = max(sh, key=sh.get)
            ch = {int(float(k)): v for k, v in (p.get("cyc") or {}).items() if float(k)}
            dc = max(ch, key=ch.get) if ch else None
            mp = medmap.get(cidn, {})
            tr = {s: (k2d(mp[f"c{dc}_{k}"]) if dc and mp.get(f"c{dc}_{k}") is not None else None)
                  for s, k in zip(SHORT, ("inc", "max", "dec", "min"))}
            d["phenology"] = {
                "stage_at_flood": STAGE[dom],
                "stage_agreement": round(sh[dom] / n, 3),
                "stage_fractions": {STAGE[k]: round(v / n, 3) for k, v in sorted(sh.items())},
                "flood_date": d["event_date_start"][:10],
                "growth_cycle": (None if not dc else
                                 f"{'previous' if dc < 10 else 'event'}-year cycle {dc % 10}"),
                "transitions": tr,
                "n_valid_pixels": int(round(n)),
            }
            d.pop("phenology_note", None)
            done += 1
        jf.write_text(json.dumps(d, indent=2))       # in place -> hardlink stays shared
    return done, nulls


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Stage C1 — add phenology to chip sidecars")
    ap.add_argument("--force", action="store_true", help="relabel chips that already have it")
    ap.add_argument("--limit", type=int, default=0, help="only the first N events")
    a = ap.parse_args(argv)

    ee = init_ee()
    by_ev, skipped = collect(a.force, a.limit)
    print(f"{sum(len(v) for v in by_ev.values())} chips to label across {len(by_ev)} events "
          f"({skipped} already done)", flush=True)

    t0 = time.time(); done = nulls = errs = 0
    for i, (eid, items) in enumerate(sorted(by_ev.items()), 1):
        d0 = items[0][1]
        yr = int(d0["event_date_start"][:4]); F = d2k(d0["event_date_start"])

        if yr > PHEN_MAX_YEAR:                       # outside VNP22Q2 coverage: null, not a proxy
            for jf, d in items:
                d["phenology"] = None
                d["phenology_note"] = ("VNP22Q2 land-surface phenology covers 2013-2024; "
                                       f"no data for event year {yr}")
                jf.write_text(json.dumps(d, indent=2))
                nulls += 1
            continue

        hist, meds = label_event(ee, items, yr, F)
        if hist is None:
            errs += len(items)
            print(f"  [{i}/{len(by_ev)}] {eid}: GEE error {type(meds).__name__}", flush=True)
            continue
        dn, nl = write_phenology(items, hist, meds)
        done += dn; nulls += nl
        if i % 25 == 0:
            el = time.time() - t0
            print(f"  [{i}/{len(by_ev)} events] labelled={done} null={nulls} err={errs} "
                  f"| {el/60:.1f} min elapsed, ~{el/i*(len(by_ev)-i)/60:.0f} min left", flush=True)

    print(f"\nDONE: {done} chips labelled, {nulls} null (no coherent cycle or out of coverage), "
          f"{errs} errors in {(time.time()-t0)/60:.1f} min")
    return 0


if __name__ == "__main__":
    sys.exit(main())
