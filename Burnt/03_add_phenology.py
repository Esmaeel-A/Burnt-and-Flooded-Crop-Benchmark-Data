#!/usr/bin/env python3
"""
Stage C1: add a phenology information to every chip sidecar.

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
import os
import argparse, datetime as dt, json, sys, time
from collections import defaultdict
from pathlib import Path

import ee

PROJ = Path.cwd()
CHIPS = PROJ / "dataset" / "processed" / "chips"

COLL = "NASA/VIIRS/002/VNP22Q2"
EPOCH = dt.date(2000, 1, 1)      # VNP22Q2 dates are days since 2000-01-01, NOT the Unix epoch
FILL = 32000
PHEN_MAX_YEAR = 2024             # coverage ends here; later burns -> null + note

STAGE = {0: "dormant", 1: "green-up", 2: "peak/maturity", 3: "senescence"}
DATE_BANDS = ["Onset_Greenness_Increase", "Onset_Greenness_Maximum",
              "Onset_Greenness_Decrease", "Onset_Greenness_Minimum"]
SHORT = ["green_up_onset", "maximum", "senescence_onset", "dormancy_onset"]
SRC = "VIIRS VNP22Q2.002 land-surface phenology (500 m, yearly)"
METHOD = ("per-pixel stage from that pixel's own transition dates across 4 candidate "
          "cycles (granules Y-1,Y x cycles 1,2); majority vote over the chip")


def days_since_epoch(s):
    return (dt.date.fromisoformat(str(s)[:10]) - EPOCH).days


def epoch_to_iso(v):
    return (EPOCH + dt.timedelta(days=int(round(v)))).isoformat()


def stage_images(year, burn_k):
    """Stage / winning-cycle images plus the per-cycle date bands.

    Walks the four candidate cycles in order and lets the FIRST coherent one that
    brackets the burn date claim each pixel, so a pixel is never scored against a cycle
    that does not contain its burn.
    """
    anchor = None
    parts = []
    stage = has = anyok = cyc = None
    for gi, gy in enumerate((year - 1, year)):
        img = ee.ImageCollection(COLL).filterDate(f"{gy}-01-01", f"{gy+1}-01-01").first()
        if anchor is None:
            anchor = img.select(f"{DATE_BANDS[0]}_1").toFloat().multiply(0).unmask(0)
            stage = has = anyok = cyc = anchor
        for c in (1, 2):
            b = [img.select(f"{x}_{c}").toFloat() for x in DATE_BANDS]
            ok = (b[0].lt(FILL).And(b[1].lt(FILL)).And(b[2].lt(FILL)).And(b[3].lt(FILL))
                  .And(b[0].lt(b[1])).And(b[1].lt(b[2])).And(b[2].lt(b[3]))).unmask(0)
            contains = ok.And(b[0].lte(burn_k)).And(b[3].gte(burn_k)).unmask(0)
            s = (anchor.where(contains.And(b[1].gt(burn_k)), 1)
                       .where(contains.And(b[1].lte(burn_k)).And(b[2].gt(burn_k)), 2)
                       .where(contains.And(b[2].lte(burn_k)), 3))
            fresh = contains.And(has.Not())          # first coherent cycle wins the pixel
            stage = stage.where(fresh, s)
            cyc = cyc.where(fresh, gi * 10 + c)
            has = has.Or(contains)
            anyok = anyok.Or(ok)
            n = gi * 10 + c
            for nm, bb in zip(("inc", "max", "dec", "min"), b):
                parts.append(bb.updateMask(ok.selfMask()).rename(f"c{n}_{nm}"))
    return (stage.updateMask(anyok.gt(0)).rename("stage"),
            cyc.updateMask(anyok.gt(0)).rename("cyc"),
            ee.Image.cat(parts))


def collect(force, limit):
    """Sidecars needing phenology, grouped by event — chips of an event share date+year."""
    by_ev = defaultdict(list)
    skipped = 0
    for jf in sorted(CHIPS.glob("*/*.json")):
        try:
            d = json.loads(jf.read_text())
        except Exception:
            continue
        if not force and "phenology" in d:
            skipped += 1
            continue
        if not d.get("event_date_start") or d.get("min_lon") is None:
            continue
        by_ev[d["event_id"]].append((jf, d))
    if limit:
        by_ev = dict(list(by_ev.items())[:limit])
    return by_ev, skipped


def label_event(items, year, burn_k):
    """Return {chip_id: phenology_dict_or_None}. Retries around GEE rate limits."""
    feats = [ee.Feature(ee.Geometry.Rectangle([float(d["min_lon"]), float(d["min_lat"]),
                                               float(d["max_lon"]), float(d["max_lat"])]),
                        {"cid": d["chip_id"]}) for _, d in items]
    fc = ee.FeatureCollection(feats)
    hist = meds = None
    for attempt in (1, 2, 3):
        try:
            st, cy, dates = stage_images(year, burn_k)
            hist = st.addBands(cy).reduceRegions(fc, ee.Reducer.frequencyHistogram(), 500).getInfo()
            meds = dates.reduceRegions(fc, ee.Reducer.median(), 500).getInfo()
            break
        except Exception:
            if attempt == 3:
                return None
            time.sleep(2.0 * attempt)
    if hist is None:
        return None
    time.sleep(0.25)                                  # pace to stay under the rate limit

    medmap = {f["properties"].get("cid"): f["properties"] for f in meds["features"]}
    out = {}
    for f in hist["features"]:
        p = f["properties"]
        cid = p.get("cid")
        counts = {int(float(k)): v for k, v in (p.get("stage") or {}).items()}
        if not counts:
            out[cid] = None
            continue
        n = sum(counts.values())
        dom = max(counts, key=counts.get)
        cyc_counts = {int(float(k)): v for k, v in (p.get("cyc") or {}).items() if float(k)}
        dc = max(cyc_counts, key=cyc_counts.get) if cyc_counts else None
        mp = medmap.get(cid, {})
        transitions = {s: (epoch_to_iso(mp[f"c{dc}_{k}"])
                           if dc and mp.get(f"c{dc}_{k}") is not None else None)
                       for s, k in zip(SHORT, ("inc", "max", "dec", "min"))}
        out[cid] = {
            "stage_at_burn": STAGE[dom],
            "stage_agreement": round(counts[dom] / n, 3),
            "stage_fractions": {STAGE[k]: round(v / n, 3) for k, v in sorted(counts.items())},
            "growth_cycle": (None if not dc else
                             f"{'previous' if dc < 10 else 'event'}-year cycle {dc % 10}"),
            "transitions": transitions,
            "n_valid_pixels": int(round(n)),
            "source": SRC,
            "method": METHOD,
        }
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="Stage C1 — add phenology to chip sidecars")
    ap.add_argument("--force", action="store_true", help="relabel chips that already have it")
    ap.add_argument("--limit", type=int, default=0, help="only the first N events")
    a = ap.parse_args(argv)

    if not os.environ.get("EE_PROJECT"):
        sys.exit("Set EE_PROJECT in the environment (export EE_PROJECT=your-project-id).")
    ee.Initialize(project=os.environ["EE_PROJECT"])

    by_ev, skipped = collect(a.force, a.limit)
    total = sum(len(v) for v in by_ev.values())
    if not total:
        print(f"nothing to do ({skipped} chips already labelled)")
        return 0
    print(f"{total} chips to label across {len(by_ev)} events ({skipped} already done)", flush=True)

    t0 = time.time()
    done = nulls = errs = 0
    for i, (eid, items) in enumerate(sorted(by_ev.items()), 1):
        d0 = items[0][1]
        year = int(d0["event_date_start"][:4])
        burn_k = days_since_epoch(d0["event_date_start"])

        if year > PHEN_MAX_YEAR:
            for jf, d in items:
                d["phenology"] = None
                d["phenology_note"] = ("VNP22Q2 land-surface phenology covers 2013-2024; "
                                       f"no data for event year {year}")
                jf.write_text(json.dumps(d, indent=2))
                nulls += 1
            continue

        res = label_event(items, year, burn_k)
        if res is None:
            errs += len(items)
            print(f"  [{i}/{len(by_ev)}] {eid}: GEE error after 3 attempts", flush=True)
            continue

        for jf, d in items:
            ph = res.get(d["chip_id"], None)
            if ph is None:
                d["phenology"] = None
                d["phenology_note"] = "no fill-free, monotonic VNP22Q2 cycle over this chip"
                nulls += 1
            else:
                ph["burn_date"] = d["event_date_start"][:10]
                d["phenology"] = ph
                d.pop("phenology_note", None)
                done += 1
            jf.write_text(json.dumps(d, indent=2))    # in place -> hardlinks stay shared

        if i % 25 == 0:
            el = time.time() - t0
            print(f"  [{i}/{len(by_ev)} events] labelled={done} null={nulls} err={errs} "
                  f"| {el/60:.1f} min elapsed, ~{el/i*(len(by_ev)-i)/60:.0f} min left", flush=True)

    print(f"\nDONE: {done} chips labelled, {nulls} null (no coherent cycle or out of coverage), "
          f"{errs} errors in {(time.time()-t0)/60:.1f} min")
    return 0


if __name__ == "__main__":
    sys.exit(main())
