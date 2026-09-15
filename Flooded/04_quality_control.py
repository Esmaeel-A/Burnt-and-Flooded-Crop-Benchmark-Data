#!/usr/bin/env python3
"""
Stage C2 : quality control over every chip raster.


    scan      every raster -> one primary status + independent defect flags
    roll up   4 raster statuses -> one chip verdict
    patch     write a `qc` block into every <chip_id>.json
    report    data/qc_chips.csv, data/qc_by_chip.csv, stats/QC_STATS.md, 2 figures

 Criteria

    missing      file absent
    unreadable   open/read raises, wrong shape or band count, or a constant raster
    empty        >= 98% nodata (S1 fully off-swath; an S2 tripwire that never fires)
    corrupt      S2 median visible reflectance > 10000 (>100%, physically impossible),
                 or S1 single-polarisation dropout
    partial      nodata wedge anchored to the border — a sensor footprint or orbit
                 swath edge cutting the chip; severity by fraction
    blank_white  S2 opaque cloud over >= 80% of the valid area
    speckled     scattered nodata dropouts, not edge-anchored — cosmetic
    cloudy       S2 partial cloud, 35-80% of valid area — advisory, still usable
    ok           none of the above

"""
from __future__ import annotations
import argparse, csv, json, os, sys, warnings
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import rasterio
from scipy import ndimage

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=RuntimeWarning)

REPO = Path.cwd()
ROOT = REPO / "release" / "chips"
REP = REPO / "release" / "report"
FIG = REP / "figures"

KEYS = ["s2_pre", "s2_post", "s1_pre", "s1_post"]
N_REFL = 11                         # B2..B12; index 11 (12th band) is NDVI, excluded
IB, IG, IR, ISWIR1 = 0, 1, 2, 9     # B2, B3, B4, B11


# ---- calibrated thresholds (full-corpus; see stats/qc_spec.json) ----
EMPTY_FRAC   = 0.98     # >= this nodata share -> no usable data
S2_PARTIAL   = 0.005    # S2 nodata >= this -> partial
S1_PARTIAL   = 0.02     # S1 nodata >= this -> partial (the empty valley of the S1 distribution)
CC_FRAC      = 0.95     # largest connected component >= this share of nodata -> single blob
PART_MOD     = 0.05     # minor < .05 <= moderate < .15 <= severe
PART_SEV     = 0.15
VMIN_WHITE   = 2000     # visible-band floor for a "white" pixel
FLATNESS     = 0.25     # (max-min)/mean over RGB below this = spectrally flat -> cloud, not sand
B11_WHITE    = 800      # cloud retains SWIR; snow insurance
W_WHITE      = 0.80     # >= this white share of valid px -> blank_white
W_CLOUD      = 0.35     # 0.35-0.80 -> cloudy (advisory)
SAT_VIS      = 10000    # median visible reflectance above this -> corrupt/saturated
SNOW_VETO_CF = 0.90     # Cloud Score+ clear_frac >= this -> never white
S1_CONST_STD = 1.0      # S1 valid-pixel dB std below this -> degenerate
# raster status -> usability class for the chip rollup
UNUSABLE = {"missing", "empty", "unreadable", "corrupt", "blank_white", "partial_severe"}
COMPROMISED = {"partial_moderate"}
# (ok, cloudy, speckled, partial_minor -> usable)

# ---- design tokens ----
SURFACE = "#fcfcfb"; INK = "#1a1a19"; INK_MUTED = "#55554f"; GRID = "#e4e4e0"
# severity ramp: clean/green -> advisory/amber -> hard-defect/red, missing = grey
STATUS_C = {
    "ok": "#3f9b52", "cloudy": "#e8a33d", "speckled": "#c9c9c4",
    "partial_minor": "#f2c14e", "partial_moderate": "#e07b39", "partial_severe": "#c0392b",
    "blank_white": "#7b8794", "corrupt": "#5a2a82", "empty": "#2c3e50", "missing": "#a0a0a0",
    "unreadable": "#000000",
}
STATUS_ORDER = ["ok", "cloudy", "speckled", "partial_minor", "partial_moderate",
                "partial_severe", "blank_white", "corrupt", "empty", "missing", "unreadable"]
plt.rcParams.update({"figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "text.color": INK, "axes.labelcolor": INK_MUTED, "axes.edgecolor": GRID,
                     "xtick.color": INK_MUTED, "ytick.color": INK_MUTED, "font.size": 10,
                     "axes.titlesize": 12, "axes.titleweight": "bold"})


def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(GRID); ax.spines["bottom"].set_color(GRID)
    ax.set_axisbelow(True)


# --------------------------------------------------------------------------- #
# scan
# --------------------------------------------------------------------------- #
def _edge_blob(mask: np.ndarray) -> tuple:
    """Return (largest_cc_frac, touches_border) for a boolean nodata mask."""
    tot = int(mask.sum())
    if tot == 0:
        return 0.0, False
    lab, n = ndimage.label(mask)                    # 4-connectivity
    if n == 0:
        return 0.0, False
    sizes = np.bincount(lab.ravel())[1:]
    k = int(sizes.argmax()) + 1
    cc = lab == k
    touches = bool(cc[0].any() or cc[-1].any() or cc[:, 0].any() or cc[:, -1].any())
    return float(sizes.max() / tot), touches


def _blank(f):
    return {"status": f, "nodata_frac": "", "edge_cc": "", "white_frac": "",
            "vis_med": "", "valid_std": "", "is_partial": 0, "is_white": 0, "is_corrupt": 0,
            "is_empty": 0, "is_speckle": 0, "is_cloudy": 0, "note": ""}


def _severity(frac):
    return "partial_minor" if frac < PART_MOD else ("partial_moderate" if frac < PART_SEV
                                                    else "partial_severe")


def _scan_raster(path: Path, kind: str, clear_frac) -> dict:
    r = _blank("ok")
    try:
        with rasterio.open(path) as ds:
            if (ds.width, ds.height) != (512, 512):
                r.update(status="unreadable", note=f"shape {ds.width}x{ds.height}"); return r
            want = 12 if kind == "s2" else 2
            if ds.count != want:
                r.update(status="unreadable", note=f"{ds.count} bands != {want}"); return r
            a = ds.read().astype("float32")
    except Exception as exc:
        r.update(status="unreadable", note=f"{type(exc).__name__}: {str(exc)[:60]}"); return r

    if not np.isfinite(a).all():
        r["note"] = f"nonfinite={float((~np.isfinite(a)).mean()):.4f}"
        a = np.nan_to_num(a, nan=0.0, posinf=0.0, neginf=0.0)

    core = a[:N_REFL] if kind == "s2" else a
    nod = (core == 0).all(axis=0)
    frac = float(nod.mean())
    valid = ~nod
    nvalid = int(valid.sum())
    r["nodata_frac"] = round(frac, 6)

    # --- degenerate / empty (highest precedence after unreadable) ---
    if frac >= EMPTY_FRAC:
        r.update(status="empty", is_empty=1); return r
    if nvalid >= 10000:
        vstd = float(core[:, valid].std())
        r["valid_std"] = round(vstd, 4)
        if (kind == "s1" and (float(a[0][valid].std()) < S1_CONST_STD
                              or float(a[1][valid].std()) < S1_CONST_STD)):
            r.update(status="unreadable", note="constant-valued S1"); return r

    # --- independent measurements: nodata geometry ---
    part_thr = S2_PARTIAL if kind == "s2" else S1_PARTIAL
    is_partial = is_speckle = False
    if frac >= part_thr:
        cc, touches = _edge_blob(nod)
        r["edge_cc"] = round(cc, 4)
        if cc >= CC_FRAC and touches:
            is_partial = True
        else:
            is_speckle = True
    elif frac > 0:
        is_speckle = True                            # tiny scattered dropout, ignored

    # --- S2 : whiteness + saturation ---
    is_white = is_corrupt = is_cloudy = False
    if kind == "s2":
        rgb = a[[IR, IG, IB]][:, valid]
        vis = rgb.mean(axis=0)
        vmed = float(np.median(vis)) if nvalid else 0.0
        r["vis_med"] = round(vmed, 1)
        vmin = rgb.min(axis=0); vmax = rgb.max(axis=0)
        flat = (vmax - vmin) / np.maximum(vis, 1.0)
        b11 = a[ISWIR1][valid]
        W = float(((vmin > VMIN_WHITE) & (flat < FLATNESS) & (b11 > B11_WHITE)).mean()) if nvalid else 0.0
        r["white_frac"] = round(W, 4)
        vetoed = (clear_frac is not None and clear_frac != "" and float(clear_frac) >= SNOW_VETO_CF)
        if vmed > SAT_VIS and not vetoed:
            is_corrupt = True
        elif W >= W_WHITE and not vetoed:
            is_white = True
        elif W >= W_CLOUD and not vetoed:
            is_cloudy = True
    else:
        # --- S1 : single-polarisation dropout ---
        fzv = float((a[0] == 0).mean()); fzh = float((a[1] == 0).mean())
        if abs(fzv - fzh) > 1e-4:
            is_corrupt = True; r["note"] = f"band dropout |dz|={abs(fzv-fzh):.4f}"

    r.update(is_partial=int(is_partial), is_white=int(is_white), is_corrupt=int(is_corrupt),
             is_speckle=int(is_speckle), is_cloudy=int(is_cloudy))

    # --- single primary status by precedence ---
    if is_corrupt:
        r["status"] = "corrupt"
    elif is_partial:
        r["status"] = _severity(frac)
    elif is_white:
        r["status"] = "blank_white"
    elif is_speckle:
        r["status"] = "speckled"
    elif is_cloudy:
        r["status"] = "cloudy"
    else:
        r["status"] = "ok"
    return r


def _scan_chip(args) -> list:
    evdir, chip_id, meta = args
    rows = []
    for k in KEYS:
        kind = "s2" if k.startswith("s2") else "s1"
        p = Path(evdir) / f"{chip_id}_{k}.tif"
        base = {"chip_id": chip_id, "event_id": meta["event_id"], "raster": k,
                "continent": meta.get("continent", ""), "country": meta.get("country", ""),
                "year": (meta.get("event_date_start") or "")[:4],
                "clear_frac": meta["clear"].get(k, "")}
        if p.exists():
            base.update(_scan_raster(p, kind, meta["clear"].get(k)))
        else:
            base.update(_blank("missing"))
        rows.append(base)
    return rows


def scan_all(workers, limit) -> int:
    """Scan every raster of every chip -> report/qc_chips.csv. Returns the chip count."""
    jobs = []
    for jf in sorted(ROOT.glob("*/*.json")):
        try:
            d = json.loads(jf.read_text())
        except Exception:
            continue
        im = d.get("imagery") or {}
        clear = {}
        for k in KEYS:
            e = im.get(k) or {}
            clear[k] = e.get("clear_frac") if e.get("clear_frac") is not None else ""
        jobs.append((str(jf.parent), d.get("chip_id", jf.stem),
                     {"event_id": d.get("event_id", jf.parent.name),
                      "continent": d.get("continent", ""), "country": d.get("country", ""),
                      "event_date_start": d.get("event_date_start", ""), "clear": clear}))
    if limit:
        jobs = jobs[:limit]
    if not jobs:
        return 0
    print(f"scanning {len(jobs)} chips x 4 rasters with {workers} workers", flush=True)

    rows, done = [], 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_scan_chip, j) for j in jobs]
        for f in as_completed(futs):
            rows.extend(f.result()); done += 1
            if done % 500 == 0:
                print(f"  {done}/{len(jobs)} chips", flush=True)

    REP.mkdir(parents=True, exist_ok=True)
    fields = ["chip_id", "event_id", "raster", "continent", "country", "year", "clear_frac",
              "status", "nodata_frac", "edge_cc", "white_frac", "vis_med", "valid_std",
              "is_partial", "is_white", "is_corrupt", "is_empty", "is_speckle", "is_cloudy", "note"]
    rows.sort(key=lambda r: (r["chip_id"], KEYS.index(r["raster"])))
    with open(REP / "qc_chips.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"  -> {REP/'qc_chips.csv'}  ({len(rows)} raster rows)")
    return len(jobs)


# --------------------------------------------------------------------------- #
# roll up
# --------------------------------------------------------------------------- #
def usability(st):
    if st in UNUSABLE:
        return "unusable"
    if st in COMPROMISED:
        return "compromised"
    return "usable"


def rollup():
    """report/qc_chips.csv -> per-chip verdicts -> report/qc_by_chip.csv."""
    rows = list(csv.DictReader(open(REP / "qc_chips.csv")))
    by_chip = defaultdict(dict)
    meta = {}
    for r in rows:
        by_chip[r["chip_id"]][r["raster"]] = r
        meta[r["chip_id"]] = (r["event_id"], r["continent"], r["year"])

    rk = defaultdict(Counter)            # raster-status by modality
    for r in rows:
        rk[r["raster"]][r["status"]] += 1

    chip_rows = []
    verdict_ct = Counter()
    anyflag = Counter()
    s2_pair = Counter(); s1_pair = Counter()
    for cid, rs in by_chip.items():
        ev, cont, yr = meta[cid]
        st = {k: rs.get(k, {}).get("status", "missing") for k in KEYS}
        fl = {f: any(rs.get(k, {}).get(f) == "1" for k in KEYS)
              for f in ["is_white", "is_partial", "is_corrupt", "is_empty", "is_speckle", "is_cloudy"]}
        fl["is_missing"] = any(st[k] == "missing" for k in KEYS)
        s2u = "usable" if all(usability(st[k]) != "unusable" for k in ("s2_pre", "s2_post")) else "unusable"
        s1u = "usable" if all(usability(st[k]) != "unusable" for k in ("s1_pre", "s1_post")) else "unusable"
        s2_pair[s2u] += 1; s1_pair[s1u] += 1
        n_bad = sum(1 for k in KEYS if st[k] in UNUSABLE or st[k] in COMPROMISED
                    or st[k] in ("cloudy",))
        # single verdict, worst-first
        if all(st[k] == "missing" or st[k] == "empty" for k in KEYS):
            v = "dead"
        elif s2u == "unusable" and s1u == "unusable":
            v = "both_unusable"
        elif s2u == "unusable":
            v = "s2_unusable"          # SAR-only chip
        elif s1u == "unusable":
            v = "s1_unusable"          # optical-only chip
        elif fl["is_partial"]:
            v = "partial_edge"         # usable but has an edge wedge somewhere
        elif fl["is_cloudy"]:
            v = "cloudy"
        elif any(st[k] in ("speckled", "partial_minor") for k in KEYS):
            v = "minor"
        else:
            v = "clean"
        verdict_ct[v] += 1
        for f, on in fl.items():
            if on:
                anyflag[f] += 1
        chip_rows.append({"chip_id": cid, "event_id": ev, "continent": cont, "year": yr,
                          **{f"{k}_status": st[k] for k in KEYS},
                          "s2_usable": s2u, "s1_usable": s1u,
                          **{f: int(v2) for f, v2 in fl.items()},
                          "n_bad_rasters": n_bad, "verdict": v})

    chip_rows.sort(key=lambda r: r["chip_id"])
    cf = ["chip_id", "event_id", "continent", "year", "s2_pre_status", "s2_post_status",
          "s1_pre_status", "s1_post_status", "s2_usable", "s1_usable", "is_missing", "is_white",
          "is_partial", "is_corrupt", "is_empty", "is_speckle", "is_cloudy", "n_bad_rasters", "verdict"]
    with open(REP / "qc_by_chip.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cf, extrasaction="ignore")
        w.writeheader(); w.writerows(chip_rows)
    print(f"chips: {len(chip_rows)} | verdicts: {dict(verdict_ct)}")
    return rows, rk, chip_rows, verdict_ct, anyflag, s2_pair, s1_pair


# --------------------------------------------------------------------------- #
# report
# --------------------------------------------------------------------------- #
VC = {"clean": "#3f9b52", "minor": "#8fce6a", "cloudy": "#e8a33d", "partial_edge": "#e07b39",
      "s1_unusable": "#7b8794", "s2_unusable": "#c0392b", "both_unusable": "#5a2a82",
      "dead": "#000000"}
VLAB = {"clean": "clean (all 4 ok)", "minor": "minor (speckle/tiny edge)",
        "cloudy": "cloudy but usable", "partial_edge": "partial edge-black",
        "s1_unusable": "S1 unusable (optical ok)", "s2_unusable": "S2 unusable (SAR ok)",
        "both_unusable": "both unusable", "dead": "dead (no imagery)"}


def fig_overview(rk, anyflag, verdict_ct, NC):
    FIG.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(15, 10.5))

    # (a) raster primary-status, stacked by modality
    a0 = ax[0, 0]
    present = [s for s in STATUS_ORDER if any(rk[k][s] for k in KEYS)]
    bottoms = [0] * len(KEYS)
    for s in present:
        vals = [rk[k][s] for k in KEYS]
        a0.bar(KEYS, vals, bottom=bottoms, color=STATUS_C[s], label=s, width=0.7,
               edgecolor=SURFACE, linewidth=0.6)
        bottoms = [b + v for b, v in zip(bottoms, vals)]
    a0.set_title("Per-raster QC status by modality", loc="left", color=INK, pad=10)
    a0.set_ylabel("rasters"); _style(a0)
    a0.legend(frameon=False, fontsize=7.5, labelcolor=INK_MUTED, ncol=2, loc="upper center",
              bbox_to_anchor=(0.5, -0.08))

    # (b) defect rate per modality (% non-ok), the headline read
    a1 = ax[0, 1]
    defect_pct = [100 * (1 - rk[k]["ok"] / sum(rk[k].values())) for k in KEYS]
    bars = a1.bar(KEYS, defect_pct, color=["#2a78d6", "#2a78d6", "#7a5195", "#7a5195"], width=0.6)
    for b, p in zip(bars, defect_pct):
        a1.text(b.get_x() + b.get_width() / 2, p + 0.4, f"{p:.1f}%", ha="center", fontsize=10,
                color=INK_MUTED)
    a1.set_title("Share of rasters with any defect", loc="left", color=INK, pad=10)
    a1.set_ylabel("% of rasters"); a1.set_ylim(0, max(defect_pct) * 1.2); _style(a1)
    a1.grid(axis="y", color=GRID, linewidth=0.8)

    # (c) chips affected by each defect (independent, chip-level)
    a2 = ax[1, 0]
    labs = [("is_partial", "partial / edge-black", "#c0392b"),
            ("is_white", "blank white (cloud)", "#7b8794"),
            ("is_cloudy", "heavy cloud (advisory)", "#e8a33d"),
            ("is_corrupt", "corrupt / saturated", "#5a2a82"),
            ("is_missing", "missing raster", "#a0a0a0"),
            ("is_empty", "empty / off-swath", "#2c3e50"),
            ("is_speckle", "speckle (ignored)", "#c9c9c4")]
    labs.sort(key=lambda x: -anyflag[x[0]])
    yp = range(len(labs))[::-1]
    a2.barh(list(yp), [anyflag[k] for k, _, _ in labs], color=[c for _, _, c in labs], height=0.66)
    a2.set_yticks(list(yp)); a2.set_yticklabels([l for _, l, _ in labs])
    for y, (k, _, _) in zip(yp, labs):
        a2.text(anyflag[k] + NC * 0.006, y, f"{anyflag[k]}  ({anyflag[k]/NC*100:.1f}%)",
                va="center", fontsize=9, color=INK_MUTED)
    a2.set_xlim(0, max(anyflag.values()) * 1.25)
    a2.set_title(f"Chips affected by each defect ({NC} chips)", loc="left", color=INK, pad=10)
    a2.set_xlabel("chips (≥1 raster with the defect)"); _style(a2)
    a2.grid(axis="x", color=GRID, linewidth=0.8); a2.invert_yaxis()

    # (d) chip usability verdict
    a3 = ax[1, 1]
    vorder = [v for v in VC if verdict_ct[v]]
    yp = range(len(vorder))[::-1]
    a3.barh(list(yp), [verdict_ct[v] for v in vorder], color=[VC[v] for v in vorder], height=0.66)
    a3.set_yticks(list(yp)); a3.set_yticklabels([VLAB[v] for v in vorder])
    for y, v in zip(yp, vorder):
        a3.text(verdict_ct[v] + NC * 0.006, y, f"{verdict_ct[v]}  ({verdict_ct[v]/NC*100:.1f}%)",
                va="center", fontsize=9, color=INK_MUTED)
    a3.set_xlim(0, max(verdict_ct.values()) * 1.25)
    a3.set_title("Chip usability verdict", loc="left", color=INK, pad=10)
    a3.set_xlabel("chips"); _style(a3)
    a3.grid(axis="x", color=GRID, linewidth=0.8); a3.invert_yaxis()

    fig.suptitle("Flooded — chip imagery quality control", fontsize=15, fontweight="bold",
                 x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    fig.savefig(FIG / "qc_overview.png", dpi=150)
    plt.close(fig)
    print(f"  -> {FIG/'qc_overview.png'}")
    return present, labs, vorder


def write_stats_md(rows, rk, anyflag, verdict_ct, s2_pair, s1_pair, present, labs, vorder, NC):
    tot_r = len(rows)
    L = ["# Flooded — imagery quality control\n",
         f"Automated QC over all **{NC} chips × 4 imagery rasters = {tot_r} raster slots**. "
         "Every threshold is calibrated on the full corpus; see `qc_spec.json` for the justification "
         "behind each number. Full per-raster detail in `qc_chips.csv`; per-chip rollup in "
         "`qc_by_chip.csv`.\n",
         "## Per-raster status\n",
         "`partial_*` = a Sentinel footprint / orbit swath edge cuts the chip (nodata wedge anchored "
         "to the border); `blank_white` = opaque cloud over ≥80% of the valid area; `corrupt` = "
         "physically-impossible reflectance (S2) or a single-polarisation dropout (S1); `empty` = "
         "≥98% nodata; `cloudy` = 35–80% cloud, still usable (advisory).\n",
         "| status | s2_pre | s2_post | s1_pre | s1_post | total |",
         "|---|--:|--:|--:|--:|--:|"]
    for s in present:
        t = sum(rk[k][s] for k in KEYS)
        L.append(f"| {s} | {rk['s2_pre'][s]} | {rk['s2_post'][s]} | {rk['s1_pre'][s]} | {rk['s1_post'][s]} | {t} |")
    L.append(f"| **any defect** | {sum(rk['s2_pre'].values())-rk['s2_pre']['ok']} | "
             f"{sum(rk['s2_post'].values())-rk['s2_post']['ok']} | "
             f"{sum(rk['s1_pre'].values())-rk['s1_pre']['ok']} | "
             f"{sum(rk['s1_post'].values())-rk['s1_post']['ok']} | "
             f"{tot_r - sum(rk[k]['ok'] for k in KEYS)} |")

    L.append("\n## Chips affected by each defect\n")
    L.append("Independent flags — a chip is counted if **any** of its rasters carries the defect, so "
             "these overlap (a chip can be both partial and cloudy).\n")
    L.append("| defect | chips | share |\n|---|--:|--:|")
    for k, lab, _ in labs:
        L.append(f"| {lab} | {anyflag[k]} | {anyflag[k]/NC*100:.1f}% |")

    L.append("\n## Chip usability verdict\n")
    L.append("Single worst-defect-wins verdict per chip. A chip is **S2-usable** if both its S2 epochs "
             "avoid missing/empty/corrupt/blank_white/severe-partial (same for S1).\n")
    L.append("| verdict | chips | share |\n|---|--:|--:|")
    for v in vorder:
        L.append(f"| {VLAB[v]} | {verdict_ct[v]} | {verdict_ct[v]/NC*100:.1f}% |")
    L.append(f"\n- **S2 change-detection usable** (both epochs): {s2_pair['usable']} chips "
             f"({s2_pair['usable']/NC*100:.1f}%)")
    L.append(f"- **S1 change-detection usable** (both epochs): {s1_pair['usable']} chips "
             f"({s1_pair['usable']/NC*100:.1f}%)")
    clean = verdict_ct["clean"]
    L.append(f"- **Fully clean** (all 4 rasters ok): {clean} chips ({clean/NC*100:.1f}%)\n")
    L.append("![qc overview](figures/qc_overview.png)\n")
    L.append("## How to filter\n")
    L.append("```python\nimport csv\nrows = list(csv.DictReader(open('qc_by_chip.csv')))\n"
             "# chips with a swath-edge wedge in any raster:\n"
             "partial = [r for r in rows if r['is_partial'] == '1']\n"
             "# chips whose optical pair is fully cloud-blank:\n"
             "white = [r for r in rows if r['s2_usable'] == 'unusable' and r['is_white'] == '1']\n"
             "# fully clean chips only:\n"
             "clean = [r for r in rows if r['verdict'] == 'clean']\n```")
    (REP / "QC_STATS.md").write_text("\n".join(L))
    print(f"  -> {REP/'QC_STATS.md'}")
    print(f"\nS2 usable pair: {s2_pair['usable']}/{NC}  S1 usable pair: {s1_pair['usable']}/{NC}  clean: {clean}")


def fig_usability():
    """Two-panel figure: chip usability ladder + the corruption that causes it."""
    SK = ["s2_pre_status", "s2_post_status", "s1_pre_status", "s1_post_status"]
    chip = list(csv.DictReader(open(REP / "qc_by_chip.csv")))
    rast = list(csv.DictReader(open(REP / "qc_chips.csv")))
    N = len(chip)

    # ---- panel A data: usability ladder ----
    def nusable(r): return sum(1 for k in SK if r[k] not in UNUSABLE)
    allok = sum(1 for r in chip if all(r[k] == "ok" for k in SK))
    tier = Counter(nusable(r) for r in chip)
    minor = tier[4] - allok
    LADDER = [
        ("Clean — all 4 rasters ok",        allok,    "#2f8f43"),
        ("Usable — only minor defects",     minor,    "#7cc45c"),
        ("3 of 4 rasters usable",           tier[3],  "#e8a33d"),
        ("2 of 4 rasters usable",           tier[2],  "#e07b39"),
        ("1 of 4 rasters usable",           tier[1],  "#c0392b"),
        ("None usable (dead)",              tier[0],  "#3a3a38"),
    ]
    usable4 = allok + minor

    # ---- panel B data: raster defect composition ----
    rc = Counter(r["status"] for r in rast)
    TOT_R = len(rast)
    nonok = TOT_R - rc["ok"]
    DEFECTS = [  # (status, label, color, hard?)
        ("partial_severe",   "swath-edge black (severe)", "#c0392b", True),
        ("blank_white",      "cloud-blank (white)",       "#7b8794", True),
        ("missing",          "missing file",              "#a0a0a0", True),
        ("corrupt",          "corrupt / saturated",       "#5a2a82", True),
        ("empty",            "off-swath black (empty)",   "#2c3e50", True),
        ("cloudy",           "heavy cloud (usable)",      "#e8a33d", False),
        ("partial_moderate", "swath-edge black (mod.)",   "#e07b39", False),
        ("speckled",         "speckle (usable)",          "#c9c9c4", False),
        ("partial_minor",    "swath-edge black (minor)",  "#f2c14e", False),
    ]
    DEFECTS.sort(key=lambda d: -rc[d[0]])
    hard = sum(rc[s] for s, _, _, h in DEFECTS if h)
    soft = sum(rc[s] for s, _, _, h in DEFECTS if not h)

    # ---- draw (this figure was drawn with a 13 pt title) ----
    with plt.rc_context({"axes.titlesize": 13}):
        fig, (aA, aB) = plt.subplots(1, 2, figsize=(15.5, 6.2))

        yp = list(range(len(LADDER)))[::-1]
        for y, (lab, val, col) in zip(yp, LADDER):
            aA.barh(y, val, color=col, height=0.68)
            aA.text(val + N * 0.008, y, f"{val:,}  ({val/N*100:.1f}%)", va="center", fontsize=10,
                    color=INK_MUTED)
        aA.set_yticks(yp); aA.set_yticklabels([l for l, _, _ in LADDER])
        aA.set_xlim(0, max(v for _, v, _ in LADDER) * 1.2)
        aA.set_title("Chip usability", loc="left", color=INK, pad=24)
        aA.set_xlabel(f"chips (of {N:,})"); _style(aA)
        aA.grid(axis="x", color=GRID, linewidth=0.8)
        aA.annotate(f"all 4 rasters usable: {usable4:,} ({usable4/N*100:.1f}%)",
                    xy=(0, 1.02), xycoords="axes fraction", fontsize=10.5, color="#2f8f43",
                    fontweight="bold")
        aA.annotate(f"S2 pair usable {sum(1 for r in chip if r['s2_usable']=='usable')/N*100:.0f}%  ·  "
                    f"S1 pair usable {sum(1 for r in chip if r['s1_usable']=='usable')/N*100:.0f}%",
                    xy=(0, 0.965), xycoords="axes fraction", fontsize=9, color=INK_MUTED)

        yp = list(range(len(DEFECTS)))[::-1]
        for y, (s, lab, col, h) in zip(yp, DEFECTS):
            aB.barh(y, rc[s], color=col, height=0.68,
                    hatch="" if h else "///", edgecolor=SURFACE, linewidth=0)
            aB.text(rc[s] + TOT_R * 0.004, y, f"{rc[s]:,}", va="center", fontsize=9.5, color=INK_MUTED)
        aB.set_yticks(yp); aB.set_yticklabels([l for _, l, _, _ in DEFECTS])
        aB.set_xlim(0, max(rc[s] for s, *_ in DEFECTS) * 1.18)
        aB.set_title("What corrupts a raster", loc="left", color=INK, pad=24)
        aB.set_xlabel(f"rasters ({nonok:,} of {TOT_R:,} carry a defect)"); _style(aB)
        aB.grid(axis="x", color=GRID, linewidth=0.8)
        aB.annotate(f"hard — makes the raster unusable: {hard:,}",
                    xy=(0, 1.02), xycoords="axes fraction", fontsize=10.5, color="#c0392b",
                    fontweight="bold")
        aB.annotate(f"soft — degraded but still usable: {soft:,}   (hatched)",
                    xy=(0, 0.965), xycoords="axes fraction", fontsize=9, color=INK_MUTED)

        fig.suptitle("Flooded — chip usability & corruption", fontsize=15, fontweight="bold",
                     x=0.02, ha="left")
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        fig.savefig(FIG / "qc_usability.png", dpi=150)
    plt.close(fig)
    print(f"clean {allok} | usable-with-minor {minor} | all-4-usable {usable4} ({usable4/N*100:.1f}%)")
    print(f"  -> {FIG/'qc_usability.png'}")


# --------------------------------------------------------------------------- #
# patch
# --------------------------------------------------------------------------- #
def _f(v):
    return None if v in (None, "") else round(float(v), 6)


def patch_sidecars():
    """Write the `qc` block into every <chip_id>.json in place, from the two QC tables."""
    rast = {}
    for r in csv.DictReader(open(REP / "qc_chips.csv")):
        ent = {"status": r["status"], "nodata_frac": _f(r["nodata_frac"])}
        if r["raster"].startswith("s2"):
            ent["white_frac"] = _f(r["white_frac"])
        if r["note"]:
            ent["note"] = r["note"]
        rast[(r["chip_id"], r["raster"])] = ent

    chip = {r["chip_id"]: r for r in csv.DictReader(open(REP / "qc_by_chip.csv"))}

    patched = missing_row = 0
    for jf in ROOT.glob("*/*.json"):
        try:
            d = json.loads(jf.read_text())
        except Exception:
            continue
        cid = d.get("chip_id", jf.stem)
        c = chip.get(cid)
        if c is None:
            missing_row += 1
            continue
        rasters = {k: rast.get((cid, k), {"status": "missing", "nodata_frac": None}) for k in KEYS}
        statuses = [rasters[k]["status"] for k in KEYS]
        defects = sorted({s for s in statuses if s != "ok"})
        n_usable = sum(1 for s in statuses if s not in UNUSABLE)
        d["qc"] = {
            "all_rasters_ok": all(s == "ok" for s in statuses),
            "all_rasters_usable": c["s2_usable"] == "usable" and c["s1_usable"] == "usable",
            "s2_usable": c["s2_usable"] == "usable",
            "s1_usable": c["s1_usable"] == "usable",
            "n_usable_rasters": n_usable,
            "verdict": c["verdict"],
            "defects": defects,
            "rasters": rasters,
        }
        jf.write_text(json.dumps(d, indent=2))
        patched += 1

    print(f"patched {patched} chips"
          + (f" | {missing_row} chips had no QC row (skipped)" if missing_row else ""))


# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(description="Stage C2 — quality control over chip rasters")
    ap.add_argument("--workers", type=int, default=max(2, (os.cpu_count() or 8) - 2))
    ap.add_argument("--limit", type=int, default=0, help="only the first N chips")
    ap.add_argument("--no-patch", action="store_true", help="do not write the qc block into sidecars")
    ap.add_argument("--no-figures", action="store_true", help="skip figures and QC_STATS.md")
    a = ap.parse_args(argv)

    if not scan_all(a.workers, a.limit):
        sys.exit(f"no chips found under {ROOT} — run 02b_consolidate.py first")

    rows, rk, chip_rows, verdict_ct, anyflag, s2_pair, s1_pair = rollup()
    NC = len(chip_rows)
    if not a.no_figures:
        present, labs, vorder = fig_overview(rk, anyflag, verdict_ct, NC)
        write_stats_md(rows, rk, anyflag, verdict_ct, s2_pair, s1_pair, present, labs, vorder, NC)
    if not a.no_patch:
        patch_sidecars()
    if not a.no_figures:
        fig_usability()
    return 0


if __name__ == "__main__":
    sys.exit(main())
