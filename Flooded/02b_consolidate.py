"""Consolidate every batch into one release folder.


"""
import csv, os, sys
from pathlib import Path

REPO = Path.cwd()
OUT = REPO / "release"
DRY = "--dry-run" in sys.argv

BATCHES = sorted((int(p.name[6:]), p.name) for p in REPO.glob("batch_*")
                 if p.is_dir() and p.name[6:].isdigit())
BATCHES = [(f"b{n:02d}", d) for n, d in BATCHES]

# ---------- 1. master events table ----------
rows, seen, dupes = [], set(), []
for bid, d in BATCHES:
    ev_csv = REPO / d / "manifests" / "events.csv"
    if not ev_csv.exists():
        continue
    for r in csv.DictReader(open(ev_csv, newline="")):
        eid = r["event_id"]
        chipdir = REPO / d / "processed" / "chips" / eid
        n_chips = len(list(chipdir.glob("*_label.tif"))) if chipdir.is_dir() else 0
        n_imaged = len(list(chipdir.glob("*_s2_pre.tif"))) if chipdir.is_dir() else 0
        if eid in seen:
            dupes.append((eid, bid)); continue
        seen.add(eid)
        r["batch"] = bid
        r["n_chips"] = n_chips
        r["n_chips_imaged"] = n_imaged
        r["flooded_crops"] = 1 if n_chips > 0 else 0
        r["_srcdir"] = str(chipdir)          # internal, stripped before write
        rows.append(r)

if not rows:
    sys.exit(f"no batch_NN/manifests/events.csv under {REPO}")
fields = [k for k in rows[0].keys() if k != "_srcdir"]
tot_chips = sum(int(r["n_chips"]) for r in rows)
prod = sum(int(r["flooded_crops"]) for r in rows)
print(f"master: {len(rows)} events | flooded-crop {prod} | 0-chip {len(rows)-prod} | chips {tot_chips}")
if dupes:
    print(f"  !! {len(dupes)} duplicate event_ids skipped: {dupes[:5]}")

if not DRY:
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "events_master.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"  -> {OUT/'events_master.csv'}")

# ---------- 2. merge per-event chip folders (hardlink) ----------
linked_ev = linked_files = copied = 0
for r in rows:
    if int(r["n_chips"]) == 0:
        continue
    src = Path(r["_srcdir"])
    dst = OUT / "chips" / r["event_id"]
    if DRY:
        linked_ev += 1; linked_files += len(list(src.iterdir())); continue
    dst.mkdir(parents=True, exist_ok=True)
    for f in src.iterdir():
        if not f.is_file():
            continue
        t = dst / f.name
        if t.exists():
            continue
        try:
            os.link(f, t)                     # hardlink
            linked_files += 1
        except OSError:
            import shutil; shutil.copy2(f, t); copied += 1
    linked_ev += 1
print(f"chips merged: {linked_ev} event folders, {linked_files} files hardlinked"
      + (f", {copied} copied (cross-device fallback)" if copied else ""))
print("DRY RUN — nothing written" if DRY else f"  -> {OUT/'chips'}")
