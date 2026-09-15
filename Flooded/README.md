# Flooded

Builds 512 × 512 px chips at 10 m of cropland hit by floods. Each chip has a flood label from the
Global Flood Monitoring (GFM) product, Sentinel-2 and Sentinel-1 before and after the flood, the
crop growth stage at the flood date, and a quality-control verdict.

## Requirements

- The conda environment in `environment.yml` at the repository root
- Google Earth Engine: run `earthengine authenticate` once, then set a registered Cloud project:
  ```bash
  export EE_PROJECT=<your-project-id>
  ```
- Optional: `export EODC_TOKEN=<token>`, sent with the GFM reads in Stage B
- Internet access to the EODC STAC API, Zenodo, the Copernicus EMS and IFRC GO APIs, and Earth
  Engine

The scripts read and write in the folder they are run from. Run every step from the same working
folder; every path below is relative to it, and `$REPO` stands for the path to this repository.

## Pipeline

| step | script | writes |
|---|---|---|
| A : events | `01_build_events.py` | the event list of one sampling run, split into batches |
| B–D : chips | `02_build_chips.py` | per batch: flood extent, labels, chips, imagery |
| merge | `02b_consolidate.py` | `release/`: every batch in one folder |
| C1 : phenology | `03_add_phenology.py` | a `phenology` block in every chip sidecar |
| C2 : quality control | `04_quality_control.py` | QC tables and figures, a `qc` block in every sidecar |

Run order: `01` → `02` for each batch → `02b` → `03` → `04`.

---

## 01 — events

```bash
python $REPO/Flooded/01_build_events.py --root <run folder> --n 200 [--append] \
    [--seed <earlier run>/manifests/events_ledger.csv] [--first-batch N] [--shard-size 50]
```

For example, two sampling runs of 200 events each, sharded into batches 1–4 and 5–8:

```bash
python $REPO/Flooded/01_build_events.py --root run_01 --n 200 --first-batch 1
python $REPO/Flooded/01_build_events.py --root run_02 --n 200 --append \
    --seed run_01/manifests/events_ledger.csv --first-batch 5
```

1. **Candidates.** Groundsource news-reported floods (2020–2026, 10–2,000 km²) and DFO flood
   records (severity ≥ 1.5), both downloaded from Zenodo into `data/raw/` on first use.
   Groundsource events are checked against Copernicus EMS activations (within 150 km), DFO
   polygons and IFRC GO flood reports, each within ±15 days; the EMS and IFRC responses are cached
   in `data/raw/crossverify/`. Groundsource events confirmed by EMS or DFO (tier A) and the DFO
   events (tier C) are merged, keeping one event per flood (within 60 km and 15 days).
2. **Sample.** Continents are visited in a seeded round-robin, with years spread within each.
   Each candidate is screened against GFM on the EODC STAC API: its area padded by 8 km, from 4
   days before to 14 days after the event, read at 80 m. It passes with at least one Sentinel-1
   scene and 300 flooded pixels, and its area is then tightened to the flooded extent. Sampling
   stops when `--n` events have passed.
3. **Record and shard.** `--seed` first copies an earlier cumulative list into
   `<root>/manifests/events.csv`. `--append` skips the events already in that file and those in
   `data/raw/screened_noflood.txt` (candidates read cleanly with no flood; the list grows every
   run).

Written to `<root>/manifests/`: `events.csv`, `events.geojson` and a run report (plus
`candidates.csv` without `--append`); `events_new<N>.csv` and `rows.txt`, this run's new events
with and without the header; `events_ledger.csv`, the cumulative list and the `--seed` of the next
run; and `shard_aa`, `shard_ab`, … of up to `--shard-size` events each. With `--first-batch N`,
the shards become batches N, N + 1, …, each written as `batch_NN/manifests/events.csv`. Log:
`<root>/logs/stageA.log`.

## 02 — chips (Stages B, C, D)

```bash
python $REPO/Flooded/02_build_chips.py --root batch_01 [--stages BCD|B|CD] [--workers 8] [--max-attempts 5]
```

Reads `<root>/manifests/events.csv` and runs each stage in its own process.

- **B — flood extent.** Searches GFM on the EODC STAC API over each event's area and window, and
  downloads the flood extent, reference water, exclusion mask and likelihood onto the event's UTM
  grid at 20 m, one file per overpass (`<root>/raw/gfm/<event>/scenes/`). These are combined in
  `<root>/interim/gfm/<event>/`: flood extent (flooded in any overpass), flood persistence (% of
  valid overpasses flooded), reference water, exclusion, peak likelihood, and `flood_series.csv`
  (flooded area per overpass). An event that runs longer than one hour is skipped.
- **C — labels and chips.** Reprojects the layers to 10 m (`<root>/interim/labels_reproj/`) and
  adds cropland from Earth Engine: USDA CDL inside the conterminous US, ESA WorldCover elsewhere.
  The label `<root>/processed/labels/<event>/label_10m.tif` has 5 bands:

  | band | content |
  |---|---|
  | 1 | class: 1 flooded cropland, 2 dry cropland, 3 excluded cropland (GFM exclusion mask), 4 non-cropland, 255 never observed |
  | 2 | GFM likelihood (0–100) |
  | 3 | GFM exclusion mask |
  | 4 | cropland (0/1) |
  | 5 | flood persistence (0–100) |

  The label is cut into 512 px tiles. A tile qualifies with at least 50% observed pixels, 20%
  cropland and 2% flooded cropland (both as shares of the observed pixels). The 5 tiles with the
  most flooded cropland are written to `<root>/processed/chips/<event>/` as `<chip>_label.tif`
  with a `<chip>.json` sidecar, and listed in `<root>/manifests/chips.jsonl` and `chips_gee.csv`.
- **D — imagery.** For every chip, from Earth Engine: Sentinel-2 from the 30 days before the event
  and the 30 days after, each the scene with the highest Cloud Score+ clear fraction over the chip
  (12 bands: B2 B3 B4 B5 B6 B7 B8 B8A B9 B11 B12 NDVI); Sentinel-1 VV + VH, the nearest scene
  before, and the nearest scene after on the same pass and relative orbit. They are written next
  to the label as `<chip>_s2_pre.tif`, `_s2_post.tif`, `_s1_pre.tif` and `_s1_post.tif`, with
  their dates in the sidecar's `imagery` block and in `<root>/manifests/imagery.jsonl`.
  `--workers` chips are fetched in parallel.

If Stage B is killed (exit code 137 or 143, e.g. out of memory), the event it was working on is
removed from `events.csv` (backup `events.csv.bak_<ID>`) and Stage B restarts, up to
`--max-attempts` times. `--stages B` and `--stages CD` run the two halves separately. Re-runs
reuse downloaded scenes (B) and existing labels (C), and skip chips that already have all four
images (D). Logs: `<root>/logs/stageB.log` and `stageCD.log`; run reports in
`<root>/manifests/`.

## 02b — merge

```bash
python $REPO/Flooded/02b_consolidate.py [--dry-run]
```

Reads every `batch_NN` folder that has a `manifests/events.csv`, in batch order. Writes
`release/events_master.csv`, every event with its `batch` (`b01`, `b02`, …), `n_chips`,
`n_chips_imaged` and `flooded_crops` (an id seen twice keeps its first batch), and hard-links each
event's chip files into `release/chips/<event>/`, copying them where hard links are not possible.
`--dry-run` prints the counts and writes nothing.

## 03 — phenology

```bash
python $REPO/Flooded/03_add_phenology.py [--force] [--limit N]
```

Adds a `phenology` block to every sidecar in `release/chips/`, in place. Each pixel's transition
dates come from VIIRS VNP22Q2 across four candidate growth cycles (cycles 1 and 2 of the event
year and of the year before). A cycle counts only if it is complete and in order (greenness
increase < maximum < decrease < minimum). The flood date's position in the first such cycle that
contains it gives the pixel's stage, and the chip takes the majority:

| stage | the flood date falls |
|---|---|
| green-up | between greenness increase and maximum |
| peak/maturity | between maximum and decrease |
| senescence | between decrease and minimum |
| dormant | outside every valid cycle |

The block holds `stage_at_flood`, `stage_agreement`, `stage_fractions`, `flood_date`,
`growth_cycle`, `transitions` and `n_valid_pixels`. VNP22Q2 covers 2013–2024: chips from later
years, and chips with no valid cycle, get `phenology: null` and a `phenology_note`. Chips that
already have the block are skipped unless `--force`; `--limit` processes the first N events.

## 04 — quality control

```bash
python $REPO/Flooded/04_quality_control.py [--workers N] [--limit N] [--no-patch] [--no-figures]
```

Scans the four images of every chip in `release/chips/` and gives each image one status:

| status | meaning |
|---|---|
| `ok` | none of the below |
| `cloudy` | Sentinel-2 cloud over 35–80% of the valid area; still usable |
| `speckled` | scattered no-data pixels; still usable |
| `partial_minor` / `partial_moderate` / `partial_severe` | a no-data wedge from a swath or footprint edge, covering < 5% / 5–15% / ≥ 15% of the chip |
| `blank_white` | Sentinel-2 opaque cloud over ≥ 80% of the valid area |
| `corrupt` | Sentinel-2 median visible reflectance above 10,000, or Sentinel-1 missing one polarisation |
| `empty` | ≥ 98% no-data |
| `unreadable` | cannot be read, wrong size or band count, or constant Sentinel-1 |
| `missing` | file absent |

`missing`, `empty`, `unreadable`, `corrupt`, `blank_white` and `partial_severe` make an image
unusable. Each chip gets one verdict, worst first: `dead` (no imagery), `both_unusable`,
`s2_unusable`, `s1_unusable`, `partial_edge`, `cloudy`, `minor`, `clean`.

Writes `release/report/qc_chips.csv` (one row per image), `qc_by_chip.csv` (one row per chip),
`QC_STATS.md`, `figures/qc_overview.png` and `figures/qc_usability.png`, and a `qc` block in every
sidecar (`all_rasters_ok`, `all_rasters_usable`, `s2_usable`, `s1_usable`, `n_usable_rasters`,
`verdict`, `defects`, `rasters`). `--no-patch` leaves the sidecars untouched, `--no-figures` skips
the figures and `QC_STATS.md`, and `--limit` scans only the first N chips.
