# Burnt

Builds 512 × 512 px chips at 10 m of cropland burned by agricultural fires. Each chip has a burn
label from Sentinel-2 dNBR, Sentinel-2 and Sentinel-1 before and after the fire, the crop growth
stage at the burn date, and a quality-control verdict.

## Requirements

- The conda environment in `environment.yml` at the repository root
- Google Earth Engine: run `earthengine authenticate` once, then set a registered Cloud project:
  ```bash
  export EE_PROJECT=<your-project-id>
  ```

The scripts read and write in the folder they are run from. Run every step from the same working
folder; every path below is relative to it, and `$REPO` stands for the path to this repository.

## Pipeline

| step | script | writes |
|---|---|---|
| A : events | `01_build_events.py` | `dataset/manifests/events.csv`, the burn event list |
| B : chips | `02_build_chips.py` | `dataset/processed/chips/<event>/`: labels, imagery, sidecars |
| C1 : phenology | `03_add_phenology.py` | a `phenology` block in every chip sidecar |
| C2 : quality control | `04_quality_control.py` | QC tables and figures, a `qc` block in every sidecar |

Run order: `01` → `02` → `03` → `04`.

---

## 01 — events

```bash
python $REPO/Burnt/01_build_events.py [--n 1600] [--years 2020 2021 2022 2023 2024 2025] \
    [--per-call 40] [--per-cell 8] [--fresh]
```

Scans the world in cells of one month × one year × one continent, month by month, each month
across every year and continent. In each cell:

1. **Fire patches.** The month's MODIS MCD64A1 burned area is split into connected patches of at
   least 20 pixels (500 m) and turned into polygons. Patches of 50–20,000 km² with at least 15%
   cropland (ESA WorldCover) are kept; up to `--per-call` are fetched per cell.
2. **Confirmation.** Within the patch's burn dates ± 45 days, all three must agree: MODIS MCD64A1
   burned area over ≥ 20% of the patch, VIIRS VNP64A1 burned area over ≥ 20%, and ≥ 3 VIIRS
   VNP14A1 active-fire detections.
3. **De-duplication.** A patch is dropped if its centre lies within 25 km, and its first burn day
   within 30 days, of an event already in the list.

At most `--per-cell` events are taken from a cell, until `--n` new events are found; a cell that
fails (for example, times out) is skipped. Events are appended to an existing list unless
`--fresh`. Each event records its first and last burn day (`event_date_start`, `event_date_end`),
its bounding box, country, continent, cropland share and the three confirmation values.

Writes `dataset/manifests/events.csv`, `events.geojson` and `run_report_events.json`.

## 02 — chips

```bash
python $REPO/Burnt/02_build_chips.py [--start N] [--limit N] [--cont-cap 700]
```

Reads `dataset/manifests/events.csv` in order, taking at most `--cont-cap` events per continent;
`--start` and `--limit` select a slice of that list. For each event, in its UTM zone:

- **Sentinel-2 before and after.** Before: 110 to 8 days before the first burn day. After: 10 to
  75 days after the last burn day. Up to 25 of the least cloudy scenes are masked for cloud
  (Cloud Score+ clear ≥ 0.60) and snow, and mosaicked with the clearest scene on top.
- **Label layers.** dNBR = (NBR before − NBR after) × 1000, with NBR = (B8 − B12) / (B8 + B12);
  severity classes start at dNBR 100, 270, 440 and 660. The burn perimeter is the MCD64A1 burned
  area within the burn dates ± 20 days, grown by one pixel.
- **Sentinel-1 before and after.** VV + VH (IW mode), mosaicked with the nearest date on top.
  Before: up to 60 days before (first burn day − 5 days). After: up to 75 days after (last burn
  day + 8 days).
- **Chip placement.** A 100 m map of burned cropland (observed before and after, dNBR ≥ 100,
  cropland, inside the perimeter) gives up to 5 non-overlapping 5.12 km windows with the most
  burned cropland.

A chip is kept if cropland covers at least 20% of it and burned cropland at least 2% of that
cropland. Chip ids are `<event>_r<row>_c<col>`. Files per chip, in
`dataset/processed/chips/<event>/`:

| file | content |
|---|---|
| `<chip>_label.tif` | 5 bands, uint8: class (1 burned cropland, 2 unburned cropland, 3 excluded cropland — no clear view before and after, 4 non-cropland); dNBR scaled as (dNBR + 500) / 8; excluded mask; cropland mask; severity 0–4. 255 = not observed |
| `<chip>_s2_pre.tif`, `_s2_post.tif` | 12 bands, float32: B2 B3 B4 B5 B6 B7 B8 B8A B9 B11 B12 NDVI |
| `<chip>_s1_pre.tif`, `_s1_post.tif` | 2 bands, float32: VV, VH (dB); absent where Sentinel-1 has no acquisition |
| `<chip>.json` | sidecar: event and confirmation details, dates, cropland / burned / excluded fractions, bounds, imagery dates and scene counts, an initial `qc` block |

In the sidecar, `cropland_frac` and `burned_crop_frac` are shares of the whole chip, and
`excluded_crop_frac` is a share of the cropland. An event is marked done (`_done`) once sampled
and existing sidecars are kept, so a re-run continues where it stopped.

## 03 — phenology

```bash
python $REPO/Burnt/03_add_phenology.py [--force] [--limit N]
```

Adds a `phenology` block to every sidecar in `dataset/processed/chips/`, in place. Each pixel's
transition dates come from VIIRS VNP22Q2 across four candidate growth cycles (cycles 1 and 2 of the
event year and of the year before). A cycle counts only if it is complete and in order (greenness
increase < maximum < decrease < minimum). The burn date's position in the first such cycle that
contains it gives the pixel's stage, and the chip takes the majority:

| stage | the burn date falls |
|---|---|
| green-up | between greenness increase and maximum |
| peak/maturity | between maximum and decrease |
| senescence | between decrease and minimum |
| dormant | outside every valid cycle |

The block holds `stage_at_burn`, `stage_agreement`, `stage_fractions`, `growth_cycle`,
`transitions`, `n_valid_pixels`, `burn_date`, `source` and `method`. VNP22Q2 covers 2013–2024:
chips from later years, and chips with no valid cycle, get `phenology: null` and a
`phenology_note`. Chips that already have the block are skipped unless `--force`; `--limit`
processes the first N events.

## 04 — quality control

```bash
python $REPO/Burnt/04_quality_control.py [--workers N] [--limit N] [--no-patch] [--no-figures]
```

Scans the four images of every chip in `dataset/processed/chips/` and gives each image one status:

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

Writes `data/qc_chips.csv` (one row per image), `data/qc_by_chip.csv` (one row per chip),
`stats/QC_STATS.md`, `stats/figures/qc_overview.png` and `stats/figures/qc_usability.png`, and
replaces the `qc` block in every sidecar (`all_rasters_ok`, `all_rasters_usable`, `s2_usable`,
`s1_usable`, `n_usable_rasters`, `verdict`, `defects`, `rasters`). `--no-patch` leaves the
sidecars untouched, `--no-figures` skips the figures and `QC_STATS.md`, and `--limit` scans only
the first N chips.
