#!/usr/bin/env python3
"""
Stages B, C and D 


  B  download + aggregate   GFM STAC search, AOI-clipped download onto the event's UTM grid,
                            temporal max flood extent + per-overpass flood persistence
  C  label + chips          cropland (USDA CDL in CONUS, ESA WorldCover elsewhere),
                            4 class label, top 5 cropland-flood 512 px chips
  D  imagery                before/after Sentinel-2 (clearest scene over the chip) + Sentinel-1


"""
from __future__ import annotations
import os

for _k, _v in {"GDAL_HTTP_TIMEOUT": "30", "GDAL_HTTP_CONNECTTIMEOUT": "15",
               "GDAL_HTTP_MAX_RETRY": "4", "GDAL_HTTP_RETRY_DELAY": "2",
               "CPL_VSIL_CURL_USE_HEAD": "NO"}.items():
    os.environ[_k] = _v

import argparse, csv, json, logging, re, shutil, signal, subprocess, sys, time, urllib.request
import datetime as dt
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import rasterio
import rioxarray  # noqa: F401  (registers the .rio accessor)
from odc.geo.geobox import GeoBox
from odc.geo.geom import BoundingBox
from odc.stac import load as odc_load
from pyproj import Transformer
from pystac_client import Client
from rasterio.enums import Resampling
from rasterio.windows import Window, transform as win_transform

ROOT = Path.cwd()

# --------------------------------------------------------------------------- #
# configuration: the values flood build ran with
# --------------------------------------------------------------------------- #
CONFIG = {
    "gfm": {"stac_url": "https://stac.eodc.eu/api/v1", "collection": "GFM",
            "assets": ["ensemble_flood_extent", "reference_water_mask", "exclusion_mask",
                       "ensemble_likelihood"]},
    "labeling": {"target_crs": "utm", "target_res_m": 10,
                 "final_classes": {"flooded_cropland": 1, "dry_cropland": 2, "excluded_cropland": 3,
                                   "non_cropland": 4, "nodata": 255},
                 "provenance": "observed", "quality_tier": "T1-auto"},
    "cdl": {"non_crop_codes": [0, 63, 64, 65, 81, 82, 83, 87, 88, 111, 112, 121, 122, 123, 124,
                               131, 141, 142, 143, 152, 176, 190, 195]},
    "chipping": {"size": 512, "stride": 512, "min_valid_fraction": 0.5,
                 "selection": {"chips_per_event": 5, "min_flood_frac": 0.02, "min_cropland_frac": 0.20},
                 "imagery": {"tol_days": 30,
                             "s2_cloud_score_plus": "GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED",
                             "s2_clear_thresh": 0.60, "s2_min_coverage": 0.60}},
}


class Config(dict):
    """Dict subclass with attribute access; resolve_path() places outputs under ROOT/<root>/."""

    def __getattr__(self, name: str) -> Any:
        try:
            return self[name]
        except KeyError as exc:  # pragma: no cover - defensive
            raise AttributeError(name) from exc

    def resolve_path(self, key: str) -> Path:
        """Resolve one of the entries under ``paths:`` against the repository root."""
        rel = self["paths"][key]
        p = Path(rel)
        return p if p.is_absolute() else (ROOT / p)

    def ensure_dirs(self) -> None:
        for key in self.get("paths", {}):
            self.resolve_path(key).mkdir(parents=True, exist_ok=True)


def _wrap(obj: Any) -> Any:
    """Recursively wrap dicts as Config for attribute access."""
    if isinstance(obj, dict):
        return Config({k: _wrap(v) for k, v in obj.items()})
    if isinstance(obj, list):
        return [_wrap(v) for v in obj]
    return obj


def make_cfg(root: str) -> Config:
    cfg = _wrap(json.loads(json.dumps(CONFIG)))
    cfg["paths"] = _wrap({k: f"{root}/{k}" for k in ("raw", "interim", "processed", "manifests", "report")})
    return cfg


# --------------------------------------------------------------------------- #
# logging, raster writers, run report
# --------------------------------------------------------------------------- #
_LOG_FORMAT = "%(asctime)s %(levelname)-7s %(name)s | %(message)s"
_CONFIGURED = False


def get_logger(name: str = "flooded", level: int = logging.INFO) -> logging.Logger:
    global _CONFIGURED
    if not _CONFIGURED:
        logging.basicConfig(level=level, format=_LOG_FORMAT, stream=sys.stderr)
        _CONFIGURED = True
    return logging.getLogger(name)


def write_geotiff(path: str | Path, array: np.ndarray, transform, crs,
                  nodata: float | int | None = None, dtype: str | None = None) -> Path:
    """Write a plain (non-tiled) GeoTIFF. Accepts 2-D (H,W) or 3-D (bands,H,W)."""
    import rasterio
    arr = array if array.ndim == 3 else array[np.newaxis, ...]
    count, height, width = arr.shape
    profile = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": count,
        "dtype": dtype or str(arr.dtype),
        "crs": crs,
        "transform": transform,
        "compress": "deflate",
    }
    if nodata is not None:
        profile["nodata"] = nodata
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(arr.astype(profile["dtype"]))
    return path


def write_cog(path: str | Path, array: np.ndarray, transform, crs,
              nodata: float | int | None = None, dtype: str | None = None, blocksize: int = 256) -> Path:
    """ Cloud-Optimized GeoTIFF (tiled + overviews) """
    from rasterio.io import MemoryFile
    from rio_cogeo.cogeo import cog_translate
    from rio_cogeo.profiles import cog_profiles

    arr = array if array.ndim == 3 else array[np.newaxis, ...]
    count, height, width = arr.shape
    out_dtype = dtype or str(arr.dtype)
    src_profile = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": count,
        "dtype": out_dtype,
        "crs": crs,
        "transform": transform,

        "photometric": "MINISBLACK",
    }
    if nodata is not None:
        src_profile["nodata"] = nodata

    dst_profile = cog_profiles.get("deflate")
    dst_profile.update({"blockxsize": blocksize, "blockysize": blocksize})

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with MemoryFile() as mem:
        with mem.open(**src_profile) as tmp:      # write pass
            tmp.write(arr.astype(out_dtype))
        with mem.open() as src:                    # reopen read-only for translate
            cog_translate(src, path, dst_profile, nodata=nodata, quiet=True,
                          in_memory=True)
    return path


def assert_crs_res(path: str | Path, expected_res_m: float, expected_epsg: int | None = None,
                   tol: float = 1e-3) -> None:
    """Assert a raster's pixel size (and optionally EPSG) after reprojection."""
    import rasterio
    with rasterio.open(path) as src:
        xres, yres = abs(src.transform.a), abs(src.transform.e)
        if abs(xres - expected_res_m) > tol or abs(yres - expected_res_m) > tol:
            raise AssertionError(
                f"{path}: resolution {xres:.4f}x{yres:.4f} m != expected {expected_res_m} m")
        if expected_epsg is not None:
            epsg = src.crs.to_epsg() if src.crs else None
            if epsg != expected_epsg:
                raise AssertionError(f"{path}: CRS EPSG:{epsg} != expected EPSG:{expected_epsg}")


def _utc_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


@dataclass
class RunReport:
    """Accumulates per-run provenance and writes a JSON report."""

    stage: str
    started_at: str = field(default_factory=_utc_stamp)
    params: dict[str, Any] = field(default_factory=dict)
    events: list[dict[str, Any]] = field(default_factory=list)
    failures: list[dict[str, Any]] = field(default_factory=list)
    counters: dict[str, int] = field(default_factory=dict)

    def add_event(self, event_id: str, **provenance: Any) -> None:
        self.events.append({"event_id": event_id, **provenance})

    def add_failure(self, event_id: str, error: str, **context: Any) -> None:
        self.failures.append({"event_id": event_id, "error": error, **context})

    def bump(self, key: str, n: int = 1) -> None:
        self.counters[key] = self.counters.get(key, 0) + n

    def to_dict(self) -> dict[str, Any]:
        return {
            "stage": self.stage,
            "started_at": self.started_at,
            "finished_at": _utc_stamp(),
            "params": self.params,
            "counters": self.counters,
            "n_events": len(self.events),
            "n_failures": len(self.failures),
            "events": self.events,
            "failures": self.failures,
        }

    def write(self, manifests_dir: str | Path) -> Path:
        manifests_dir = Path(manifests_dir)
        manifests_dir.mkdir(parents=True, exist_ok=True)
        out = manifests_dir / f"run_report_{self.stage}_{self.started_at}.json"
        out.write_text(json.dumps(self.to_dict(), indent=2, default=str))
        return out


# --------------------------------------------------------------------------- #
# Earth Engine + grid
# --------------------------------------------------------------------------- #
log_gee = get_logger("flooded.gee")
_INITIALIZED = False


def ee_available() -> bool:
    """True if EE_PROJECT is set"""
    return bool(os.environ.get("EE_PROJECT", "").strip())


def init_ee(project: str | None = None):
    """Initialise Earth Engine with the Cloud project in EE_PROJECT; returns the ee module."""
    global _INITIALIZED
    import ee

    if _INITIALIZED:
        return ee
    project = project or os.environ.get("EE_PROJECT", "").strip()
    if not project:
        raise RuntimeError(
            "Set EE_PROJECT in the environment (export EE_PROJECT=your-project-id).")
    try:
        ee.Initialize(project=project)
    except Exception as exc:  
        raise RuntimeError(
            f"ee.Initialize(project={project!r}) failed: {exc}. Confirm the project is registered "
            "for Earth Engine and that `earthengine authenticate` has been run.") from exc
    _INITIALIZED = True
    log_gee.info("Earth Engine initialised (project=%s)", project)
    return ee


GFM_NA_EPSG = 27705            # native GFM CRS for the North America Equi7 sub-grid
GFM_NATIVE_RES_M = 20


def event_geobox(bbox_wsen, crs=f"EPSG:{GFM_NA_EPSG}", resolution: float = GFM_NATIVE_RES_M) -> GeoBox:
    """GeoBox for an event bbox (lon/lat W,S,E,N) in ``crs`` at ``resolution`` m."""
    west, south, east, north = bbox_wsen
    bb = BoundingBox(west, south, east, north, crs="EPSG:4326").to_crs(crs)
    return GeoBox.from_bbox(bb, resolution=resolution, tight=True)


def utm_epsg_for_bbox(bbox_wsen) -> int:
    """Global UTM EPSG for a bbox centroid (both hemispheres)."""
    lon = (bbox_wsen[0] + bbox_wsen[2]) / 2.0
    lat = (bbox_wsen[1] + bbox_wsen[3]) / 2.0
    zone = int((lon + 180.0) / 6.0) + 1
    return (32600 if lat >= 0 else 32700) + zone


def event_geobox_utm(bbox_wsen, resolution: float = 10.0) -> GeoBox:
    """Per-event geobox in the AOI's UTM zone at ``resolution`` m (global)."""
    epsg = utm_epsg_for_bbox(bbox_wsen)
    return event_geobox(bbox_wsen, crs=f"EPSG:{epsg}", resolution=resolution)


# --------------------------------------------------------------------------- #
# Stage B: GFM search, download, aggregation
# --------------------------------------------------------------------------- #
log_stac = get_logger("flooded.gfm_stac")
log_dl = get_logger("flooded.gfm_download")
log_agg = get_logger("flooded.aggregate")
log_b = get_logger("flooded.run_download")
GFM_NODATA = 255
NODATA = 255
EVENT_TIMEOUT_SEC = 3600   # drop an event that takes long


@dataclass
class SearchResult:
    """Outcome of a look-before-download STAC search for one event."""

    event_id: str
    bbox: list[float]
    date_start: str
    date_end: str
    items: list = field(default_factory=list)

    @property
    def n_scenes(self) -> int:
        return len(self.items)

    @property
    def solar_days(self) -> list[str]:
        days = sorted({it.datetime.strftime("%Y-%m-%d") for it in self.items})
        return days

    @property
    def item_ids(self) -> list[str]:
        return [it.id for it in self.items]

    def summary(self) -> dict:
        return {
            "event_id": self.event_id,
            "n_scenes": self.n_scenes,
            "n_overpass_days": len(self.solar_days),
            "solar_days": self.solar_days,
            "item_ids": self.item_ids,
        }


def open_client(stac_url: str) -> Client:
    return Client.open(stac_url)


def search_event(event_id: str, bbox, date_start: str, date_end: str,
                 stac_url: str = "https://stac.eodc.eu/api/v1", collection: str = "GFM",
                 client: Client | None = None) -> SearchResult:
    """Search GFM by bbox + [date_start, date_end] and log what matched."""
    client = client or open_client(stac_url)
    datetime_range = f"{date_start}T00:00:00Z/{date_end}T23:59:59Z"
    search = client.search(collections=[collection], bbox=list(bbox), datetime=datetime_range)
    items = list(search.items())
    result = SearchResult(event_id, list(bbox), date_start, date_end, items)

    if result.n_scenes == 0:
        log_stac.warning("[%s] 0 scenes for %s / %s..%s — SKIP (S1 timing miss?)",
                         event_id, bbox, date_start, date_end)
    else:
        # per-day scene counts (multiple Equi7 tiles per overpass)
        per_day = defaultdict(int)
        for it in items:
            per_day[it.datetime.strftime("%Y-%m-%d")] += 1
        log_stac.info("[%s] %d scenes across %d overpass-days: %s",
                      event_id, result.n_scenes, len(result.solar_days),
                      ", ".join(f"{d}({per_day[d]})" for d in result.solar_days))
    return result


def _gdal_env(token: str | None):
    """rasterio.Env tuned for windowed COG reads, with optional EODC auth."""
    opts = {
        "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
        "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif",
        "GDAL_HTTP_MAX_RETRY": "3",
        "GDAL_HTTP_RETRY_DELAY": "2",
    }
    if token:
        opts["GDAL_HTTP_HEADERS"] = f"Authorization: Bearer {token}"
    return rasterio.Env(**opts)


def _items_by_day(items) -> dict[str, list[str]]:
    """Map YYYYMMDD -> [item_id, ...] for provenance."""
    out: dict[str, list[str]] = defaultdict(list)
    for it in items:
        out[it.datetime.strftime("%Y%m%d")].append(it.id)
    return dict(out)


def load_event_stack(items, assets, geobox, token: str | None = None, retries: int = 4,
                     backoff: float = 3.0):
    """odc-stac load of ``assets`` onto ``geobox``, grouped by solar day."""
    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with _gdal_env(token):
                ds = odc_load(
                    items,
                    bands=list(assets),
                    geobox=geobox,
                    groupby="solar_day",
                    resampling="nearest",   
                    chunks={"x": 2048, "y": 2048},
                    dtype="uint8",
                    nodata=GFM_NODATA,
                    fail_on_error=False,
                )
            return ds
        except Exception as exc:  # noqa: BLE001 : retry transient network errors
            last_err = exc
            wait = backoff * attempt
            log_dl.warning("odc load attempt %d/%d failed (%s); retry in %.0fs",
                           attempt, retries, type(exc).__name__, wait)
            time.sleep(wait)
    raise RuntimeError(f"odc-stac load failed after {retries} attempts: {last_err}")


def cache_scenes(ds, items, out_dir: Path, assets, geobox, skip_existing: bool = True) -> dict:
    """Write each (overpass-day, asset) slice to a cached GeoTIFF."""
    scenes_dir = Path(out_dir) / "scenes"
    scenes_dir.mkdir(parents=True, exist_ok=True)
    day_items = _items_by_day(items)

    manifest: dict[str, dict] = {}
    times = [np.datetime64(t, "s").astype("datetime64[D]").astype(str) for t in ds["time"].values]
    for i, tstr in enumerate(times):
        day = tstr.replace("-", "")
        entry = manifest.setdefault(day, {"item_ids": day_items.get(day, []), "files": {}})
        for asset in assets:
            out = scenes_dir / f"{day}__{asset}.tif"
            if skip_existing and out.exists():
                entry["files"][asset] = str(out)
                continue
            arr = ds[asset].isel(time=i).values.astype("uint8")
            write_geotiff(out, arr, geobox.transform, geobox.crs,
                          nodata=GFM_NODATA, dtype="uint8")
            entry["files"][asset] = str(out)
        log_dl.info("cached overpass %s (%d assets, %d source items)",
                    day, len(entry["files"]), len(entry["item_ids"]))
    return manifest


def download_event(result, assets, out_dir, geobox, token: str | None = None,
                   skip_existing: bool = True) -> dict:
    """Full per-event download: load stack -> cache per-overpass clips."""
    token = token or os.environ.get("EODC_TOKEN") or None
    if result.n_scenes == 0:
        log_dl.warning("[%s] nothing to download (0 scenes)", result.event_id)
        return {}
    log_dl.info("[%s] downloading %d assets over %d scenes -> %s",
                result.event_id, len(assets), result.n_scenes, out_dir)
    ds = load_event_stack(result.items, assets, geobox, token=token)
    manifest = cache_scenes(ds, result.items, Path(out_dir), assets, geobox,
                            skip_existing=skip_existing)
    return manifest


def _stack_asset(scenes_dir: Path, asset: str):
    """Load all cached per-overpass clips for one asset -> (days, stack, profile)."""
    files = sorted(scenes_dir.glob(f"*__{asset}.tif"))
    if not files:
        return [], None, None
    days, arrs, profile = [], [], None
    for f in files:
        day = f.name.split("__")[0]
        with rasterio.open(f) as src:
            arrs.append(src.read(1))
            if profile is None:
                profile = src.profile.copy()
        days.append(day)
    return days, np.stack(arrs, axis=0), profile


def _pixel_area_km2(profile) -> float:
    xres = abs(profile["transform"].a)
    yres = abs(profile["transform"].e)
    return (xres * yres) / 1e6


def aggregate_event(event_dir, assets, out_dir=None) -> dict:
    """Aggregate one event's cached overpasses. Returns a summary dict."""
    event_dir = Path(event_dir)
    scenes_dir = event_dir / "scenes"
    out_dir = Path(out_dir) if out_dir else event_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- flood extent  ------------------------------
    days, extent, profile = _stack_asset(scenes_dir, "ensemble_flood_extent")
    if extent is None:
        raise FileNotFoundError(f"no ensemble_flood_extent clips under {scenes_dir}")

    valid = extent != NODATA                     # observed in that overpass
    flood = (extent > 0) & valid                 # flooded in that overpass
    any_valid = valid.any(axis=0)
    # temporal max agreement, treating no-data as 0 so it never wins the max
    event_max = np.where(valid, extent, 0).max(axis=0).astype("uint8")
    event_extent = event_max.copy()
    event_extent[~any_valid] = NODATA            # never observed -> no-data
    write_geotiff(out_dir / "event_flood_extent.tif", event_extent,
                  profile["transform"], profile["crs"], nodata=NODATA, dtype="uint8")

    # --- flood persistence (per-pixel severity proxy) ---------------------
    # Fraction of a pixel's VALID overpasses in which it was flooded, 0-100.
    n_valid = valid.sum(axis=0)
    n_flood = flood.sum(axis=0)
    persistence = np.zeros(n_valid.shape, dtype="uint8")
    obs = n_valid > 0
    persistence[obs] = np.round(100.0 * n_flood[obs] / n_valid[obs]).astype("uint8")
    persistence[~any_valid] = NODATA
    write_geotiff(out_dir / "event_flood_persistence.tif", persistence,
                  profile["transform"], profile["crs"], nodata=NODATA, dtype="uint8")

    px_km2 = _pixel_area_km2(profile)
    event_flood_mask = (event_extent > 0) & (event_extent != NODATA)
    pers_flood = persistence[event_flood_mask]
    summary = {
        "n_overpasses": len(days),
        "overpass_days": days,
        "pixel_area_km2": px_km2,
        "event_flood_px": int(event_flood_mask.sum()),
        "event_flood_km2": round(float(event_flood_mask.sum()) * px_km2, 3),
        "event_valid_px": int(any_valid.sum()),
        "event_valid_km2": round(float(any_valid.sum()) * px_km2, 3),
        "flood_persistence_mean": round(float(pers_flood.mean()), 1) if pers_flood.size else 0.0,
    }

    # --- per-overpass flood series (flood evolution) ----------------------
    series_path = out_dir / "flood_series.csv"
    with open(series_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["overpass_day", "flood_px", "valid_px", "flood_km2"])
        for i, day in enumerate(days):
            fpx = int(flood[i].sum())
            vpx = int(valid[i].sum())
            w.writerow([day, fpx, vpx, round(fpx * px_km2, 3)])
    summary["flood_series_csv"] = str(series_path)

    # --- reference water: union over valid overpasses ---------------------
    _, ref, ref_prof = _stack_asset(scenes_dir, "reference_water_mask")
    if ref is not None:
        ref_valid = ref != NODATA
        ref_union = np.where(ref_valid, ref, 0).max(axis=0).astype("uint8")
        ref_union[~ref_valid.any(axis=0)] = NODATA
        write_geotiff(out_dir / "event_reference_water.tif", ref_union,
                      ref_prof["transform"], ref_prof["crs"], nodata=NODATA, dtype="uint8")
        summary["reference_water_px"] = int(((ref_union > 0) & (ref_union != NODATA)).sum())

    # --- exclusion footprint: flagged excluded in ANY overpass ------------
    _, exc, exc_prof = _stack_asset(scenes_dir, "exclusion_mask")
    if exc is not None:
        in_swath = (exc != NODATA).any(axis=0)          # inside a swath at least once
        event_excluded = (exc == 1).any(axis=0)          # flagged excluded somewhere
        out = np.zeros(event_excluded.shape, dtype="uint8")
        out[event_excluded] = 1
        out[~in_swath] = NODATA
        write_geotiff(out_dir / "event_exclusion.tif", out,
                      exc_prof["transform"], exc_prof["crs"], nodata=NODATA, dtype="uint8")
        summary["event_excluded_px"] = int(event_excluded.sum())
        summary["event_excluded_km2"] = round(float(event_excluded.sum()) * px_km2, 3)
        summary["event_inswath_km2"] = round(float(in_swath.sum()) * px_km2, 3)
        summary["excluded_frac_of_inswath"] = round(
            float(event_excluded.sum()) / max(int(in_swath.sum()), 1), 4)

    # --- likelihood: peak over valid overpasses ---------------------------
    _, lik, lik_prof = _stack_asset(scenes_dir, "ensemble_likelihood")
    if lik is not None:
        lik_valid = lik != NODATA
        lik_max = np.where(lik_valid, lik, 0).max(axis=0).astype("uint8")
        lik_max[~lik_valid.any(axis=0)] = NODATA
        write_geotiff(out_dir / "event_likelihood.tif", lik_max,
                      lik_prof["transform"], lik_prof["crs"], nodata=NODATA, dtype="uint8")

    log_agg.info("[%s] aggregated %d overpasses -> flood %.2f km2 (valid %.2f km2)",
                 event_dir.name, summary["n_overpasses"], summary["event_flood_km2"],
                 summary["event_valid_km2"])
    return summary


def _read_manifest(path: Path) -> list[dict]:
    events = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            events.append({
                "event_id": row["event_id"],
                "bbox": [float(row["bbox_w"]), float(row["bbox_s"]),
                         float(row["bbox_e"]), float(row["bbox_n"])],
                "date_start": row["date_start"],
                "date_end": row["date_end"],
            })
    return events


def process_event(ev: dict, cfg, client, report: RunReport) -> None:
    eid = ev["event_id"]
    raw_root = cfg.resolve_path("raw") / "gfm"
    interim_root = cfg.resolve_path("interim") / "gfm"
    assets = cfg["gfm"]["assets"]

    # 1. look before download
    result = search_event(
        eid, ev["bbox"], ev["date_start"], ev["date_end"],
        stac_url=cfg["gfm"]["stac_url"], collection=cfg["gfm"]["collection"], client=client)
    if result.n_scenes == 0:
        report.add_failure(eid, "0 scenes (S1 timing miss)", bbox=ev["bbox"],
                           window=[ev["date_start"], ev["date_end"]])
        report.bump("events_zero_scene")
        return

    # 2. idempotent AOI-clipped download onto the per-event UTM geobox
    gbox = event_geobox_utm(ev["bbox"], cfg["labeling"].get("gfm_native_res_m", 20.0))
    event_raw = raw_root / eid
    manifest = download_event(result, assets, event_raw, gbox)
    report.bump("scenes_downloaded", result.n_scenes)

    # 3. aggregate -> event max-flood-extent + series
    summary = aggregate_event(event_raw, assets, out_dir=interim_root / eid)

    # 4. provenance
    report.add_event(
        eid,
        bbox=ev["bbox"],
        window=[ev["date_start"], ev["date_end"]],
        n_scenes=result.n_scenes,
        search_solar_days=result.solar_days,
        stac_item_ids=result.item_ids,
        gfm_collection=cfg["gfm"]["collection"],
        outputs=str(interim_root / eid),
        **summary,
    )
    report.bump("events_processed")


def stage_b(cfg) -> int:
    """Stage B for every event in the batch's manifests/events.csv."""
    cfg.ensure_dirs()
    events = _read_manifest(cfg.resolve_path("manifests") / "events.csv")
    log_b.info("Stage B: %d event(s) to process", len(events))

    report = RunReport(stage="download", params={
        "n_events": len(events),
        "assets": cfg["gfm"]["assets"],
        "stac_url": cfg["gfm"]["stac_url"],
    })
    client = open_client(cfg["gfm"]["stac_url"])

    # per-event wall-clock budget: an event stuck past EVENT_TIMEOUT_SEC is dropped
    def _on_timeout(signum, frame):
        raise TimeoutError(f"event exceeded {EVENT_TIMEOUT_SEC}s budget")
    signal.signal(signal.SIGALRM, _on_timeout)

    for ev in events:
        try:
            signal.alarm(EVENT_TIMEOUT_SEC)
            process_event(ev, cfg, client, report)
        except Exception as exc:  # noqa: BLE001 - isolate per-event failures (incl. timeout)
            log_b.warning("[%s] SKIPPED (%s): %s", ev["event_id"], type(exc).__name__, exc)
            report.add_failure(ev["event_id"], f"{type(exc).__name__}: {exc}")
        finally:
            signal.alarm(0)

    out = report.write(cfg.resolve_path("manifests"))
    log_b.info("run report -> %s | processed=%d failures=%d scenes=%d",
               out, report.counters.get("events_processed", 0),
               len(report.failures), report.counters.get("scenes_downloaded", 0))
    return 0


# --------------------------------------------------------------------------- #
# Stage C: reproject, cropland, 4-class label, chip selection
# --------------------------------------------------------------------------- #
log_rep = get_logger("flooded.reproject")
log_crop = get_logger("flooded.cropland")
log_lab = get_logger("flooded.build_label")
log_sel = get_logger("flooded.select_chips")
log_c = get_logger("flooded.run_chips")


def target_crs_for(bbox, cfg) -> str:
    tc = str(cfg["labeling"]["target_crs"]).lower()
    if tc == "utm":
        return f"EPSG:{utm_epsg_for_bbox(bbox)}"
    return cfg["labeling"]["target_crs"]   


def _open(path):
    da = rioxarray.open_rasterio(path, masked=False).squeeze("band", drop=True)
    da.rio.write_nodata(NODATA, inplace=True)
    return da


def reproject_event(interim_dir, out_dir, bbox, cfg) -> dict:
    """Reproject all event layers to target CRS @ target_res_m. Returns paths+crs."""
    interim_dir, out_dir = Path(interim_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dst_crs = target_crs_for(bbox, cfg)
    res = float(cfg["labeling"]["target_res_m"])

    # 1. flood extent defines the reference grid 
    fe = _open(interim_dir / "event_flood_extent.tif")
    fe_r = fe.rio.reproject(dst_crs, resolution=res, resampling=Resampling.nearest,
                            nodata=NODATA)
    fe_path = out_dir / "flood_extent_10m.tif"
    fe_r.rio.to_raster(fe_path, dtype="uint8", compress="deflate")

    outputs = {"flood_extent": str(fe_path), "crs": dst_crs, "res_m": res}

    # 2. match the remaining layers onto the flood-extent grid
    layer_resampling = {
        "event_reference_water.tif": ("reference_water", Resampling.nearest),
        "event_exclusion.tif": ("exclusion", Resampling.nearest),
        "event_likelihood.tif": ("likelihood", Resampling.bilinear),
        "event_flood_persistence.tif": ("persistence", Resampling.bilinear),
    }
    for fname, (key, method) in layer_resampling.items():
        src = interim_dir / fname
        if not src.exists():
            continue
        da = _open(src)
        matched = da.rio.reproject_match(fe_r, resampling=method)
        dst = out_dir / f"{key}_10m.tif"
        matched.rio.to_raster(dst, dtype="uint8", compress="deflate")
        outputs[key] = str(dst)

    # 3. assert grid correctness (fail loud if the reprojection is wrong)
    epsg = int(dst_crs.split(":")[1])
    assert_crs_res(fe_path, expected_res_m=res, expected_epsg=epsg)
    log_rep.info("[%s] reprojected -> %s @ %.0f m, grid %s",
                 interim_dir.name, dst_crs, res, tuple(fe_r.shape))
    return outputs


CDL_MIN_YEAR = 2008
_TILE_PX = 2000        # getDownloadURL per-request pixel cap is ~1e7; stay well under


def _download_tiled(img, dst_crs, transform, width, height, res_m, out_dir) -> Path | None:
    """Download a GEE image over the label grid via tiled getDownloadURL + mosaic."""
    import rasterio
    from rasterio.merge import merge as rio_merge

    out_tif = out_dir / "cropland_download.tif"
    if out_tif.exists():                       # idempotent: reuse a prior mosaic
        log_crop.info("cropland download cached: %s", out_tif.name)
        return out_tif

    ee = init_ee()
    proj = ee.Projection(dst_crs)
    tiles, td = [], out_dir / "cropland_tiles"
    td.mkdir(parents=True, exist_ok=True)
    for r0 in range(0, height, _TILE_PX):
        for c0 in range(0, width, _TILE_PX):
            bw, bh = min(_TILE_PX, width - c0), min(_TILE_PX, height - r0)
            x0, y0 = transform.c + c0 * transform.a, transform.f + r0 * transform.e
            x1, y1 = transform.c + (c0 + bw) * transform.a, transform.f + (r0 + bh) * transform.e
            region = ee.Geometry.Rectangle([min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)],
                                           proj, geodesic=False, evenOdd=True)
            tp = td / f"tile_{r0}_{c0}.tif"
            if tp.exists() and tp.stat().st_size > 0:      # idempotent: reuse fetched tiles
                tiles.append(tp)
                continue
            for attempt in range(1, 5):                    # retry transient GEE/network errors
                try:
                    url = img.getDownloadURL({"region": region, "scale": res_m,
                                              "crs": dst_crs, "format": "GEO_TIFF"})
                    with urllib.request.urlopen(url, timeout=120) as resp, open(tp, "wb") as fh:
                        shutil.copyfileobj(resp, fh)   # bounded: a hung read raises -> retry loop
                    break
                except Exception as exc:  # noqa: BLE001
                    if attempt == 4:
                        raise
                    import time
                    log_crop.warning("tile %s_%s attempt %d failed (%s); retrying",
                                     r0, c0, attempt, type(exc).__name__)
                    time.sleep(2.0 * attempt)
            tiles.append(tp)
    if not tiles:
        return None
    srcs = [rasterio.open(t) for t in tiles]
    mosaic, mtransform = rio_merge(srcs)
    prof = srcs[0].profile.copy()
    for s in srcs:
        s.close()
    prof.update(height=mosaic.shape[1], width=mosaic.shape[2], transform=mtransform, count=1)
    out_tif = out_dir / "cropland_download.tif"
    with rasterio.open(out_tif, "w", **prof) as dst:
        dst.write(mosaic[0], 1)
    log_crop.info("cropland tiled download: %d tiles -> %s", len(tiles), out_tif.name)
    return out_tif


# CONUS bbox — inside it we use CDL (crop-type), else WorldCover (binary crop).
_CONUS = (-125.0, 24.0, -66.5, 50.0)


def _is_conus(bbox) -> bool:
    lon = (bbox[0] + bbox[2]) / 2.0
    lat = (bbox[1] + bbox[3]) / 2.0
    w, s, e, n = _CONUS
    return w <= lon <= e and s <= lat <= n


def _ee_cropland_image(bbox, year):
    """Return (ee.Image crop-type, source_tag, kind) for the AOI/year."""
    ee = init_ee()
    if _is_conus(bbox):
        yr = max(int(year), CDL_MIN_YEAR)
        img = (ee.ImageCollection("USDA/NASS/CDL")
               .filter(ee.Filter.calendarRange(yr, yr, "year"))
               .first().select("cropland").toInt())
        return img, f"CDL{yr}", "cdl"
    # ESA WorldCover: pick the version nearest the event year
    asset = "ESA/WorldCover/v100" if int(year) <= 2020 else "ESA/WorldCover/v200"
    img = ee.ImageCollection(asset).first().select("Map").toInt()
    return img, asset.split("/")[-1] + "-WorldCover", "worldcover"


def fetch_cropland(event_id, bbox, year, label_path, profile, cfg, out_dir):
    """Return (cropland_bool, crop_type_int, meta) aligned to the label grid."""
    out_dir = Path(out_dir)
    if not ee_available():
        log_crop.warning("[%s] EE_PROJECT not set — cropland mask unavailable", event_id)
        return None, None, {"cropland_source": "pending_gee"}

    try:
        from rasterio.enums import Resampling
        import rioxarray
        img, source, kind = _ee_cropland_image(bbox, year)
        out_dir.mkdir(parents=True, exist_ok=True)
        raw_tif = _download_tiled(img, str(profile["crs"]), profile["transform"],
                                  profile["width"], profile["height"],
                                  cfg["labeling"]["target_res_m"], out_dir)
        if raw_tif is None:
            raise RuntimeError("no tiles downloaded")
        ref = rioxarray.open_rasterio(label_path).isel(band=0)
        ct = (rioxarray.open_rasterio(raw_tif, masked=False).squeeze("band", drop=True)
              .rio.reproject_match(ref, resampling=Resampling.nearest))
        crop_type = ct.values.astype("int32")
    except Exception as exc:  # noqa: BLE001
        log_crop.warning("[%s] cropland fetch failed (%s) — cropland unknown", event_id, exc)
        return None, None, {"cropland_source": "gee_download_failed"}

    if kind == "cdl":
        non_crop = set(cfg["cdl"]["non_crop_codes"])
        cropland = ~np.isin(crop_type, list(non_crop)) & (crop_type > 0)
    else:                                   # WorldCover class 40 == cropland
        cropland = crop_type == 40

    write_geotiff(out_dir / "cropland_mask_10m.tif", cropland.astype("uint8"),
                  profile["transform"], profile["crs"], nodata=0, dtype="uint8")
    write_geotiff(out_dir / "crop_type_10m.tif", crop_type.astype("uint16"),
                  profile["transform"], profile["crs"], nodata=0, dtype="uint16")
    meta = {"cropland_source": source, "cropland_px": int(cropland.sum())}
    log_crop.info("[%s] cropland: %s -> %d crop px (%.1f%% of AOI)",
                  event_id, source, meta["cropland_px"], 100 * cropland.mean())
    return cropland, crop_type, meta


def _read(path, band=1):
    with rasterio.open(path) as s:
        return s.read(band), s.profile.copy()


def build_four_class(layers: dict, cropland, cfg) -> tuple:
    """Construct (label, likelihood, exclusion, cropland_u8, persistence, profile)."""
    fe, profile = _read(layers["flood_extent"])
    exc = _read(layers["exclusion"])[0] if "exclusion" in layers else np.full_like(fe, NODATA)
    lik = _read(layers["likelihood"])[0] if "likelihood" in layers else np.full_like(fe, NODATA)
    pers = _read(layers["persistence"])[0] if "persistence" in layers else np.full_like(fe, NODATA)

    F = cfg["labeling"]["final_classes"]
    flood = (fe > 0) & (fe != NODATA)
    observed = fe != NODATA
    excluded = exc == 1
    crop = np.zeros(fe.shape, dtype=bool) if cropland is None else cropland.astype(bool)

    label = np.full(fe.shape, F["non_cropland"], dtype="uint8")      # 4 = non-cropland
    label[crop] = F["dry_cropland"]                                  # 2
    label[crop & flood] = F["flooded_cropland"]                     # 1
    label[crop & excluded] = F["excluded_cropland"]                 # 3 (exclusion beats flood)
    label[crop & ~observed] = F["nodata"]                          # 255 (never observed)

    lik = lik.astype("uint8").copy()
    lik[~observed] = NODATA
    pers = pers.astype("uint8").copy()
    pers[~observed] = NODATA                                   # persistence only where observed
    return label, lik, exc.astype("uint8"), crop.astype("uint8"), pers, profile


def class_histogram(label, cfg) -> dict:
    inv = {v: k for k, v in cfg["labeling"]["final_classes"].items()}
    vals, counts = np.unique(label, return_counts=True)
    total = int(label.size)
    return {inv.get(int(v), str(int(v))): {"px": int(c), "pct": round(100 * c / total, 3)}
            for v, c in zip(vals, counts)}


def _crop_stats(label, exc, crop, profile, cfg) -> dict:
    F = cfg["labeling"]["final_classes"]
    px_km2 = (abs(profile["transform"].a) * abs(profile["transform"].e)) / 1e6
    excluded = exc == 1
    inswath = exc != NODATA
    crop_px = int(crop.sum())
    crop_excluded = int((crop.astype(bool) & excluded).sum())
    flooded_crop = int((label == F["flooded_cropland"]).sum())
    return {
        "aoi_excluded_frac_of_inswath": round(float(excluded.sum()) / max(int(inswath.sum()), 1), 4),
        "cropland_km2": round(crop_px * px_km2, 3),
        "cropland_in_exclusion_px": crop_excluded,
        "cropland_in_exclusion_pct": round(100 * crop_excluded / max(crop_px, 1), 3) if crop_px else None,
        "cropland_flooded_km2": round(flooded_crop * px_km2, 3),
    }


def build_label_for_event(event_id, bbox, year, cfg, report: RunReport) -> None:
    interim = cfg.resolve_path("interim") / "gfm" / event_id
    if not (interim / "event_flood_extent.tif").exists():
        report.add_failure(event_id, "no Stage-B aggregate found (run download first)")
        return

    reproj_dir = cfg.resolve_path("interim") / "labels_reproj" / event_id
    layers = reproject_event(interim, reproj_dir, bbox, cfg)

    out_dir = cfg.resolve_path("processed") / "labels" / event_id
    out_dir.mkdir(parents=True, exist_ok=True)

    # cropland on the reprojected flood-extent grid (GEE: WorldCover / CDL)
    fe_profile = _read(layers["flood_extent"])[1]
    cropland, _crop_type, crop_meta = fetch_cropland(
        event_id, bbox, year, layers["flood_extent"], fe_profile, cfg, out_dir)

    label, lik, exc, crop_u8, pers, profile = build_four_class(layers, cropland, cfg)
    hist = class_histogram(label, cfg)
    log_lab.info("[%s] 4-class label: %s", event_id, {k: v["pct"] for k, v in hist.items()})

    stack = np.stack([label, lik, exc, crop_u8, pers], axis=0)   # band5 = flood persistence
    label_path = out_dir / "label_10m.tif"
    write_cog(label_path, stack, profile["transform"], profile["crs"], nodata=NODATA, dtype="uint8")

    stats = _crop_stats(label, exc, crop_u8, profile, cfg)
    report.add_event(
        event_id, target_crs=layers["crs"], res_m=layers["res_m"],
        label_histogram=hist, label_path=str(label_path),
        provenance=cfg["labeling"]["provenance"], quality_tier=cfg["labeling"]["quality_tier"],
        **crop_meta, **stats)


def _as_bool(v):
    return str(v).strip().lower() in ("true", "1", "yes")


def event_provenance(event_meta) -> dict:
    """Confirming-source provenance for a chip's flood event (independent of the news report)."""
    confirmed = ["GFM"]
    if _as_bool(event_meta.get("match_ems")):
        confirmed.append("EMS")
    if _as_bool(event_meta.get("match_dfo")):
        confirmed.append("DFO")
    if _as_bool(event_meta.get("match_ifrc")):
        confirmed.append("IFRC")
    gfm_km2 = event_meta.get("gfm_flood_km2")
    return {
        "origin": "Groundsource (news)" if event_meta.get("source") == "groundsource" else event_meta.get("source"),
        "confirmed_by": confirmed,
        "n_independent": len(confirmed),
        "tier": event_meta.get("tier"),
        "gfm_flood_km2": float(gfm_km2) if gfm_km2 not in (None, "") else None,
    }


def _chip_bbox_lonlat(c0, r0, size, transform, to_ll):
    wt = win_transform(Window(c0, r0, size, size), transform)
    xs = [wt.c, wt.c + size * wt.a]
    ys = [wt.f, wt.f + size * wt.e]
    lons = [to_ll.transform(x, y)[0] for x in xs for y in ys]
    lats = [to_ll.transform(x, y)[1] for x in xs for y in ys]
    return dict(min_lon=round(min(lons), 6), max_lon=round(max(lons), 6),
                min_lat=round(min(lats), 6), max_lat=round(max(lats), 6))


def select_and_write_chips(event_id, event_meta, cfg) -> list[dict]:
    """Score, select and write the top-K cropland-flood chips for one event."""
    F = cfg["labeling"]["final_classes"]
    proc = cfg.resolve_path("processed")
    label_path = proc / "labels" / event_id / "label_10m.tif"
    if not label_path.exists():
        log_sel.warning("[%s] no 4-class label at %s — skip", event_id, label_path)
        return []

    size = int(cfg["chipping"]["size"])
    stride = int(cfg["chipping"]["stride"])
    sel = cfg["chipping"]["selection"]
    k = int(sel["chips_per_event"])
    min_flood = float(sel["min_flood_frac"])
    min_crop = float(sel["min_cropland_frac"])
    min_valid = float(cfg["chipping"]["min_valid_fraction"])

    out_dir = proc / "chips" / event_id
    out_dir.mkdir(parents=True, exist_ok=True)

    with rasterio.open(label_path) as src:
        label = src.read(1)
        stack = src.read()                 # (5, H, W): label, lik, exc, cropland, persistence
        transform, crs = src.transform, src.crs
    to_ll = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    H, W = label.shape

    cands = []
    for r0 in range(0, H - size + 1, stride):
        for c0 in range(0, W - size + 1, stride):
            lab = label[r0:r0 + size, c0:c0 + size]
            valid = lab != NODATA                          # observed pixels
            nvalid = int(valid.sum())
            if nvalid < min_valid * size * size:
                continue
            crop = np.isin(lab, [F["flooded_cropland"], F["dry_cropland"], F["excluded_cropland"]])
            flooded_crop = lab == F["flooded_cropland"]
            crop_frac = crop.sum() / nvalid
            fcrop_frac = flooded_crop.sum() / nvalid
            if crop_frac < min_crop or fcrop_frac < min_flood:
                continue
            rec = {
                "row": r0, "col": c0, "n_valid": nvalid,
                "cropland_frac": round(float(crop_frac), 4),
                "flooded_crop_frac": round(float(fcrop_frac), 4),
                "flooded_crop_px": int(flooded_crop.sum()),
                "excluded_crop_frac": round(float((lab == F["excluded_cropland"]).sum() / nvalid), 4),
            }
            rec.update(_chip_bbox_lonlat(c0, r0, size, transform, to_ll))
            cands.append(rec)

    cands.sort(key=lambda r: r["flooded_crop_frac"], reverse=True)
    kept = cands[:k]

    records = []
    for i, rec in enumerate(kept):
        chip_id = f"{event_id}_r{rec['row']:04d}_c{rec['col']:04d}"
        r0, c0 = rec["row"], rec["col"]
        chip = stack[:, r0:r0 + size, c0:c0 + size]
        wt = win_transform(Window(c0, r0, size, size), transform)
        chip_path = out_dir / f"{chip_id}_label.tif"
        write_cog(chip_path, chip, wt, crs, nodata=NODATA, dtype="uint8")
        # mean flood persistence (band 5) over this chip's flooded-crop pixels — severity proxy
        pers_mean = None
        if chip.shape[0] >= 5:
            fcmask = chip[0] == F["flooded_cropland"]
            if fcmask.any():
                pers_mean = round(float(chip[4][fcmask].mean()), 1)
        meta = {
            "chip_id": chip_id, "event_id": event_id, "rank": i,
            "source": event_meta.get("source"), "tier": event_meta.get("tier"),
            "continent": event_meta.get("continent"),
            "country": event_meta.get("country") or event_meta.get("place"),
            "place": event_meta.get("place") or event_meta.get("country"),
            "provenance": event_provenance(event_meta),
            # the event window drives the before/after imagery fetch (Stage D)
            "event_date_start": event_meta.get("event_date_start") or event_meta.get("date_start"),
            "event_date_end": event_meta.get("event_date_end") or event_meta.get("date_end"),
            "crs": str(crs), "res_m": cfg["labeling"]["target_res_m"], "size": size,
            "label_path": str(chip_path),
            "flood_persistence_mean": pers_mean,
            **{key: rec[key] for key in ("cropland_frac", "flooded_crop_frac",
                                         "flooded_crop_px", "excluded_crop_frac",
                                         "min_lon", "min_lat", "max_lon", "max_lat")},
        }
        side = out_dir / f"{chip_id}.json"
        if side.exists():                       # preserve a prior imagery block (resumable re-runs)
            prev = json.loads(side.read_text())
            if "imagery" in prev:
                meta["imagery"] = prev["imagery"]
        side.write_text(json.dumps(meta, indent=2))
        records.append(meta)

    log_sel.info("[%s] %d candidate chips -> kept %d (flooded_crop_frac %.3f..%.3f, cropland-required)",
                 event_id, len(cands), len(kept),
                 kept[-1]["flooded_crop_frac"] if kept else 0.0,
                 kept[0]["flooded_crop_frac"] if kept else 0.0)
    return records


GEE_FIELDS = ["chip_id", "event_id", "min_lon", "min_lat", "max_lon", "max_lat",
              "event_date_start", "event_date_end", "flooded_crop_frac",
              "flooded_crop_px", "cropland_frac", "continent", "tier"]


def _read_events(cfg, only=None) -> list[dict]:
    path = cfg.resolve_path("manifests") / "events.csv"
    rows = []
    for row in csv.DictReader(open(path, newline="")):
        row["bbox"] = [float(row["bbox_w"]), float(row["bbox_s"]),
                       float(row["bbox_e"]), float(row["bbox_n"])]
        row["year"] = int(str(row.get("event_date_start") or row["date_start"])[:4])
        rows.append(row)
    return [r for r in rows if not only or r["event_id"] == only]


def stage_c(cfg) -> int:
    """Stage C for every event: 4-class label (if not built yet) + chip selection."""
    cfg.ensure_dirs()
    events = _read_events(cfg)
    manifests = cfg.resolve_path("manifests")
    log_c.info("chips: %d event(s)", len(events))

    report = RunReport(stage="chips", params={
        "selection": dict(cfg["chipping"]["selection"]),
        "size": cfg["chipping"]["size"]})

    def _on_timeout(signum, frame):
        raise TimeoutError(f"event exceeded {EVENT_TIMEOUT_SEC}s budget")
    signal.signal(signal.SIGALRM, _on_timeout)

    all_chips: list[dict] = []
    for ev in events:
        try:
            signal.alarm(EVENT_TIMEOUT_SEC)   # per-event budget (huge-AOI cropland can run for hours)
            label_path = cfg.resolve_path("processed") / "labels" / ev["event_id"] / "label_10m.tif"
            if not label_path.exists():
                build_label_for_event(ev["event_id"], ev["bbox"], ev["year"], cfg, report)
            chips = select_and_write_chips(ev["event_id"], ev, cfg)
            all_chips.extend(chips)
            report.add_event(ev["event_id"], n_chips=len(chips),
                             mean_flooded_crop_frac=round(
                                 sum(c["flooded_crop_frac"] for c in chips) / max(len(chips), 1), 4))
            report.bump("chips_written", len(chips))
        except Exception as exc:  # noqa: BLE001 - isolate per-event failures (incl. timeout)
            log_c.warning("[%s] chips SKIPPED (%s)", ev["event_id"], type(exc).__name__)
            report.add_failure(ev["event_id"], f"{type(exc).__name__}: {exc}")
        finally:
            signal.alarm(0)

    with open(manifests / "chips.jsonl", "w") as fh:
        for c in all_chips:
            fh.write(json.dumps(c) + "\n")
    with open(manifests / "chips_gee.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=GEE_FIELDS, extrasaction="ignore")
        w.writeheader()
        for c in all_chips:
            w.writerow({k: c.get(k, "") for k in GEE_FIELDS})

    out = report.write(manifests)
    log_c.info("chips done | events=%d chips=%d failures=%d | report=%s",
               len(report.events), len(all_chips), len(report.failures), out)
    return 0


# --------------------------------------------------------------------------- #
# Stage D: before/after Sentinel-2 + Sentinel-1 per chip
# --------------------------------------------------------------------------- #
log_ba = get_logger("flooded.before_after")
log_d = get_logger("flooded.run_imagery")
S2_BANDS = ["B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B9", "B11", "B12", "NDVI"]
S1_BANDS = ["VV", "VH"]
_TAGS = ("s2_pre", "s2_post", "s1_pre", "s1_post")


def _shift(d: str, days: int) -> str:
    return (dt.date.fromisoformat(str(d)[:10]) + dt.timedelta(days=days)).isoformat()


def _s2_clearest(ee, region, d0, d1, cfg):
    """S2 scene in [d0,d1] with the highest CLEAR FRACTION over the chip, + NDVI."""
    icfg = cfg["chipping"]["imagery"]
    thresh = float(icfg["s2_clear_thresh"])
    min_cov = float(icfg["s2_min_coverage"])
    expected_px = region.area(1).divide(20 * 20)
    min_valid = expected_px.multiply(min_cov)

    s2 = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
          .filterBounds(region).filterDate(d0, d1)
          .linkCollection(ee.ImageCollection(icfg["s2_cloud_score_plus"]), ["cs"]))

    def calc(img):
        clear_frac = (img.select("cs").gte(thresh)          # mean of 0/1 = clear fraction over chip
                      .reduceRegion(ee.Reducer.mean(), region, 20, bestEffort=True).get("cs"))
        valid = (img.select("B4").gt(0)
                 .reduceRegion(ee.Reducer.sum(), region, 20, bestEffort=True).get("B4"))
        ndvi = img.normalizedDifference(["B8A", "B4"]).rename("NDVI")
        return (img.addBands(ndvi)
                .set("clear_frac", ee.Number(clear_frac))
                .set("covers_chip", ee.Number(valid).gte(min_valid)).toFloat())

    coll = (s2.map(calc).filter(ee.Filter.eq("covers_chip", 1))
            .sort("clear_frac", False))                    # clearest-over-chip first
    return ee.Image(coll.first())


def _s1_nearest(ee, region, d0, d1, ascending_first, pass_dir=None, rel_orbit=None):
    """Nearest S1_GRD scene (VV+VH, IW) in [d0,d1]; optionally fixed pass+orbit."""
    coll = (ee.ImageCollection("COPERNICUS/S1_GRD")
            .filterBounds(region).filterDate(d0, d1)
            .filter(ee.Filter.eq("instrumentMode", "IW"))
            .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VV"))
            .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VH")))
    if pass_dir is not None:
        coll = coll.filter(ee.Filter.eq("orbitProperties_pass", pass_dir))
    if rel_orbit is not None:
        coll = coll.filter(ee.Filter.eq("relativeOrbitNumber_start", rel_orbit))
    return ee.Image(coll.sort("system:time_start", ascending_first).first())


def _download_image(ee, img, bands, region, dst_crs, res_m, ref_label_path, out_tif) -> bool:
    """Download an EE image over the chip region, aligned to the chip label grid."""
    import time
    import rioxarray
    from rasterio.enums import Resampling
    if out_tif.exists() and out_tif.stat().st_size > 0:
        return True
    raw = out_tif.with_suffix(".raw.tif")
    for attempt in range(1, 4):
        try:
            url = img.select(bands).getDownloadURL(
                {"region": region, "scale": res_m, "crs": dst_crs, "format": "GEO_TIFF"})
        except Exception as exc:  # noqa: BLE001 - typically 'no image found' (empty collection)
            log_ba.info("  no scene (%s)", type(exc).__name__)
            return False
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(raw, "wb") as f:
                while chunk := r.read(1 << 20):
                    f.write(chunk)
            ref = rioxarray.open_rasterio(ref_label_path).isel(band=0)
            da = (rioxarray.open_rasterio(raw, masked=False)
                  .rio.reproject_match(ref, resampling=Resampling.bilinear))
            import numpy as np
            vals = da.values.astype("float32")
            vals[~np.isfinite(vals)] = 0.0                 # cloud/edge -inf/NaN -> 0 (never broken)
            write_cog(out_tif, vals, ref.rio.transform(), ref.rio.crs,
                      nodata=0, dtype="float32")
            raw.unlink(missing_ok=True)
            return True
        except Exception as exc:  # noqa: BLE001 - retry transient download/read failures
            raw.unlink(missing_ok=True)
            if attempt == 3:
                log_ba.warning("  download failed after %d attempts (%s)", attempt, type(exc).__name__)
                raise
            time.sleep(3.0 * attempt)
    return False


def _img_date(ee, img):
    """Acquisition date (YYYY-MM-DD) of an EE image, or None (retried transiently)."""
    import time
    for attempt in range(3):
        try:
            return ee.Date(img.get("system:time_start")).format("YYYY-MM-dd").getInfo()
        except Exception:  # noqa: BLE001
            time.sleep(1.5 * (attempt + 1))
    return None


def _img_prop(ee, img, name):
    """Numeric property of an EE image (e.g. clear_frac), rounded, or None."""
    import time
    for attempt in range(3):
        try:
            v = img.get(name).getInfo()
            return round(float(v), 4) if v is not None else None
        except Exception:  # noqa: BLE001
            time.sleep(1.5 * (attempt + 1))
    return None


def fetch_chip_imagery(chip_meta, cfg, tol_days=None) -> dict:
    """Download S2/S1 before+after for one chip. Returns provenance dict."""
    ee = init_ee()
    tol = int(tol_days if tol_days is not None else cfg["chipping"]["imagery"].get("tol_days", 30))
    res_m = int(cfg["labeling"]["target_res_m"])
    dst_crs = chip_meta["crs"]
    label_path = Path(chip_meta["label_path"])
    out_dir = label_path.parent
    cid = chip_meta["chip_id"]

    region = ee.Geometry.Rectangle(
        [chip_meta["min_lon"], chip_meta["min_lat"], chip_meta["max_lon"], chip_meta["max_lat"]])
    start, end = chip_meta["event_date_start"], chip_meta["event_date_end"]
    prov = {"chip_id": cid, "flood_start": start, "flood_end": end, "imagery": {}}

    # --- S2 before / after (clearest OVER THE CHIP via Cloud Score+) ---
    for tag, d0, d1 in [("s2_pre", _shift(start, -tol), start),
                        ("s2_post", end, _shift(end, tol))]:
        img = _s2_clearest(ee, region, d0, d1, cfg)
        out = out_dir / f"{cid}_{tag}.tif"
        if _download_image(ee, img, S2_BANDS, region, dst_crs, res_m, label_path, out):
            prov["imagery"][tag] = {"date": _img_date(ee, img),
                                    "clear_frac": _img_prop(ee, img, "clear_frac"),
                                    "path": str(out)}

    # --- S1 before (nearest before) then after on same pass+orbit ---
    s1b = _s1_nearest(ee, region, _shift(start, -tol), start, ascending_first=False)
    pass_dir = rel_orbit = None
    try:
        pass_dir = s1b.get("orbitProperties_pass").getInfo()
        rel_orbit = s1b.get("relativeOrbitNumber_start").getInfo()
    except Exception:  # noqa: BLE001
        pass
    out = out_dir / f"{cid}_s1_pre.tif"
    if pass_dir and _download_image(ee, s1b, S1_BANDS, region, dst_crs, res_m, label_path, out):
        prov["imagery"]["s1_pre"] = {"date": _img_date(ee, s1b), "path": str(out),
                                     "pass": pass_dir, "rel_orbit": rel_orbit}

    s1a = _s1_nearest(ee, region, end, _shift(end, tol), ascending_first=True,
                      pass_dir=pass_dir, rel_orbit=rel_orbit)
    out = out_dir / f"{cid}_s1_post.tif"
    if _download_image(ee, s1a, S1_BANDS, region, dst_crs, res_m, label_path, out):
        prov["imagery"]["s1_post"] = {"date": _img_date(ee, s1a), "path": str(out),
                                      "pass": pass_dir, "rel_orbit": rel_orbit}

    log_ba.info("[%s] imagery: %s", cid,
                ", ".join(f"{k} {v.get('date')}" for k, v in prov["imagery"].items()) or "none found")
    return prov


def _read_chips(cfg, event=None, limit=None) -> list[dict]:
    path = cfg.resolve_path("manifests") / "chips.jsonl"
    chips = [json.loads(line) for line in open(path)]
    if event:
        chips = [c for c in chips if c["event_id"] == event]
    return chips[:limit] if limit else chips


def _already_done(chip) -> bool:
    d = Path(chip["label_path"]).parent
    cid = chip["chip_id"]
    return all((d / f"{cid}_{t}.tif").exists() for t in _TAGS)


def _process_chip(chip, cfg, tol_days, force):
    """Fetch one chip's imagery + refresh its sidecar (touches only this chip's files)."""
    if not force and _already_done(chip):
        return "skipped", None
    prov = fetch_chip_imagery(chip, cfg, tol_days=tol_days)
    side_p = Path(Path(chip["label_path"]).with_suffix("").as_posix().replace("_label", "") + ".json")
    if side_p.exists():
        meta = json.loads(side_p.read_text())
        meta["imagery"] = prov["imagery"]
        side_p.write_text(json.dumps(meta, indent=2))
    return "done", prov


def stage_d(cfg, workers: int) -> int:
    """Stage D for every chip in manifests/chips.jsonl."""
    if not ee_available():
        log_d.error("EE_PROJECT not set — set it in the environment to fetch imagery")
        return 1

    chips = _read_chips(cfg)
    log_d.info("imagery: %d chip(s)", len(chips))
    report = RunReport(stage="imagery", params={"tol_days": cfg["chipping"]["imagery"].get("tol_days", 30)})

    prov_all = []
    init_ee()  # initialise Earth Engine once in the main thread before fan-out
    log_d.info("imagery: fetching with %d parallel workers", workers)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_process_chip, chip, cfg, None, False): chip for chip in chips}
        # consume completions in the MAIN thread only, so report/prov_all mutation
        # stays single-threaded (no lock needed); workers touch only per-chip files
        for fut in as_completed(futs):
            cid = futs[fut]["chip_id"]
            try:
                status, prov = fut.result()
                if status == "skipped":
                    report.bump("skipped_existing")
                    continue
                prov_all.append(prov)
                got = list(prov["imagery"].keys())
                report.add_event(cid, scenes=got)
                report.bump("chips_with_imagery")
                report.bump("scenes_downloaded", len(got))
            except Exception as exc:  # noqa: BLE001 - isolate per-chip failures
                log_d.exception("[%s] imagery FAILED", cid)
                report.add_failure(cid, f"{type(exc).__name__}: {exc}")

    manifests = cfg.resolve_path("manifests")
    with open(manifests / "imagery.jsonl", "w") as fh:
        for p in prov_all:
            fh.write(json.dumps(p) + "\n")
    out = report.write(manifests)
    log_d.info("imagery done | chips=%d scenes=%d failures=%d skipped=%d | report=%s",
               report.counters.get("chips_with_imagery", 0),
               report.counters.get("scenes_downloaded", 0), len(report.failures),
               report.counters.get("skipped_existing", 0), out)
    return 0


# --------------------------------------------------------------------------- #
# batch driver: each stage in its own process; Stage B drops a faulyu event and resumes
# --------------------------------------------------------------------------- #
KILLED = {137: "out of memory (SIGKILL)", 143: "terminated (SIGTERM)"}
EVENT_TAG = re.compile(r"\[([A-Za-z]+[A-Za-z0-9_]*)\]")


def banner(log, msg, stamp=True):
    line = f"== {msg} {dt.datetime.now():%Y-%m-%d %H:%M:%S} ==" if stamp else f"== {msg} =="
    print(line)
    log.write(line + "\n")
    log.flush()


def run_stage(root: str, stage: str, log, workers: int = 8) -> int:
    cmd = [sys.executable, str(Path(__file__).resolve()), "--root", root, "--run-stage", stage,
           "--workers", str(workers)]
    return subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT).returncode


def event_ids(events_csv: Path) -> set:
    return {ln.split(",", 1)[0] for ln in events_csv.read_text().splitlines()[1:] if ln}


def in_flight_event(log_text: str, ids) -> str | None:
    
    last = None
    for m in EVENT_TAG.finditer(log_text):
        if m.group(1) in ids:
            last = m.group(1)
    return last


def raw_lines(path: Path) -> list[bytes]:
   
    return re.findall(rb"[^\n]*\n|[^\n]+$", path.read_bytes())


def drop_event(events_csv: Path, eid: str) -> None:
    """Back up events.csv as events.csv.bak_<ID>, then remove that event's row."""
    shutil.copyfile(events_csv, events_csv.with_name(f"events.csv.bak_{eid}"))
    lines = raw_lines(events_csv)
    events_csv.write_bytes(lines[0] + b"".join(ln for ln in lines[1:] if ln.split(b",", 1)[0] != eid.encode()))


def drive_b(root: str, mdir: Path, logs: Path, max_attempts: int) -> int:
    events_csv = mdir / "events.csv"
    path = logs / "stageB.log"
    with open(path, "a") as log:
        for attempt in range(1, max_attempts + 1):
            banner(log, f"B attempt {attempt} {root}")
            start = path.stat().st_size
            rc = run_stage(root, "B", log)
            if rc == 0:
                banner(log, f"B DONE {root}")
                return 0
            if rc not in KILLED:
                banner(log, f"B FAILED rc={rc} {root}")
                return rc
            seg = path.read_bytes()[start:].decode("utf-8", "replace")
            eid = in_flight_event(seg, event_ids(events_csv))
            if eid is None:
                banner(log, f"B died rc={rc}; in-flight event not identified; stopping", stamp=False)
                return rc
            drop_event(events_csv, eid)
            banner(log, f"B died rc={rc}; dropped poison [{eid}]; resuming", stamp=False)
        banner(log, f"B gave up after {max_attempts} attempts {root}")
    return 1


def drive_cd(root: str, logs: Path, workers: int) -> int:
    with open(logs / "stageCD.log", "a") as log:
        banner(log, f"C+D START {root}")
        for stage in ("C", "D"):
            rc = run_stage(root, stage, log, workers)
            if rc != 0:
                banner(log, f"{stage} FAILED rc={rc} {root}")
                return rc
        banner(log, f"C+D DONE {root}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Stages B, C, D — flood chips for one batch")
    ap.add_argument("--root", required=True, help="batch output folder, e.g. batch_01")
    ap.add_argument("--stages", choices=["BCD", "B", "CD"], default="BCD")
    ap.add_argument("--workers", type=int, default=8, help="Stage D parallel imagery workers")
    ap.add_argument("--max-attempts", type=int, default=5,
                    help="Stage B launches before giving up (the build needed up to 5)")
    ap.add_argument("--run-stage", choices=["B", "C", "D"], help=argparse.SUPPRESS)
    a = ap.parse_args(argv)

    os.chdir(ROOT)
    if a.run_stage:                       
        cfg = make_cfg(a.root)
        if a.run_stage == "B":
            return stage_b(cfg)
        if a.run_stage == "C":
            return stage_c(cfg)
        return stage_d(cfg, a.workers)

    mdir = ROOT / a.root / "manifests"
    events_csv = mdir / "events.csv"
    if not events_csv.exists():
        sys.exit(f"no batch manifest at {events_csv} — run 01_build_events.py with --first-batch first")
    logs = ROOT / a.root / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    if a.stages in ("BCD", "CD") and not ee_available():
        print("warning: EE_PROJECT is not set — stages C and D need Earth Engine")
    print(f"{a.root}: {len(event_ids(events_csv))} events")

    if "B" in a.stages:
        rc = drive_b(a.root, mdir, logs, a.max_attempts)
        if rc:
            return rc
    if a.stages in ("BCD", "CD"):
        return drive_cd(a.root, logs, a.workers)
    return 0


if __name__ == "__main__":
    sys.exit(main())
