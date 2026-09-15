#!/usr/bin/env python3
"""
Stage A 

Groundsource
news + DFO  checked against EMS / IFRC / DFO 
continent-stratified
sample from GFM

"""
from __future__ import annotations
import os

for _k, _v in {"GDAL_HTTP_TIMEOUT": "30", "GDAL_HTTP_CONNECTTIMEOUT": "15",
               "GDAL_HTTP_MAX_RETRY": "4", "GDAL_HTTP_RETRY_DELAY": "2",
               "CPL_VSIL_CURL_USE_HEAD": "NO"}.items():
    os.environ[_k] = _v

import argparse, csv, json, logging, math, random, re, shutil, string, sys, time, urllib.request
import datetime as dt
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import rasterio
from odc.geo.geobox import GeoBox
from odc.geo.geom import BoundingBox
from odc.stac import load as odc_load
from pystac_client import Client
from shapely import wkt
from shapely.geometry import box, mapping

ROOT = Path.cwd()

# --------------------------------------------------------------------------- #
# configuration
# --------------------------------------------------------------------------- #
CONFIG = {
    "scope": {"region": "GLOBAL", "year_min": 2020, "year_max": 2026,
              "bbox": [-180.0, -60.0, 180.0, 75.0]},
    "sampling": {
        "seed": 42,
        "n_events": 10,
        "sources": ["dfo", "groundsource"],
        "strata": ["year", "continent"],
        "groundsource": {
            "parquet": "data/raw/groundsource/groundsource_2026.parquet",
            "zenodo_url": "https://zenodo.org/records/18647054/files/groundsource_2026.parquet?download=1",
            "start_year": 2020, "min_area_km2": 10.0, "max_area_km2": 2000.0},
        "dfo": {"gpkg": "data/raw/dfo/Global_Flood_Records.gpkg", "zenodo_record": 19288171,
                "min_severity": 1.5},
        "crossverify": {
            "ems_url": "https://rapidmapping.emergency.copernicus.eu/backend/dashboard-api/public-activations-info/",
            "ifrc_url": "https://goadmin.ifrc.org/api/v2/event/",
            "ifrc_dtypes": [12, 27], "match_days": 15, "ems_match_km": 150,
            "cache_dir": "data/raw/crossverify"},
        "dedup": {"cluster_km": 60, "cluster_days": 15},
        "gfm_screen": {"aoi_pad_km": 8, "coarse_res_m": 80, "min_flood_pixels": 300,
                       "min_scenes": 1, "window_back_days": 4, "window_fwd_days": 14},
        "tiers": ["A", "B", "C"],
    },
    "gfm": {"stac_url": "https://stac.eodc.eu/api/v1", "collection": "GFM"},
}


class Config(dict):
    

    def __getattr__(self, name: str) -> Any:
        try:
            return self[name]
        except KeyError as exc:  # pragma: no cover - defensive
            raise AttributeError(name) from exc

    def resolve_path(self, key: str) -> Path:
      
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
    random.seed(cfg["sampling"]["seed"])          
    return cfg


# --------------------------------------------------------------------------- #
# logging +  report
# --------------------------------------------------------------------------- #
_LOG_FORMAT = "%(asctime)s %(levelname)-7s %(name)s | %(message)s"
_CONFIGURED = False


def get_logger(name: str = "flooded", level: int = logging.INFO) -> logging.Logger:
    global _CONFIGURED
    if not _CONFIGURED:
        logging.basicConfig(level=level, format=_LOG_FORMAT, stream=sys.stderr)
        _CONFIGURED = True
    return logging.getLogger(name)


def _utc_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


@dataclass
class RunReport:
    """Accumulates and writes a JSON report."""

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
# geo helpers
# --------------------------------------------------------------------------- #
log_geo = get_logger("flooded.geo")

# Coarse continent boxes as (lon_min, lat_min, lon_max, lat_max)
_CONTINENT_BOXES = [
    ("NorthAmerica", (-170, 12, -50, 84)),
    ("CentralAmericaCaribbean", (-95, 5, -58, 27)),
    ("SouthAmerica", (-93, -57, -32, 13)),
    ("Europe", (-25, 34, 45, 72)),
    ("Africa", (-20, -37, 52, 38)),
    ("MiddleEast", (34, 12, 63, 42)),
    ("SouthAsia", (60, 5, 98, 38)),
    ("EastAsia", (98, 18, 150, 55)),
    ("SoutheastAsia", (92, -11, 142, 21)),
    ("CentralAsiaRussia", (45, 40, 180, 78)),
    ("Oceania", (110, -50, 180, -10)),
]


def continent_of(lon: float, lat: float) -> str:
    """Coarse continent label for a lon/lat centroid (stratification only)."""
    for name, (w, s, e, n) in _CONTINENT_BOXES:
        if w <= lon <= e and s <= lat <= n:
            return name
    return "Other"


def pad_bbox_km(bbox, pad_km: float) -> list[float]:
    """Expand a lon/lat bbox by ~pad_km on every side (approximate)."""
    w, s, e, n = bbox
    lat = (s + n) / 2.0
    import math
    dlat = pad_km / 111.0
    dlon = pad_km / (111.0 * max(math.cos(math.radians(lat)), 0.1))
    return [w - dlon, s - dlat, e + dlon, n + dlat]


def download_if_missing(url: str, dest: Path, expected_bytes: int | None = None) -> Path:
    """Idempotent download: skip if present (and size matches when known)."""
    dest = Path(dest)
    if dest.exists() and (expected_bytes is None or dest.stat().st_size == expected_bytes):
        log_geo.info("cached: %s (%d bytes)", dest.name, dest.stat().st_size)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    log_geo.info("downloading %s -> %s", url, dest)
    req = urllib.request.Request(url, headers={"User-Agent": "flooded/0.2"})
    with urllib.request.urlopen(req, timeout=600) as r, open(dest, "wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)
    log_geo.info("downloaded %s (%d bytes)", dest.name, dest.stat().st_size)
    return dest


def utm_epsg_for_bbox(bbox_wsen) -> int:
    """Global UTM EPSG for a bbox centroid (both hemispheres)."""
    lon = (bbox_wsen[0] + bbox_wsen[2]) / 2.0
    lat = (bbox_wsen[1] + bbox_wsen[3]) / 2.0
    zone = int((lon + 180.0) / 6.0) + 1
    return (32600 if lat >= 0 else 32700) + zone


# --------------------------------------------------------------------------- #
# sources: Groundsource (news) and DFO
# --------------------------------------------------------------------------- #
log_gs = get_logger("flooded.groundsource")
log_dfo = get_logger("flooded.dfo")
_ZENODO_FILE = "https://zenodo.org/records/{rec}/files/Global_Flood_Records.gpkg?download=1"


def load_groundsource(cfg) -> gpd.GeoDataFrame:
    """Load Groundsource events filtered to scope + area band, as candidate rows."""
    gcfg = cfg["sampling"]["groundsource"]
    path = Path(gcfg["parquet"])
    download_if_missing(gcfg["zenodo_url"], path)

    pf = pq.ParquetFile(path)
    scalars = pf.read(columns=["uuid", "area_km2", "start_date", "end_date"]).to_pandas()
    scalars["year"] = scalars["start_date"].str[:4].astype(int)

    y0 = max(gcfg["start_year"], cfg["scope"]["year_min"])
    y1 = cfg["scope"]["year_max"]
    mask = (
        scalars["year"].between(y0, y1)
        & scalars["area_km2"].between(gcfg["min_area_km2"], gcfg["max_area_km2"])
    )
    idx = scalars.index[mask]
    log_gs.info("groundsource: %d/%d events pass year %d-%d + area [%.0f,%.0f] km2",
                len(idx), len(scalars), y0, y1, gcfg["min_area_km2"], gcfg["max_area_km2"])

    # Materialize geometry only for surviving rows (WKB -> shapely).
    full = gpd.read_parquet(path, columns=["uuid", "area_km2", "start_date", "end_date", "geometry"])
    g = full.iloc[idx].copy()
    if g.crs is None:
        g.set_crs(4326, inplace=True)

    b = g.geometry.bounds
    lon_c = (b["minx"].values + b["maxx"].values) / 2.0
    lat_c = (b["miny"].values + b["maxy"].values) / 2.0
    rows = {
        "event_id": ["GS" + u[:12] for u in g["uuid"]],
        "source": "groundsource",
        "date_start": g["start_date"].values,
        "date_end": g["end_date"].values,
        "bbox_w": b["minx"].values,
        "bbox_s": b["miny"].values,
        "bbox_e": b["maxx"].values,
        "bbox_n": b["maxy"].values,
        "year": g["start_date"].str[:4].astype(int).values,
        "area_km2_reported": g["area_km2"].values,
        "continent": [continent_of(x, y) for x, y in zip(lon_c, lat_c)],
        "uuid": g["uuid"].values,
        "geometry": g.geometry.values,
    }
    out = gpd.GeoDataFrame(rows, geometry="geometry", crs=4326)
    for col in ("country", "cause", "glide", "magnitude", "notes"):
        out[col] = None
    return out


def load_dfo(cfg) -> gpd.GeoDataFrame:
    """Load DFO events filtered to scope + Severity>=min, as candidate rows."""
    dcfg = cfg["sampling"]["dfo"]
    path = Path(dcfg["gpkg"])
    download_if_missing(_ZENODO_FILE.format(rec=dcfg["zenodo_record"]), path)

    g = gpd.read_file(path)
    g = g.to_crs(4326)
    g["begin"] = pd.to_datetime(g["BeginDate"], errors="coerce")
    g["end"] = pd.to_datetime(g["EndDate"], errors="coerce")
    g["severity"] = pd.to_numeric(g["Severity"], errors="coerce")
    g["year"] = g["begin"].dt.year

    y0, y1 = cfg["scope"]["year_min"], cfg["scope"]["year_max"]
    keep = (
        g["begin"].notna()
        & g["year"].between(y0, y1)
        & (g["severity"] >= dcfg["min_severity"])
        & g.geometry.notna()
        & ~g.geometry.is_empty
    )
    g = g[keep].copy()

    rows = []
    for _, r in g.iterrows():
        w, s, e, n = r.geometry.bounds
        lon, lat = (w + e) / 2.0, (s + n) / 2.0
        rows.append(
            {
                "event_id": f"DFO{int(r['ReportNumber'])}",
                "source": "dfo",
                "date_start": r["begin"].strftime("%Y-%m-%d"),
                "date_end": (r["end"] if pd.notna(r["end"]) else r["begin"]).strftime("%Y-%m-%d"),
                "bbox_w": w, "bbox_s": s, "bbox_e": e, "bbox_n": n,
                "country": r.get("Country"),
                "continent": continent_of(lon, lat),
                "year": int(r["year"]),
                "magnitude": float(r["severity"]),
                "cause": r.get("MainCause"),
                "glide": r.get("GlideNumber") or None,
                "area_km2_reported": pd.to_numeric(str(r.get("Area", "")).replace(",", ""), errors="coerce"),
                "notes": f"DFO Severity {r['Severity']}; {r.get('SubdivisionName') or ''}".strip("; "),
                "geometry": r.geometry,
            }
        )
    out = gpd.GeoDataFrame(rows, geometry="geometry", crs=4326)
    log_dfo.info("dfo: %d events (year %d-%d, severity>=%.1f)", len(out), y0, y1, dcfg["min_severity"])
    return out


# --------------------------------------------------------------------------- #
# independent-catalog corroboration (EMS, DFO, IFRC)
# --------------------------------------------------------------------------- #
log_cv = get_logger("flooded.crossverify")
_UA = {"User-Agent": "flooded/0.2"}


def _get_json(url: str, params: dict | None = None, timeout: int = 60):
    if params:
        from urllib.parse import urlencode
        url = f"{url}?{urlencode(params)}"
    req = urllib.request.Request(url, headers=_UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def fetch_ems(cfg) -> gpd.GeoDataFrame:
    """Copernicus EMS flood activations as centroid points (cached)."""
    ccfg = cfg["sampling"]["crossverify"]
    cache = Path(ccfg["cache_dir"]) / "ems_flood_activations.json"
    if cache.exists():
        results = json.loads(cache.read_text())
    else:
        data = _get_json(ccfg["ems_url"], {"category": "Flood", "limit": 200})
        results = data.get("results", data if isinstance(data, list) else [])
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(results))
    rows = []
    for a in results:
        try:
            pt = wkt.loads(a["centroid"])
        except Exception:  # noqa: BLE001
            continue
        rows.append({
            "code": a.get("code"),
            "date": str(a.get("eventTime", ""))[:10],
            "countries": ";".join(a.get("countries") or []),
            "geometry": pt,
        })
    g = gpd.GeoDataFrame(rows, geometry="geometry", crs=4326)
    g["date"] = pd.to_datetime(g["date"], errors="coerce")
    log_cv.info("ems: %d flood activations (%s..%s)", len(g),
                g["date"].min().date() if len(g) else "-", g["date"].max().date() if len(g) else "-")
    return g


def fetch_ifrc(cfg) -> pd.DataFrame:
    """IFRC GO flood events (country iso3 + start date), paginated + cached."""
    ccfg = cfg["sampling"]["crossverify"]
    y0 = cfg["scope"]["year_min"]
    cache = Path(ccfg["cache_dir"]) / "ifrc_flood_events.json"
    if cache.exists():
        events = json.loads(cache.read_text())
    else:
        events = []
        for dtype in ccfg["ifrc_dtypes"]:
            offset = 0
            while True:
                data = _get_json(ccfg["ifrc_url"], {
                    "dtype": dtype, "limit": 200, "offset": offset,
                    "disaster_start_date__gte": f"{y0}-01-01",
                })
                res = data.get("results", [])
                for e in res:
                    events.append({
                        "id": e.get("id"),
                        "date": str(e.get("disaster_start_date", ""))[:10],
                        "iso3": ";".join(c.get("iso3", "") for c in (e.get("countries") or []) if c.get("iso3")),
                        "name": e.get("name"),
                    })
                if not data.get("next") or not res:
                    break
                offset += 200
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(events))
    df = pd.DataFrame(events)
    if len(df):
        df["date"] = pd.to_datetime(df["date"], errors="coerce")
    log_cv.info("ifrc: %d flood events since %d", len(df), y0)
    return df


def assign_iso3(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Assign iso_a3 to each event by centroid spatial-join to Natural Earth."""
    ne = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))[["iso_a3", "geometry"]]
    cent = gpd.GeoDataFrame(
        {"__i": np.arange(len(gdf))},
        geometry=gpd.points_from_xy(
            (gdf["bbox_w"].values + gdf["bbox_e"].values) / 2.0,
            (gdf["bbox_s"].values + gdf["bbox_n"].values) / 2.0,
        ),
        crs=4326,
    )
    joined = gpd.sjoin(cent, ne, how="left", predicate="within")
    joined = joined.drop_duplicates("__i").set_index("__i").sort_index()
    out = gdf.copy()
    out["iso3"] = joined["iso_a3"].values
    return out


def _within_days(a: pd.Timestamp, start: pd.Timestamp, end: pd.Timestamp, days: int) -> bool:
    lo, hi = start - pd.Timedelta(days=days), end + pd.Timedelta(days=days)
    return (a >= lo) and (a <= hi)


def corroborate(gs: gpd.GeoDataFrame, ems: gpd.GeoDataFrame, dfo: gpd.GeoDataFrame,
                ifrc: pd.DataFrame, cfg) -> gpd.GeoDataFrame:
    """Flag each Groundsource event with independent-catalog matches.

    Adds columns match_ems / match_dfo / match_ifrc (bool) and n_indep (int).
    """
    ccfg = cfg["sampling"]["crossverify"]
    days = ccfg["match_days"]
    gs = gs.reset_index(drop=True).copy()
    gs["ds"] = pd.to_datetime(gs["date_start"], errors="coerce")
    gs["de"] = pd.to_datetime(gs["date_end"], errors="coerce")
    match_ems = np.zeros(len(gs), dtype=bool)
    match_dfo = np.zeros(len(gs), dtype=bool)
    sindex = gs.sindex

    # EMS points: candidate GS polygons containing the point (buffered by match_km).
    buf_deg = ccfg["ems_match_km"] / 111.0
    for _, e in ems.iterrows():
        if pd.isna(e["date"]):
            continue
        pt = e.geometry
        for j in sindex.query(pt.buffer(buf_deg), predicate="intersects"):
            if _within_days(e["date"], gs.at[j, "ds"], gs.at[j, "de"], days):
                match_ems[j] = True

    # DFO polygons: candidate GS polygons intersecting, with temporal overlap.
    dfo = dfo.copy()
    dfo["ds"] = pd.to_datetime(dfo["date_start"], errors="coerce")
    dfo["de"] = pd.to_datetime(dfo["date_end"], errors="coerce")
    for _, d in dfo.iterrows():
        for j in sindex.query(d.geometry, predicate="intersects"):
            # temporal overlap (windowed)
            if (d["ds"] - pd.Timedelta(days=days)) <= gs.at[j, "de"] and \
               (d["de"] + pd.Timedelta(days=days)) >= gs.at[j, "ds"]:
                match_dfo[j] = True

    gs["match_ems"] = match_ems
    gs["match_dfo"] = match_dfo

    # IFRC: country + date (weak). Assign iso3 to GS, then check per-country dates.
    gs["match_ifrc"] = False
    if len(ifrc):
        gs = assign_iso3(gs)
        ifrc_valid = ifrc.dropna(subset=["date"])
        by_iso = {}
        for _, r in ifrc_valid.iterrows():
            for iso in str(r["iso3"]).split(";"):
                if iso:
                    by_iso.setdefault(iso, []).append(r["date"])
        mi = np.zeros(len(gs), dtype=bool)
        for j in range(len(gs)):
            iso = gs.at[j, "iso3"]
            if iso and iso in by_iso:
                for dt in by_iso[iso]:
                    if _within_days(dt, gs.at[j, "ds"], gs.at[j, "de"], days):
                        mi[j] = True
                        break
        gs["match_ifrc"] = mi

    gs["n_indep"] = (gs["match_ems"].astype(int) + gs["match_dfo"].astype(int)
                     + gs["match_ifrc"].astype(int))
    log_cv.info("corroborate: %d GS events with EMS=%d DFO=%d IFRC=%d (any-independent=%d)",
                len(gs), int(match_ems.sum()), int(match_dfo.sum()),
                int(gs["match_ifrc"].sum()), int((gs["n_indep"] > 0).sum()))
    return gs


# --------------------------------------------------------------------------- #
# dedup + candidate pool
# --------------------------------------------------------------------------- #
log_dedup = get_logger("flooded.dedup")
log_frame = get_logger("flooded.event_frame")


def cluster_events(gdf: gpd.GeoDataFrame, cfg) -> gpd.GeoDataFrame:
    """Return one representative row per space-time cluster (adds cluster_id)."""
    dcfg = cfg["sampling"]["dedup"]
    cell = dcfg["cluster_km"] / 111.0                     # deg per cluster cell
    days = dcfg["cluster_days"]

    g = gdf.reset_index(drop=True).copy()
    lon = (g["bbox_w"] + g["bbox_e"]) / 2.0
    lat = (g["bbox_s"] + g["bbox_n"]) / 2.0
    ds = pd.to_datetime(g["date_start"], errors="coerce")
    epoch_bucket = (ds.astype("int64") // (86_400_000_000_000 * days))  # ns -> day-bucket

    g["cluster_id"] = (
        (lon / cell).round().astype("Int64").astype(str) + "_"
        + (lat / cell).round().astype("Int64").astype(str) + "_"
        + epoch_bucket.astype("Int64").astype(str)
    )

    # rank within cluster: n_indep desc, source (dfo first), area desc, onset asc
    g["_area"] = pd.to_numeric(g.get("area_km2_reported"), errors="coerce").fillna(0.0)
    g["_nind"] = pd.to_numeric(g.get("n_indep"), errors="coerce").fillna(0).astype(int)
    g["_dfo"] = (g["source"] == "dfo").astype(int)
    g["_onset"] = ds
    g = g.sort_values(["cluster_id", "_nind", "_dfo", "_area", "_onset"],
                      ascending=[True, False, False, False, True])
    reps = g.drop_duplicates("cluster_id", keep="first").copy()
    reps["cluster_size"] = g.groupby("cluster_id").size().reindex(reps["cluster_id"]).values
    reps = reps.drop(columns=["_area", "_nind", "_dfo", "_onset"])
    log_dedup.info("dedup: %d events -> %d distinct floods (%d clusters)",
                   len(gdf), len(reps), reps["cluster_id"].nunique())
    return reps.reset_index(drop=True)


_SCHEMA = ["event_id", "source", "tier", "date_start", "date_end",
           "bbox_w", "bbox_s", "bbox_e", "bbox_n", "country", "continent", "year",
           "magnitude", "cause", "glide", "area_km2_reported",
           "match_ems", "match_dfo", "match_ifrc", "n_indep", "notes", "geometry"]


def _normalize(g: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    for col in _SCHEMA:
        if col not in g.columns:
            g[col] = None
    return g[_SCHEMA].copy()


def build_candidates(cfg, corroborated_gs: gpd.GeoDataFrame | None = None):
    """Return (candidates, corroborated_gs) — deduped union of Tier-A GS + Tier-C DFO."""
    dfo = load_dfo(cfg)
    dfo["tier"] = "C"

    if corroborated_gs is None:
        gs = load_groundsource(cfg)
        ems = fetch_ems(cfg)
        ifrc = fetch_ifrc(cfg)
        corroborated_gs = corroborate(gs, ems, dfo, ifrc, cfg)

    # Tier A = news events spatially corroborated (EMS point OR DFO polygon).
    tierA = corroborated_gs[corroborated_gs["match_ems"] | corroborated_gs["match_dfo"]].copy()
    tierA["tier"] = "A"
    log_frame.info("candidates: DFO(TierC)=%d + Groundsource-corroborated(TierA)=%d", len(dfo), len(tierA))

    union = gpd.GeoDataFrame(
        pd.concat([_normalize(tierA), _normalize(dfo)], ignore_index=True),
        geometry="geometry", crs=4326,
    )
    union["year"] = union["date_start"].str[:4].astype(int)
    candidates = cluster_events(union, cfg)
    return candidates, corroborated_gs


# --------------------------------------------------------------------------- #
# stratified sampling with lazy GFM screening
# --------------------------------------------------------------------------- #
log_strat = get_logger("flooded.stratify")


def _rank_key(row) -> tuple:
    # within a bucket: prefer more independent matches, then Tier A, then bigger.
    tier_rank = {"A": 0, "B": 1, "C": 0}.get(row.get("tier"), 2)
    nind = int(row.get("n_indep") or 0)
    area = float(row.get("area_km2_reported") or 0.0)
    return (-nind, tier_rank, -area)


def _spread_years(rows, rng) -> list:
    """Order a continent's candidates to round-robin across years (temporal spread)."""
    by_year: dict[int, list] = {}
    for r in rows:
        by_year.setdefault(int(r["year"]), []).append(r)
    for yr in by_year.values():
        yr.sort(key=_rank_key)
    years = sorted(by_year)
    rng.shuffle(years)
    out, i = [], 0
    while any(by_year[y] for y in years):
        y = years[i % len(years)]
        if by_year[y]:
            out.append(by_year[y].pop(0))
        i += 1
    return out


def _plan_order(candidates, seed: int):
    """Continent-primary round-robin (max geographic spread), years spread within."""
    rng = random.Random(seed)
    groups: dict[str, list] = {}
    for _, r in candidates.iterrows():
        groups.setdefault(str(r["continent"]), []).append(r.to_dict())
    for cont in groups:
        groups[cont] = _spread_years(groups[cont], rng)
    order = sorted(groups.keys())
    rng.shuffle(order)                  
    return groups, order


def sample_and_screen(candidates, cfg, screen_fn, tierB_pool=None) -> tuple[list, list]:
    """Draw up to n_events GFM-passing events, spread across strata."""
    seed = cfg["sampling"]["seed"]
    n_events = cfg["sampling"]["n_events"]
    groups, order = _plan_order(candidates, seed)

    kept, screen_log, seen = [], [], set()
    i = 0
    # round-robin across strata, screening lazily until n_events pass
    while len(kept) < n_events and any(groups[k] for k in order):
        k = order[i % len(order)]
        i += 1
        if not groups[k]:
            continue
        row = groups[k].pop(0)
        if row["event_id"] in seen:
            continue
        seen.add(row["event_id"])
        try:
            res = screen_fn(row)
        except Exception as exc:  # noqa: BLE001 - isolate per-candidate screen failures
            log_strat.warning("[%s] screen errored (%s) — skipping candidate",
                              row["event_id"], type(exc).__name__)
            screen_log.append({"event_id": row["event_id"], "continent": k,
                               "year": row.get("year"), "error": str(exc)})
            continue
        screen_log.append({"event_id": row["event_id"], "continent": k, "year": row.get("year"),
                           "n_scenes": res.get("n_scenes"), "flood_km2": res.get("flood_km2"),
                           "passed": res.get("passed"), "outcome": res.get("outcome")})
        if res.get("passed"):
            merged = {**row, **{f"gfm_{key}": res.get(key) for key in
                                ("n_scenes", "n_overpass_days", "flood_px", "flood_km2", "epsg")}}
            # tighten AOI to the actual flood footprint when available (DFO fix)
            fb = res.get("flood_bbox")
            if fb:
                merged["bbox_w"], merged["bbox_s"], merged["bbox_e"], merged["bbox_n"] = fb
                merged["aoi_tightened"] = True
            # keep original news/event dates for provenance; the download window is
            # the padded search window that actually captured the flood.
            merged["event_date_start"], merged["event_date_end"] = row["date_start"], row["date_end"]
            merged["date_start"] = res.get("search_start", row["date_start"])
            merged["date_end"] = res.get("search_end", row["date_end"])
            merged["gfm_solar_days"] = res.get("solar_days")
            kept.append(merged)

    if len(kept) < n_events:
        log_strat.warning("only %d/%d events passed the GFM gate from the corroborated pool",
                          len(kept), n_events)
    log_strat.info("sample_and_screen: %d events kept (screened %d candidates across %d strata)",
                   len(kept), len(screen_log), len(groups))
    return kept, screen_log


# --------------------------------------------------------------------------- #
# GFM: STAC search + coarse read 
# --------------------------------------------------------------------------- #
log_stac = get_logger("flooded.gfm_stac")
log_dl = get_logger("flooded.gfm_download")
log_screen = get_logger("flooded.gfm_screen")
GFM_NODATA = 255


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
                    resampling="nearest",   # class rasters must not interpolate
                    chunks={"x": 2048, "y": 2048},
                    dtype="uint8",
                    nodata=GFM_NODATA,
                    fail_on_error=False,
                )
            return ds
        except Exception as exc:  # noqa: BLE001 - retry transient network errors
            last_err = exc
            wait = backoff * attempt
            log_dl.warning("odc load attempt %d/%d failed (%s); retry in %.0fs",
                           attempt, retries, type(exc).__name__, wait)
            time.sleep(wait)
    raise RuntimeError(f"odc-stac load failed after {retries} attempts: {last_err}")


def pad_window(date_start: str, date_end: str, back_days: int, fwd_days: int) -> tuple[str, str]:
    """Widen a news/event window to catch the nearest post-onset S1 overpass."""
    d0 = dt.date.fromisoformat(str(date_start)[:10]) - dt.timedelta(days=back_days)
    d1 = dt.date.fromisoformat(str(date_end)[:10]) + dt.timedelta(days=fwd_days)
    return d0.isoformat(), d1.isoformat()


def screen_event(event_id, bbox, date_start, date_end, cfg, client=None) -> dict:
    """Return GFM-signal screen result for one event (does NOT download scenes)."""
    scfg = cfg["sampling"]["gfm_screen"]
    aoi = pad_bbox_km(bbox, scfg["aoi_pad_km"])
    res = float(scfg["coarse_res_m"])
    px_km2 = (res / 1000.0) ** 2

    win_start, win_end = pad_window(date_start, date_end,
                                    scfg["window_back_days"], scfg["window_fwd_days"])

    result = {
        "event_id": event_id, "n_scenes": 0, "n_overpass_days": 0,
        "flood_px": 0, "flood_km2": 0.0, "gfm_epsg": None, "passed": False,
        # outcome: 'transient' (STAC/read error or 0-scene timing miss -> retryable) by
        # default; set to 'pass'/'no_flood' only after a clean GFM read below.
        "outcome": "transient",
        "solar_days": [], "item_ids": [], "flood_bbox": None,
        "search_start": win_start, "search_end": win_end,
    }

    # STAC search with retry on transient network errors (EODC occasionally drops).
    sr = None
    for attempt in range(1, 5):
        try:
            sr = search_event(event_id, aoi, win_start, win_end,
                              stac_url=cfg["gfm"]["stac_url"],
                              collection=cfg["gfm"]["collection"], client=client)
            break
        except Exception as exc:  # noqa: BLE001
            if attempt == 4:
                log_screen.warning("[%s] STAC search failed after %d attempts (%s)",
                                   event_id, attempt, type(exc).__name__)
                return result
            time.sleep(2.0 * attempt)
    if sr is None:
        return result
    result.update(n_scenes=sr.n_scenes, n_overpass_days=len(sr.solar_days),
                  solar_days=sr.solar_days, item_ids=sr.item_ids)
    if sr.n_scenes < scfg["min_scenes"]:
        return result

    epsg = utm_epsg_for_bbox(bbox)
    result["gfm_epsg"] = epsg
    bb = BoundingBox(*aoi, crs="EPSG:4326").to_crs(f"EPSG:{epsg}")
    geobox = GeoBox.from_bbox(bb, resolution=res, tight=True)

    try:
        ds = load_event_stack(sr.items, ["ensemble_flood_extent"], geobox)
        arr = ds["ensemble_flood_extent"]
        # temporal union: flooded if flooded (==1) in ANY overpass; 255 = nodata.
        flooded = (arr == 1).any("time").compute()
        flood_px = int(flooded.sum())
    except Exception as exc:  # noqa: BLE001
        log_screen.warning("[%s] screen read failed (%s) — treating as no signal",
                           event_id, type(exc).__name__)
        return result

    result["flood_px"] = flood_px
    result["flood_km2"] = round(flood_px * px_km2, 3)
    result["passed"] = flood_px >= scfg["min_flood_pixels"]
    result["outcome"] = "pass" if result["passed"] else "no_flood"   # clean read = definitive

    # Tighten the AOI to the actual flooded extent (DFO footprints are region-wide).
    if flood_px > 0:
        ys, xs = flooded.values.nonzero()
        xcoords = flooded["x"].values[xs]
        ycoords = flooded["y"].values[ys]
        fb = BoundingBox(xcoords.min(), ycoords.min(), xcoords.max(), ycoords.max(),
                         crs=f"EPSG:{epsg}").to_crs("EPSG:4326")
        result["flood_bbox"] = [round(v, 5) for v in (fb.left, fb.bottom, fb.right, fb.top)]
    log_screen.info("[%s] screen: %d scenes / %d days -> %d flood px (%.1f km2) %s",
                    event_id, sr.n_scenes, len(sr.solar_days), flood_px, result["flood_km2"],
                    "PASS" if result["passed"] else "fail")
    return result


# --------------------------------------------------------------------------- #
# the sampled manifest
# --------------------------------------------------------------------------- #
log_bm = get_logger("flooded.build_manifest")

# Persistent, cross-cycle ledger of candidates that got a CLEAN GFM read but no flood
#  Excluded from future --append sampling so they aren't re-screened
# every cycle. Transient/error fails are NOT recorded here, so they stay re-tryable.
NOFLOOD_LEDGER = ROOT / "data" / "raw" / "screened_noflood.txt"

CSV_FIELDS = ["event_id", "source", "tier", "bbox_w", "bbox_s", "bbox_e", "bbox_n",
              "date_start", "date_end", "event_date_start", "event_date_end",
              "year", "continent", "country", "magnitude",
              "cause", "glide", "area_km2_reported", "match_ems", "match_dfo",
              "match_ifrc", "n_indep", "aoi_tightened", "gfm_n_scenes",
              "gfm_n_overpass_days", "gfm_flood_px", "gfm_flood_km2", "notes"]


def _read_ledger(path: Path) -> set:
    if not path.exists():
        return set()
    return {ln.strip() for ln in path.read_text().splitlines() if ln.strip()}


def _append_ledger(path: Path, ids) -> int:
    """Add new ids to the ledger (dedup); return how many were newly added."""
    fresh = sorted(set(ids) - _read_ledger(path))
    if fresh:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "a") as fh:
            fh.write("\n".join(fresh) + "\n")
    return len(fresh)


def write_csv(rows, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in CSV_FIELDS})


def write_geojson(rows, path: Path) -> None:
    feats = []
    for r in rows:
        geom = box(float(r["bbox_w"]), float(r["bbox_s"]),
                   float(r["bbox_e"]), float(r["bbox_n"]))
        props = {k: r.get(k) for k in ("event_id", "source", "tier", "date_start",
                 "date_end", "continent", "country", "n_indep", "gfm_flood_km2")}
        feats.append({"type": "Feature", "geometry": mapping(geom), "properties": props})
    path.write_text(json.dumps({"type": "FeatureCollection", "features": feats}, indent=2))


def _read_existing(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return list(csv.DictReader(open(path, newline="")))


def sample_events(cfg, append: bool) -> int:
    """Candidates -> stratified sample with the GFM gate -> manifests/events.csv."""
    cfg.ensure_dirs()
    manifests_dir = cfg.resolve_path("manifests")

    report = RunReport(stage="sampling", params={
        "sources": cfg["sampling"]["sources"], "scope": dict(cfg["scope"]),
        "seed": cfg["sampling"]["seed"], "n_events": cfg["sampling"]["n_events"],
        "continent": None, "append": append,
        "gfm_screen": dict(cfg["sampling"]["gfm_screen"]),
    })

    candidates, _ = build_candidates(cfg)
    if append:                             # exclude already-sampled + known no-flood -> next batch
        existing_ids = {r["event_id"] for r in _read_existing(manifests_dir / "events.csv")}
        noflood_ids = _read_ledger(NOFLOOD_LEDGER)
        skip = existing_ids | noflood_ids
        candidates = candidates[~candidates["event_id"].isin(skip)].reset_index(drop=True)
        log_bm.info("append: excluded %d already-sampled + %d known no-flood -> %d candidates remain",
                    len(existing_ids), len(noflood_ids), len(candidates))
    report.bump("candidate_floods", len(candidates))

    client = open_client(cfg["gfm"]["stac_url"])
    screen_fn = lambda row: screen_event(  # noqa: E731
        row["event_id"], [row["bbox_w"], row["bbox_s"], row["bbox_e"], row["bbox_n"]],
        row["date_start"], row["date_end"], cfg, client=client)

    kept, screen_log = sample_and_screen(candidates, cfg, screen_fn)
    report.bump("events_sampled", len(kept))
    report.bump("candidates_screened", len(screen_log))

    # grow the no-flood ledger with candidates that got a clean read but no flood
    noflood_new = [r["event_id"] for r in screen_log if r.get("outcome") == "no_flood"]
    added = _append_ledger(NOFLOOD_LEDGER, noflood_new)
    log_bm.info("no-flood ledger: +%d new / %d screened no-flood this run -> %d total",
                added, len(noflood_new), len(_read_ledger(NOFLOOD_LEDGER)))

    events_path = manifests_dir / "events.csv"
    if append:
        existing = _read_existing(events_path)
        have = {r["event_id"] for r in existing}
        new = [r for r in kept if r["event_id"] not in have]
        merged = existing + new
        write_csv(merged, events_path)
        (manifests_dir / "last_batch_events.txt").write_text("\n".join(r["event_id"] for r in new))
        log_bm.info("appended %d new event(s) to %d existing -> %d total", len(new), len(existing), len(merged))
        kept = merged
    else:
        write_csv(candidates.to_dict("records"), manifests_dir / "candidates.csv")
        write_csv(kept, events_path)
    write_geojson(kept, manifests_dir / "events.geojson")
    for r in kept:
        report.add_event(r["event_id"], source=r["source"], tier=r["tier"],
                         bbox=[r["bbox_w"], r["bbox_s"], r["bbox_e"], r["bbox_n"]],
                         window=[r["date_start"], r["date_end"]],
                         continent=r["continent"], n_indep=r.get("n_indep"),
                         gfm_flood_km2=r.get("gfm_flood_km2"))
    report.params["screen_log"] = screen_log

    out = report.write(manifests_dir)
    log_bm.info("manifest -> %s (%d events, %d candidates screened) | report=%s",
                manifests_dir / "events.csv", len(kept), len(screen_log), out)
    return 0


# --------------------------------------------------------------------------- #
# cycle records + batch shards
# --------------------------------------------------------------------------- #
def raw_lines(path: Path) -> list[bytes]:
    """The file's lines, split at LF with their endings untouched (events.csv is written by the
    csv module, so CRLF) — the way split(1) and grep read them."""
    return re.findall(rb"[^\n]*\n|[^\n]+$", path.read_bytes())


def row_id(line: bytes) -> str:
    return line.split(b",", 1)[0].decode()


def new_rows(events_csv: Path, new_ids):
    """(header, rows): this cycle's rows, in events.csv order, as raw lines."""
    lines = raw_lines(events_csv)
    ids = set(new_ids)
    return lines[0], [ln for ln in lines[1:] if row_id(ln) in ids]


def shard_names(n):
    """split(1) suffixes: shard_aa, shard_ab, ..., shard_az, shard_ba, ..."""
    ab = string.ascii_lowercase
    return [f"shard_{a}{b}" for a in ab for b in ab][:n]


def batch_root(nn: int) -> str:
    """Output folder of batch NN."""
    return f"batch_{nn:02d}"


def record_and_shard(mdir: Path, n: int, new_ids, shard_size: int, first_batch):
    header, rows = new_rows(mdir / "events.csv", new_ids)
    (mdir / f"events_new{n}.csv").write_bytes(header + b"".join(rows))
    (mdir / "rows.txt").write_bytes(b"".join(rows))
    shutil.copyfile(mdir / "events.csv", mdir / "events_ledger.csv")
    chunks = [rows[i:i + shard_size] for i in range(0, len(rows), shard_size)]
    for name, chunk in zip(shard_names(len(chunks)), chunks):
        (mdir / name).write_bytes(b"".join(chunk))

    written = []
    if first_batch is not None:
        for k, chunk in enumerate(chunks):
            nn = first_batch + k
            bdir = ROOT / batch_root(nn) / "manifests"
            bdir.mkdir(parents=True, exist_ok=True)
            ev = bdir / "events.csv"
            cycle_dir = bdir.resolve() == mdir.resolve()
            if ev.exists() and not cycle_dir:
                print(f"  b{nn:02d}: {ev.relative_to(ROOT)} already exists — left as is")
                continue
            # when the batch shares the cycle's folder, this subsets the cumulative
            # events.csv to just its new events; events_ledger.csv keeps the cumulative set
            ev.write_bytes(header + b"".join(chunk))
            written.append((nn, len(chunk)))
    return rows, chunks, written


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Stage A — flood event layer (one sampling cycle)")
    ap.add_argument("--root", required=True, help="output folder of this cycle, e.g. run_01")
    ap.add_argument("--n", type=int, required=True, help="events to sample this cycle")
    ap.add_argument("--append", action="store_true",
                    help="exclude the already-sampled set and the no-flood ledger")
    ap.add_argument("--seed", default=None,
                    help="cumulative sampled set to seed this cycle with (previous cycle's events_ledger.csv)")
    ap.add_argument("--shard-size", type=int, default=50)
    ap.add_argument("--first-batch", type=int, default=None,
                    help="batch number of the first shard; writes each batch's manifests/events.csv")
    a = ap.parse_args(argv)

    os.chdir(ROOT)                          # input catalogue paths are relative to the working folder
    cfg = make_cfg(a.root)
    cfg["sampling"]["n_events"] = a.n
    mdir = cfg.resolve_path("manifests")
    mdir.mkdir(parents=True, exist_ok=True)
    logs = ROOT / a.root / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(logs / "stageA.log")
    handler.setFormatter(logging.Formatter(_LOG_FORMAT))
    logging.getLogger().addHandler(handler)
    log_a = get_logger("flooded.stageA")

    if a.seed:
        target = mdir / "events.csv"
        if target.exists():
            log_a.info("%s already exists — not re-seeded", target.relative_to(ROOT))
        else:
            shutil.copyfile(a.seed, target)
            log_a.info("===== seeded with %d already-sampled events from %s =====",
                       sum(1 for _ in open(target)) - 1, a.seed)

    log_a.info("===== Stage A: sample %d %sevents (%s) =====", a.n, "NEW " if a.append else "", a.root)
    rc = sample_events(cfg, a.append)
    if rc != 0:
        return rc

    if a.append:
        new_ids = (mdir / "last_batch_events.txt").read_text().split()
    else:
        new_ids = [row_id(ln) for ln in raw_lines(mdir / "events.csv")[1:]]
    rows, chunks, written = record_and_shard(mdir, a.n, new_ids, a.shard_size, a.first_batch)
    log_a.info("===== recorded %d new events -> %d shard(s) of <= %d =====",
               len(rows), len(chunks), a.shard_size)
    for nn, k in written:
        log_a.info("  b%02d: %d events", nn, k)
    return 0


if __name__ == "__main__":
    sys.exit(main())
