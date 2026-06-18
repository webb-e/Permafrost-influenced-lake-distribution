"""
future_permafrost_analysis.py

Processes three projected permafrost warming scenarios alongside the present-day
IPA permafrost map. All data are aligned to the IPA's native North Pole Lambert
Azimuthal Equal Area CRS.

Outputs
-------
future_pf.tiff   — 3-band GeoTIFF (NP LAEA):
                   Band 1: 3degC classified (C=4, D=3, S=2, I=1, none=0)
                   Band 2: 2degC classified
                   Band 3: 1.5degC classified
future_PLD.gpkg  — Lakes with modernPF, 3degC, 2degC, 1.5degC attributes
"""

import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.warp import reproject as warp_reproject, Resampling, calculate_default_transform
from rasterio.crs import CRS
import geopandas as gpd
import xarray as xr
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────

DRIVE  = Path()
NC_DIR = DRIVE / estimated_permafrost_maps_NetCDF"  #Chadburn et al., 2017

IPA_SHP    = DRIVE / "permafrost_clean_reprojected.shp"  ##Brown et al., 2002
LAKES_GPKG = DRIVE / "PF_PLD/PF_PLD.gpkg"
OUT_DIR    = DRIVE / "future_PF"
OUT_RASTER = OUT_DIR / "future_pf.tiff"
OUT_LAKES  = OUT_DIR / "future_PLD.gpkg"

OUT_DIR.mkdir(parents=True, exist_ok=True)

# Scenarios: (column/band label, NetCDF filename)
# Labels must be valid SQL identifiers: no leading digits, no periods.
SCENARIOS = [
    ("pf_3degC",   NC_DIR / "map_3mean.nc"),
    ("pf_2degC",   NC_DIR / "map_2mean.nc"),
    ("pf_1p5degC", NC_DIR / "map_1.5mean.nc"),
]

# ── Constants ─────────────────────────────────────────────────────────────────

EXTENT_TO_INT = {"C": 4, "D": 3, "S": 2, "I": 1}  # no permafrost = 0
WGS84 = CRS.from_epsg(4326)

# ── Helpers ───────────────────────────────────────────────────────────────────

def classify_prob(arr: np.ndarray) -> np.ndarray:
    """
    Map continuous permafrost probability [0, 1] to integer class.
        [0.9, 1.0]  → 4 (Continuous)
        [0.5, 0.9)  → 3 (Discontinuous)
        [0.1, 0.5)  → 2 (Sporadic)
        [0.01, 0.1) → 1 (Isolated)
        < 0.01 / NaN → 0 (No permafrost)
    """
    safe = np.where(np.isnan(arr), 0.0, arr)
    out = np.zeros(safe.shape, dtype=np.int16)
    out = np.where((safe >= 0.9),                  4, out)
    out = np.where((safe >= 0.5)  & (safe < 0.9),  3, out)
    out = np.where((safe >= 0.1)  & (safe < 0.5),  2, out)
    out = np.where((safe >= 0.01) & (safe < 0.1),  1, out)
    return out


def read_nc_as_wgs84(path: Path):
    """Return (data float32 array N→S, src_transform, src_width, src_height)."""
    ds = xr.open_dataset(path)
    var_name = list(ds.data_vars)[0]
    da = ds[var_name].squeeze()

    def _coord(da, *names):
        for n in names:
            if n in da.coords:
                return da.coords[n].values
        return None

    lat_vals = _coord(da, "lat", "latitude", "y", "Y")
    lon_vals = _coord(da, "lon", "longitude", "x", "X")
    if lat_vals is None or lon_vals is None:
        raise ValueError(f"Cannot detect lat/lon in {path}. Coords: {list(da.coords)}")

    data = da.values.astype(np.float32)
    res_lon = float(lon_vals[1] - lon_vals[0])
    res_lat = float(lat_vals[1] - lat_vals[0])

    if res_lat < 0:
        north = float(lat_vals[0])  - res_lat / 2
    else:
        data = np.flipud(data)
        north = float(lat_vals[-1]) + res_lat / 2

    west = float(lon_vals[0]) - res_lon / 2
    h, w = data.shape
    transform = rasterio.transform.from_origin(west, north, abs(res_lon), abs(res_lat))
    return data, transform, w, h


def warp_to_target(data, src_transform, dst_transform, dst_w, dst_h, target_crs):
    """Reproject a float32 probability array from WGS84 to target_crs."""
    warped = np.full((dst_h, dst_w), np.nan, dtype=np.float32)
    warp_reproject(
        source=data,
        destination=warped,
        src_transform=src_transform,
        src_crs=WGS84,
        dst_transform=dst_transform,
        dst_crs=target_crs,
        resampling=Resampling.bilinear,
        src_nodata=np.nan,
        dst_nodata=np.nan,
    )
    return warped


def sample_band(band: np.ndarray, transform, points) -> np.ndarray:
    """Vectorised raster sampling — no Python loop over 3.9M points."""
    xs = np.fromiter((p.x for p in points), dtype=np.float64, count=len(points))
    ys = np.fromiter((p.y for p in points), dtype=np.float64, count=len(points))
    rows, cols = rasterio.transform.rowcol(transform, xs, ys)
    rows = np.asarray(rows, dtype=np.intp)
    cols = np.asarray(cols, dtype=np.intp)
    h, w = band.shape
    valid = (rows >= 0) & (rows < h) & (cols >= 0) & (cols < w)
    vals = np.zeros(len(rows), dtype=np.int16)
    vals[valid] = band[rows[valid], cols[valid]]
    return vals


# ── 1. Read IPA — establishes master CRS ──────────────────────────────────────

print("── Step 1: Reading IPA shapefile ──")
ipa = gpd.read_file(IPA_SHP)
TARGET_CRS = ipa.crs
print(f"  IPA CRS (master) : {TARGET_CRS}")
print(f"  EXTENT values    : {sorted(ipa['EXTENT'].dropna().unique())}")
print(f"  Feature count    : {len(ipa)}")

ipa["pf_int"] = (
    ipa["EXTENT"]
    .str.strip()
    .str.upper()
    .map(EXTENT_TO_INT)
    .fillna(0)
    .astype(int)
)
unmapped = (
    ipa.loc[(ipa["pf_int"] == 0) & ipa["EXTENT"].notna(), "EXTENT"]
    .unique()
    .tolist()
)
if unmapped:
    print(f"  WARNING — EXTENT values not mapped (assigned 0): {unmapped}")


# ── 2. Read first scenario to establish destination grid ──────────────────────
#    All 3 NetCDFs are assumed to share the same WGS84 grid.

print("\n── Step 2: Establishing destination grid from first scenario ──")
_, _tf0, _w0, _h0 = read_nc_as_wgs84(SCENARIOS[0][1])

dst_transform, dst_width, dst_height = calculate_default_transform(
    WGS84, TARGET_CRS, _w0, _h0,
    left=_tf0.c,
    bottom=_tf0.f + _h0 * _tf0.e,
    right=_tf0.c + _w0 * _tf0.a,
    top=_tf0.f,
)
print(f"  Destination grid : {dst_width} × {dst_height}")
print(f"  Destination transform : {dst_transform}")


# ── 3. Process each scenario ──────────────────────────────────────────────────

print("\n── Step 3: Warping and classifying scenarios ──")
scenario_bands = {}  # label → int16 classified array

for label, nc_path in SCENARIOS:
    print(f"  {label} : {nc_path.name}")
    data, src_tf, _, __ = read_nc_as_wgs84(nc_path)
    warped = warp_to_target(data, src_tf,
                            dst_transform, dst_width, dst_height, TARGET_CRS)
    classified = classify_prob(warped)
    print(f"    unique values: {np.unique(classified)}")
    scenario_bands[label] = classified


# ── 4. Rasterize IPA onto the same grid (for modernPF at lake centroids) ──────

print("\n── Step 4: Rasterizing IPA ──")
shapes_ipa = (
    (geom, val)
    for geom, val in zip(ipa.geometry, ipa["pf_int"])
    if geom is not None and not geom.is_empty
)
modern_class = rasterize(
    shapes_ipa,
    out_shape=(dst_height, dst_width),
    transform=dst_transform,
    fill=0,
    dtype=np.int16,
)
print(f"  Modern class unique values: {np.unique(modern_class)}")


# ── 5. Write 3-band raster ────────────────────────────────────────────────────

profile = {
    "driver": "GTiff",
    "dtype": "int16",
    "width": dst_width,
    "height": dst_height,
    "count": len(SCENARIOS),
    "crs": TARGET_CRS,
    "transform": dst_transform,
    "compress": "lzw",
    "nodata": -9999,
}

print(f"\n── Step 5: Writing raster → {OUT_RASTER} ──")
with rasterio.open(OUT_RASTER, "w", **profile) as dst:
    for band_idx, (label, _) in enumerate(SCENARIOS, start=1):
        dst.write(scenario_bands[label], band_idx)
        dst.update_tags(band_idx, description=f"{label} permafrost class (C=4,D=3,S=2,I=1,none=0)")
print("  Raster written.")


# ── 6. Sample raster at lake centroids ────────────────────────────────────────

print("\n── Step 6: Sampling raster at lake centroids ──")
lakes = gpd.read_file(LAKES_GPKG)
print(f"  Lakes CRS  : {lakes.crs}")
print(f"  Lake count : {len(lakes)}")

lakes_proj = lakes.to_crs(TARGET_CRS)
centroids  = lakes_proj.geometry.centroid

with rasterio.open(OUT_RASTER) as src:
    tf = src.transform
    bands = {label: src.read(i + 1) for i, (label, _) in enumerate(SCENARIOS)}

lakes["modernPF"] = sample_band(modern_class, dst_transform, centroids)
for label in [s[0] for s in SCENARIOS]:
    lakes[label] = sample_band(bands[label], tf, centroids)

print(f"  Sample of first 5 rows:\n{lakes[['modernPF'] + [s[0] for s in SCENARIOS]].head()}")


# ── 7. Save output GeoPackage ──────────────────────────────────────────────────

print(f"\n── Step 7: Saving lakes → {OUT_LAKES} ──")
lakes.to_file(OUT_LAKES, driver="GPKG")
print("  Done.\n")
print(f"Outputs:")
print(f"  {OUT_RASTER}")
print(f"  {OUT_LAKES}")
