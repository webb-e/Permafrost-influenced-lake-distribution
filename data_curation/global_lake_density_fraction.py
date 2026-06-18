#!/usr/bin/env python3
"""
Global Lake Density and Fraction Calculator
============================================
For each 50 km grid cell (MODIS sinusoidal projection), computes:

  lake_area_km2            : lake polygon area intersecting the cell (km²)
  lake_count               : number of lake centroids within the cell
  lake_fraction            : lake_area_km2 / non_glacier_land_area_km2
  lake_density_per100km2   : (lake_count / non_glacier_land_area_km2) * 100

Lakes marked on_glacier = 1 (via PLD_PF_global.csv join) are excluded.

Method
------
  Follows the same single-process sjoin → overlay → intersect pattern as
  calculate_density_fraction.py.  No multiprocessing.

  Because the dataset is global (~6M lakes), the script works in batches
  of grid rows to avoid materialising one huge intersection at once.  Each
  batch result is written to a CSV checkpoint immediately so the run can be
  resumed after a crash.

Output
------
  <output_dir>/chunks/                      -- per-batch CSV checkpoints
  <output_dir>/global_lake_density_fraction.csv
  <output_dir>/global_lake_density_fraction.gpkg
"""

import os
import gc
import glob
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.validation import make_valid

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
GRID_FP = ("Grid_global_50km_non_glacier_area.gpkg")

LAKES_GDB_PATH  = r"/SWOT_PLD_v201_02042025_attributes_updated.gdb"
LAKES_GDB_LAYER = "SWOT_PLD_v201_02042025_attributes_updated"

ON_GLACIER_CSV = (r"Lake distribution/PF_PLD/PLD_PF_global.csv")

OUTPUT_DIR = (/global_lake_density_fraction")

LAND_AREA_COL  = "non_glacier_land_area_km2"
MIN_LAND_AREA  = 1.0      # km²  — skip cells smaller than this
BATCH_SIZE     = 2_000    # grid cells per checkpoint batch

# ---------------------------------------------------------------------------
# Helpers  (same pattern as example script)
# ---------------------------------------------------------------------------

def ensure_crs(gdf, target):
    if gdf.crs is None:
        return gdf.set_crs(target, allow_override=True)
    if gdf.crs.to_string() == target:
        return gdf
    return gdf.to_crs(target)


def fix_geoms(gdf):
    n_null = gdf.geometry.isna().sum()
    if n_null:
        print(f"    Dropping {n_null} null geometries")
        gdf = gdf[gdf.geometry.notna()].copy()
    invalid = ~gdf.geometry.is_valid
    if invalid.any():
        print(f"    Fixing {invalid.sum()} invalid geometries")
        gdf.loc[invalid, "geometry"] = gdf.loc[invalid, "geometry"].apply(
            lambda g: make_valid(g) if g is not None else g
        )
    return gdf


def compute_metrics_for_batch(grid_batch, lakes, centroids_gdf):
    """
    Given a subset of grid cells, compute lake_area_km2 and lake_count
    using the same sjoin → per-row intersection pattern as the example.

    Returns a DataFrame with columns:
        grid_id, lake_area_km2, lake_count
    """
    grid_batch = grid_batch.reset_index(drop=True).copy()
    grid_batch["_inter_idx"] = grid_batch.index
    inter_geom_dict = grid_batch.set_index("_inter_idx")["geometry"].to_dict()
    inter_grid_map  = grid_batch.set_index("_inter_idx")["grid_id"].to_dict()
    inter_bbox      = grid_batch.total_bounds

    # Initialise accumulators
    lakearea_by_grid  = {}
    lakecount_by_grid = {}

    # ---- Lake area (polygon intersection) ----
    lake_cands_idx = list(lakes.sindex.intersection(inter_bbox))
    if lake_cands_idx:
        lakes_sub = lakes.iloc[lake_cands_idx].copy()
        joined = gpd.sjoin(
            lakes_sub[["lake_idx", "geometry"]],
            grid_batch[["_inter_idx", "geometry"]],
            how="inner",
            predicate="intersects",
        )
        if not joined.empty:
            for _, row in joined.iterrows():
                iidx    = row["_inter_idx"]
                lake_g  = row["geometry"]
                inter_g = inter_geom_dict.get(iidx)
                if inter_g is None or inter_g.is_empty or lake_g is None or lake_g.is_empty:
                    continue
                try:
                    clipped = lake_g.intersection(inter_g)
                    if not clipped.is_empty:
                        gid = inter_grid_map[iidx]
                        lakearea_by_grid[gid] = lakearea_by_grid.get(gid, 0.0) + clipped.area
                except Exception:
                    continue

    # ---- Lake count (centroid within cell) ----
    centroid_cands_idx = list(centroids_gdf.sindex.intersection(inter_bbox))
    if centroid_cands_idx:
        centroids_sub = centroids_gdf.iloc[centroid_cands_idx].copy()
        c_joined = gpd.sjoin(
            centroids_sub[["lake_idx", "geometry"]],
            grid_batch[["_inter_idx", "geometry"]],
            how="inner",
            predicate="within",
        )
        if not c_joined.empty:
            counts = c_joined.groupby("_inter_idx")["lake_idx"].nunique().to_dict()
            for iidx, cnt in counts.items():
                gid = inter_grid_map.get(iidx)
                if gid is not None:
                    lakecount_by_grid[gid] = lakecount_by_grid.get(gid, 0) + int(cnt)

    # Assemble results
    all_gids = grid_batch["grid_id"].tolist()
    rows = []
    for gid in all_gids:
        rows.append({
            "grid_id":      gid,
            "lake_area_km2": lakearea_by_grid.get(gid, 0.0) / 1e6,
            "lake_count":    lakecount_by_grid.get(gid, 0),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    chunk_dir = os.path.join(OUTPUT_DIR, "chunks")
    os.makedirs(chunk_dir, exist_ok=True)

    # -----------------------------------------------------------------------
    print("=" * 60)
    print("Loading grid...")
    grid = gpd.read_file(GRID_FP)
    TARGET_CRS = grid.crs.to_string() if grid.crs else "ESRI:54008"
    print(f"  CRS: {TARGET_CRS}")
    grid = fix_geoms(grid)
    grid = grid.reset_index(drop=True)
    if "grid_id" not in grid.columns:
        grid["grid_id"] = grid.index.astype(int)
    if LAND_AREA_COL not in grid.columns:
        raise ValueError(
            f"Column '{LAND_AREA_COL}' not found in grid.\n"
            f"Available: {list(grid.columns)}"
        )
    grid = grid[grid[LAND_AREA_COL].fillna(0) >= MIN_LAND_AREA].copy()
    grid = grid.reset_index(drop=True)
    print(f"  Grid cells (after land-area filter): {len(grid):,}")

    # -----------------------------------------------------------------------
    print("=" * 60)
    print("Loading lakes...")
    lakes = gpd.read_file(LAKES_GDB_PATH, layer=LAKES_GDB_LAYER)
    print(f"  Raw lakes: {len(lakes):,}")

    # Normalise lake_id column name
    lid_col = next((c for c in lakes.columns if c.lower() == "lake_id"), None)
    if lid_col is None:
        raise ValueError(f"No 'lake_id' column. Columns: {list(lakes.columns)}")
    if lid_col != "lake_id":
        lakes = lakes.rename(columns={lid_col: "lake_id"})

    print("  Joining on-glacier flag...")
    og = pd.read_csv(
        ON_GLACIER_CSV,
        usecols=["lake_id", "on_glacier"],
        dtype={"lake_id": str, "on_glacier": "Int8"},
    )
    og["lake_id"]    = og["lake_id"].astype(str).str.strip()
    lakes["lake_id"] = lakes["lake_id"].astype(str).str.strip()
    lakes = lakes.merge(og, on="lake_id", how="left")
    lakes["on_glacier"] = lakes["on_glacier"].fillna(0)
    n_glacier = (lakes["on_glacier"] == 1).sum()
    lakes = lakes[lakes["on_glacier"] != 1].copy()
    print(f"  Removed {n_glacier:,} on-glacier lakes → {len(lakes):,} remaining")

    lakes = fix_geoms(lakes)
    lakes = lakes.reset_index(drop=True)
    lakes["lake_idx"] = lakes.index
    lakes = ensure_crs(lakes, TARGET_CRS)
    # Keep only what we need
    lakes = lakes[["lake_idx", "geometry"]].copy()

    # -----------------------------------------------------------------------
    print("=" * 60)
    print("Building centroids...")
    # Re-read to get lon/lat columns if present
    lakes_full = gpd.read_file(LAKES_GDB_PATH, layer=LAKES_GDB_LAYER)
    lakes_full["lake_id"] = lakes_full[
        next((c for c in lakes_full.columns if c.lower() == "lake_id"), lakes_full.columns[0])
    ].astype(str).str.strip()
    # Filter to same set as lakes
    lakes_full = lakes_full[lakes_full["lake_id"].isin(
        pd.read_csv(ON_GLACIER_CSV, usecols=["lake_id","on_glacier"],
                    dtype={"lake_id": str, "on_glacier": "Int8"})
        .query("on_glacier != 1")["lake_id"].astype(str).str.strip()
        # union with unmatched (treated as not on glacier)
    ) | ~lakes_full["lake_id"].isin(
        pd.read_csv(ON_GLACIER_CSV, usecols=["lake_id"],
                    dtype={"lake_id": str})["lake_id"].astype(str).str.strip()
    )].copy()
    lakes_full = lakes_full[lakes_full.geometry.notna() & ~lakes_full.geometry.is_empty].copy()
    lakes_full = lakes_full.reset_index(drop=True)
    lakes_full["lake_idx"] = lakes_full.index

    lon_col = next((c for c in lakes_full.columns if c.lower() in ("lon", "longitude", "x_lon")), None)
    lat_col = next((c for c in lakes_full.columns if c.lower() in ("lat", "latitude", "y_lat")), None)

    if lon_col and lat_col:
        print(f"  Using '{lon_col}'/'{lat_col}' for centroids")
        centroids_gdf = gpd.GeoDataFrame(
            {"lake_idx": lakes_full["lake_idx"]},
            geometry=gpd.points_from_xy(
                lakes_full[lon_col].astype(float),
                lakes_full[lat_col].astype(float),
            ),
            crs="EPSG:4326",
        ).to_crs(TARGET_CRS)
    else:
        print("  Computing centroids from geometry...")
        centroids_gdf = gpd.GeoDataFrame(
            {"lake_idx": lakes_full["lake_idx"]},
            geometry=ensure_crs(lakes_full, TARGET_CRS).geometry.centroid,
            crs=TARGET_CRS,
        )
    centroids_gdf = centroids_gdf[
        centroids_gdf.geometry.notna() & ~centroids_gdf.geometry.is_empty
    ].copy()
    print(f"  Centroids: {len(centroids_gdf):,}")
    del lakes_full
    gc.collect()

    # -----------------------------------------------------------------------
    print("=" * 60)
    print("Processing grid in batches (checkpointing each batch)...")

    all_gids    = grid["grid_id"].tolist()
    n_batches   = (len(all_gids) + BATCH_SIZE - 1) // BATCH_SIZE
    print(f"  Total cells: {len(all_gids):,}  |  Batch size: {BATCH_SIZE}  |  Batches: {n_batches}")

    for batch_i in range(n_batches):
        out_csv = os.path.join(chunk_dir, f"batch_{batch_i:05d}.csv")
        if os.path.exists(out_csv):
            print(f"  [{batch_i+1:>4}/{n_batches}]  skipped (already done)")
            continue

        start = batch_i * BATCH_SIZE
        end   = min(start + BATCH_SIZE, len(all_gids))
        batch_ids    = all_gids[start:end]
        grid_batch   = grid[grid["grid_id"].isin(batch_ids)][["grid_id", "geometry"]].copy()

        results = compute_metrics_for_batch(grid_batch, lakes, centroids_gdf)
        results.to_csv(out_csv, index=False)
        print(f"  [{batch_i+1:>4}/{n_batches}]  rows {start}–{end-1}  saved → {os.path.basename(out_csv)}")
        gc.collect()

    # -----------------------------------------------------------------------
    print("=" * 60)
    print("Assembling all batches...")
    chunk_files = sorted(glob.glob(os.path.join(chunk_dir, "batch_*.csv")))
    results_df  = pd.concat([pd.read_csv(f) for f in chunk_files], ignore_index=True)
    print(f"  Total result rows: {len(results_df):,}")

    # Compute fraction and density
    grid_out = grid.merge(results_df, on="grid_id", how="left")
    denom    = grid_out[LAND_AREA_COL].fillna(0)
    valid    = denom > 0
    grid_out["lake_fraction"]         = np.nan
    grid_out["lake_density_per100km2"] = np.nan
    grid_out.loc[valid, "lake_fraction"]          = grid_out.loc[valid, "lake_area_km2"].fillna(0) / denom[valid]
    grid_out.loc[valid, "lake_density_per100km2"] = (grid_out.loc[valid, "lake_count"].fillna(0) / denom[valid]) * 100

    # -----------------------------------------------------------------------
    print("=" * 60)
    print("Saving outputs...")

    out_csv_final = os.path.join(OUTPUT_DIR, "global_lake_density_fraction.csv")
    print(f"  CSV  → {out_csv_final}")
    grid_out.drop(columns=["geometry"]).to_csv(out_csv_final, index=False)

    out_gpkg = os.path.join(OUTPUT_DIR, "global_lake_density_fraction.gpkg")
    print(f"  GPKG → {out_gpkg}")
    if os.path.exists(out_gpkg):
        os.remove(out_gpkg)
    grid_out.to_file(out_gpkg, layer="global_lake_metrics", driver="GPKG")

    # -----------------------------------------------------------------------
    print("=" * 60)
    print("Summary")
    valid_f  = grid_out["lake_fraction"].notna()
    print(f"  Cells with valid lake fraction : {valid_f.sum():,}")
    print(f"  Max lake fraction              : {grid_out['lake_fraction'].max():.4f}")
    nonzero  = grid_out.loc[grid_out["lake_fraction"] > 0, "lake_fraction"]
    print(f"  Mean lake fraction (non-zero)  : {nonzero.mean():.4f}")
    print(f"  Total lake area (km²)          : {grid_out['lake_area_km2'].sum():,.1f}")
    print(f"  Total lake count               : {grid_out['lake_count'].sum():,.0f}")
    over1 = (grid_out["lake_fraction"] > 1).sum()
    if over1:
        print(f"  *** WARNING: {over1} cells have lake_fraction > 1")
    print("\nDone.")


if __name__ == "__main__":
    main()
