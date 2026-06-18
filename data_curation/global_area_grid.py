import geopandas as gpd
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from shapely.geometry import box

# ── Paths ────────────────────────────────────────────────────────────────────
GRID_PATH       = "Grid_global_50km.shp"
GLACIERS_S_PATH = "glims_polygons.shp" ##GLIMS and NSIDC 2026
GLACIERS_N_PATH = "glims_polygons.shp" ##GLIMS and NSIDC 2026
OUTPUT_PATH     = "Grid_global_50km_non_glacier_area.gpkg"

N_WORKERS  = max(1, multiprocessing.cpu_count() - 1)
CHUNK_SIZE = 500


# ── Worker function ───────────────────────────────────────────────────────────
def process_chunk(args):
    chunk_gdf, land_gdf, glaciers_gdf, chunk_idx = args

    results = []
    for _, cell in chunk_gdf.iterrows():
        cell_geom = cell.geometry
        cell_id   = cell["cell_id"]
        row = {"cell_id": cell_id, "land_area_m2": 0.0, "glacier_area_m2": 0.0}

        # Land area
        try:
            land_candidates = land_gdf[land_gdf.intersects(cell_geom)]
            if not land_candidates.empty:
                row["land_area_m2"] = land_candidates.intersection(cell_geom).area.sum()
        except Exception as e:
            print(f"  Warning: land intersection failed for cell {cell_id}: {e}")

        # Glacier area
        try:
            glac_candidates = glaciers_gdf[glaciers_gdf.intersects(cell_geom)]
            if not glac_candidates.empty:
                row["glacier_area_m2"] = glac_candidates.intersection(cell_geom).area.sum()
        except Exception as e:
            print(f"  Warning: glacier intersection failed for cell {cell_id}: {e}")

        results.append(row)

    print(f"  Chunk {chunk_idx} done ({len(chunk_gdf)} cells)")
    return pd.DataFrame(results)


# ── Spatial filter helper ─────────────────────────────────────────────────────
def spatially_filter(source_gdf, target_gdf):
    bbox_geom = box(*target_gdf.total_bounds)
    return source_gdf[source_gdf.intersects(bbox_geom)].copy()


# ── Geometry repair helper ────────────────────────────────────────────────────
def repair_geometries(gdf, name):
    print(f"  Repairing {name} geometries with buffer(0)...")
    gdf["geometry"] = gdf.geometry.buffer(0)

    invalid = ~gdf.geometry.is_valid
    if invalid.any():
        print(f"  Dropping {invalid.sum()} still-invalid {name} geometries")
        gdf = gdf[gdf.geometry.is_valid].copy()
    else:
        print(f"  All {name} geometries valid")

    return gdf


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":

    # 1. Load grid first to get its CRS
    print("Loading grid...")
    grid = gpd.read_file(GRID_PATH)
    TARGET_CRS = grid.crs
    print(f"  Grid CRS: {TARGET_CRS}")
    print(f"  Grid cells: {len(grid):,}")

    if "cell_id" not in grid.columns:
        grid["cell_id"] = range(len(grid))

    # 2. Load glaciers
    print("\nLoading glaciers...")
    glaciers = pd.concat([
        gpd.read_file(GLACIERS_S_PATH),
        gpd.read_file(GLACIERS_N_PATH)
    ], ignore_index=True)
    print(f"  Total glacier polygons: {len(glaciers):,}")

    # 3. Load land mask
    print("Loading  land mask...")
    land = gpd.read_file("/Users/elizabethwebb/Downloads/0304143/1.1/data/0-data/GSHHS_shp/f/GSHHS_f_L1.shp")
    # 4. Reproject to match grid CRS
    print(f"\nReprojecting to match grid CRS...")
    glaciers = glaciers.to_crs(TARGET_CRS)
    land     = land.to_crs(TARGET_CRS)

    # 5. Repair geometries
    print("\nRepairing geometries...")
    glaciers = repair_geometries(glaciers, "glaciers")
    land     = repair_geometries(land, "land")

    # 6. Dissolve glaciers to remove overlapping polygons, then explode to single-part
    print("\nDissolving glaciers (may take a few minutes)...")
    glaciers = glaciers.dissolve().reset_index(drop=True)
    glaciers = glaciers.explode(index_parts=False).reset_index(drop=True)
    # Repair again after dissolve as it can introduce new invalid geometries
    glaciers = repair_geometries(glaciers, "glaciers (post-dissolve)")
    print(f"  Glacier parts after dissolve+explode: {len(glaciers):,}")

    land = land.explode(index_parts=False).reset_index(drop=True)

    # 7. Split grid into chunks
    cell_ids = grid["cell_id"].tolist()
    chunks = [
        grid[grid["cell_id"].isin(cell_ids[i:i + CHUNK_SIZE])].copy()
        for i in range(0, len(cell_ids), CHUNK_SIZE)
    ]
    print(f"\nDispatching {len(grid):,} cells across {len(chunks)} chunks "
          f"using {N_WORKERS} workers...")

    # 8. Pre-filter land/glaciers to each chunk's bounding box before dispatch
    args_list = [
        (chunk, spatially_filter(land, chunk), spatially_filter(glaciers, chunk), i)
        for i, chunk in enumerate(chunks)
    ]

    # 9. Parallel execution
    all_results = []
    with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
        futures = {executor.submit(process_chunk, args): args[3] for args in args_list}
        for future in as_completed(futures):
            try:
                all_results.append(future.result())
            except Exception as e:
                print(f"  Chunk {futures[future]} failed: {e}")

    # 10. Combine and calculate final areas
    print("\nCombining results...")
    results_df = pd.concat(all_results, ignore_index=True)

    grid = grid.merge(results_df, on="cell_id", how="left")
    grid["land_area_m2"]    = grid["land_area_m2"].fillna(0)
    grid["glacier_area_m2"] = grid["glacier_area_m2"].fillna(0)

    # Clip glacier area to land area (guard against floating point artifacts)
    grid["glacier_area_m2"]           = np.minimum(grid["glacier_area_m2"], grid["land_area_m2"])
    grid["non_glacier_land_area_m2"]  = grid["land_area_m2"] - grid["glacier_area_m2"]
    grid["land_area_km2"]             = grid["land_area_m2"] / 1e6
    grid["glacier_area_km2"]          = grid["glacier_area_m2"] / 1e6
    grid["non_glacier_land_area_km2"] = grid["non_glacier_land_area_m2"] / 1e6

    # 11. Save
    print(f"\nSaving to {OUTPUT_PATH}...")
    grid.to_file(OUTPUT_PATH, driver="GPKG")

    print("\n── Summary ──────────────────────────────────────────────────────────")
    print(f"  Total grid cells:             {len(grid):,}")
    print(f"  Cells with any land:          {(grid['land_area_m2'] > 0).sum():,}")
    print(f"  Cells with any glacier:       {(grid['glacier_area_m2'] > 0).sum():,}")
    print(f"  Total land area:              {grid['land_area_km2'].sum():,.0f} km²")
    print(f"  Total glacier area:           {grid['glacier_area_km2'].sum():,.0f} km²")
    print(f"  Total non-glacier land area:  {grid['non_glacier_land_area_km2'].sum():,.0f} km²")
