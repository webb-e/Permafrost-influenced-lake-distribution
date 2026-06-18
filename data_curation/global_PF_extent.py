"""
Global PLD Glacier Flag + Permafrost EXTENT Assignment
=======================================================
Streams the SWOT PLD GDB in chunks via fiona (never loads all 6M lakes at
once) and distributes chunks across 2 worker processes.

Workflow per chunk
------------------
1. Skip lake_ids already in PLD_PF output.
2. Split remaining lakes into northern (lat >= 0) and southern (lat < 0).
3. Glacier flag:
     - Northern lakes         --> sjoin with northern glacier file (EPSG:3575)
     - Southern lakes         --> sjoin with southern glacier file (EPSG:4326)
     - Already-processed lakes --> on_glacier = 0 (confirmed)
4. Permafrost EXTENT:
     - Northern lakes         --> sjoin with PF file in EPSG:3575 (largest-
                                  intersection rule for multi-polygon matches)
     - Southern lakes         --> EXTENT = 'None'
     - Already-processed lakes --> EXTENT carried over from existing file
5. Append [lake_id, EXTENT, on_glacier] to output CSV as each chunk completes.

Output
------
  PLD_PF_global.csv
"""

import os
import gc
import csv
import warnings
import threading
import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
from shapely.geometry import shape
from shapely.validation import make_valid
from concurrent.futures import ProcessPoolExecutor, as_completed

warnings.filterwarnings('ignore', 'GeoSeries.notna', UserWarning)

# ---------------------------------------------------------------------------
# File paths -- edit these
# ---------------------------------------------------------------------------
lakes_gdb_path      = r'SWOT_PLD_v201_02042025_attributes_updated.gdb' ## Wang et al., 2025
lakes_gdb_layer     = 'SWOT_PLD_v201_02042025_attributes_updated'

existing_pf_path    = r'PLD_PF.gpkg'
existing_pf_layer   = 'PLD_PF'

permafrost_fp       = r'/permafrost_clean_reprojected.shp' ## Brown et al., 2002
glaciers_north_fp   = r'glaciers_reprojected.shp'## GLIMS and NSIDC 2026
glaciers_south_fp   = r'glims_polygons.shp' ## GLIMS and NSIDC 2026

output_dir          = r''
output_csv          = os.path.join(output_dir, "PLD_PF_global.csv")

NORTH_CRS   = "EPSG:3575"
CHUNK_SIZE  = 200_000
N_WORKERS   = 2

# ---------------------------------------------------------------------------
# Helpers (used in both main and worker processes)
# ---------------------------------------------------------------------------
def fix_geometries(gdf):
    n_null = gdf.geometry.isna().sum()
    if n_null > 0:
        gdf = gdf[gdf.geometry.notna()].copy()
    invalid = ~gdf.geometry.is_valid
    if invalid.any():
        gdf.loc[invalid, 'geometry'] = gdf.loc[invalid, 'geometry'].apply(
            lambda g: g.buffer(0) if g is not None else g
        )
    return gdf


def to_crs(gdf, crs):
    if gdf.crs is None:
        raise ValueError("GeoDataFrame has no CRS set.")
    if gdf.crs.to_epsg() != int(crs.split(':')[1]):
        return gdf.to_crs(crs)
    return gdf


# ---------------------------------------------------------------------------
# Worker initializer
# Runs once per worker process at startup; stores overlay layers as globals
# so they are never re-serialised for each chunk.
# ---------------------------------------------------------------------------
_pf         = None   # permafrost GeoDataFrame (worker-local)
_gl_north   = None   # northern glaciers (worker-local)
_gl_south   = None   # southern glaciers (worker-local)


def worker_init(permafrost_fp, glaciers_north_fp, glaciers_south_fp):
    global _pf, _gl_north, _gl_south

    _pf = gpd.read_file(permafrost_fp)
    _pf = _pf[_pf['EXTENT'].isin(['C', 'D', 'I', 'S'])].copy()
    _pf = fix_geometries(_pf)
    _pf = _pf.reset_index(drop=True)
    _pf['_pf_idx'] = _pf.index
    # Pre-project to NORTH_CRS once
    _pf = to_crs(_pf, NORTH_CRS)

    _gl_north = gpd.read_file(glaciers_north_fp)
    _gl_north = fix_geometries(_gl_north)
    _gl_north = to_crs(_gl_north, NORTH_CRS)

    _gl_south = gpd.read_file(glaciers_south_fp)
    _gl_south = fix_geometries(_gl_south)
    # southern file is already EPSG:4326; keep as-is


# ---------------------------------------------------------------------------
# Per-chunk processing (runs in worker process)
# ---------------------------------------------------------------------------
def process_chunk(features, src_crs, existing_ids):
    """
    Parameters
    ----------
    features    : list of fiona feature dicts
    src_crs     : CRS string of the source GDB
    existing_ids: frozenset of lake_ids already processed

    Returns
    -------
    DataFrame with columns [lake_id, EXTENT, on_glacier, lat, poly_area]
    """
    # Build GeoDataFrame from raw fiona features
    rows = []
    for feat in features:
        props = dict(feat['properties'])
        geom  = feat['geometry']
        props['geometry'] = shape(geom) if geom else None
        rows.append(props)

    gdf = gpd.GeoDataFrame(rows, geometry='geometry', crs=src_crs)
    gdf = fix_geometries(gdf)

    # ---- Split into already-processed / north / south ----
    mask_existing = gdf['lake_id'].isin(existing_ids)
    existing_chunk = gdf[mask_existing][['lake_id', 'lat', 'poly_area']].copy()
    remaining      = gdf[~mask_existing].copy()
    del gdf
    gc.collect()

    north = remaining[remaining['lat'] >= 0].copy().reset_index(drop=True)
    south = remaining[remaining['lat'] <  0].copy().reset_index(drop=True)
    del remaining
    gc.collect()

    results = []

    # ---- Already-processed: EXTENT from existing file, on_glacier = 0 ----
    if not existing_chunk.empty:
        existing_chunk['EXTENT']     = 'already_processed'  # placeholder; merged later
        existing_chunk['on_glacier'] = np.int8(0)
        results.append(existing_chunk[['lake_id', 'EXTENT', 'on_glacier', 'lat', 'poly_area']])

    # ---- Northern lakes ----
    if not north.empty:
        north_proj = to_crs(north[['lake_id', 'geometry']].copy(), NORTH_CRS)

        # -- Glacier flag --
        gl_join = gpd.sjoin(
            north_proj,
            _gl_north[['geometry']],
            how='left',
            predicate='intersects'
        )
        glacier_ids = set(gl_join.loc[gl_join['index_right'].notna(), 'lake_id'].unique())

        # -- EXTENT (largest intersection) --
        pf_join = gpd.sjoin(
            north_proj,
            _pf[['EXTENT', '_pf_idx', 'geometry']],
            how='left',
            predicate='intersects'
        )

        dup_mask = pf_join.duplicated(subset='lake_id', keep=False)
        singles  = pf_join[~dup_mask][['lake_id', 'EXTENT']].copy()

        if dup_mask.any():
            multi = pf_join[dup_mask].copy().reset_index(drop=True)
            pf_geoms = _pf.set_index('_pf_idx')['geometry']
            pf_geoms = pf_geoms[~pf_geoms.index.duplicated(keep='first')]
            multi['_pf_geom'] = multi['index_right'].map(pf_geoms)
            multi['_area'] = multi.apply(
                lambda r: r.geometry.intersection(r['_pf_geom']).area
                if r['_pf_geom'] is not None else 0.0,
                axis=1
            )
            best = multi.loc[
                multi.groupby('lake_id')['_area'].idxmax(), ['lake_id', 'EXTENT']
            ]
            extent_df = pd.concat([singles, best], ignore_index=True)
        else:
            extent_df = singles

        extent_map = extent_df.set_index('lake_id')['EXTENT']

        df_north = pd.DataFrame({'lake_id': north['lake_id'].values})
        df_north['EXTENT']     = df_north['lake_id'].map(extent_map).fillna('None')
        df_north['on_glacier'] = north['lake_id'].isin(glacier_ids).astype(np.int8).values
        df_north['lat']        = north['lat'].values
        df_north['poly_area']  = north['poly_area'].values
        results.append(df_north)
        del north_proj, pf_join, gl_join
        gc.collect()

    # ---- Southern lakes ----
    if not south.empty:
        gl_join = gpd.sjoin(
            south[['lake_id', 'geometry']],
            _gl_south[['geometry']],
            how='left',
            predicate='intersects'
        )
        glacier_ids_s = set(gl_join.loc[gl_join['index_right'].notna(), 'lake_id'].unique())

        df_south = pd.DataFrame({'lake_id': south['lake_id'].values})
        df_south['EXTENT']     = 'None'
        df_south['on_glacier'] = south['lake_id'].isin(glacier_ids_s).astype(np.int8).values
        df_south['lat']        = south['lat'].values
        df_south['poly_area']  = south['poly_area'].values
        results.append(df_south)

    if not results:
        return pd.DataFrame(columns=['lake_id', 'EXTENT', 'on_glacier', 'lat', 'poly_area'])

    return pd.concat(results, ignore_index=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    os.makedirs(output_dir, exist_ok=True)

    # -----------------------------------------------------------------------
    # Load existing lake_ids + their EXTENT values
    # -----------------------------------------------------------------------
    print("Loading already-processed lake_ids...")
    ext = os.path.splitext(existing_pf_path)[1].lower()
    if ext == '.csv':
        existing_df = pd.read_csv(existing_pf_path, usecols=['lake_id', 'EXTENT'])
    else:
        existing_df = gpd.read_file(existing_pf_path, layer=existing_pf_layer)[['lake_id', 'EXTENT']]
        existing_df = pd.DataFrame(existing_df.drop(columns='geometry', errors='ignore'))

    existing_ids      = frozenset(existing_df['lake_id'].unique())
    existing_extent   = existing_df.set_index('lake_id')['EXTENT'].to_dict()
    print(f"  Already-processed lakes: {len(existing_ids):,}")
    del existing_df

    # -----------------------------------------------------------------------
    # Initialize output CSV
    # -----------------------------------------------------------------------
    write_lock = threading.Lock()

    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['lake_id', 'EXTENT', 'on_glacier', 'lat', 'poly_area'])

    def append_to_csv(df):
        """Thread-safe append of a DataFrame to the output CSV."""
        if df.empty:
            return
        # Resolve placeholder EXTENT for already-processed lakes
        mask = df['EXTENT'] == 'already_processed'
        if mask.any():
            df.loc[mask, 'EXTENT'] = df.loc[mask, 'lake_id'].map(existing_extent)
        with write_lock:
            df.to_csv(output_csv, mode='a', header=False, index=False)

    # -----------------------------------------------------------------------
    # Stream GDB and dispatch chunks to workers
    # -----------------------------------------------------------------------
    print("Streaming GDB and dispatching to workers...")

    with fiona.open(lakes_gdb_path, layer=lakes_gdb_layer) as src:
        src_crs    = src.crs_wkt
        total_feat = len(src)
        print(f"  Total features in GDB: {total_feat:,}")

        futures = {}
        chunk   = []
        chunk_n = 0

        with ProcessPoolExecutor(
            max_workers=N_WORKERS,
            initializer=worker_init,
            initargs=(permafrost_fp, glaciers_north_fp, glaciers_south_fp)
        ) as executor:

            for i, feature in enumerate(src):
                chunk.append(feature)

                if len(chunk) == CHUNK_SIZE:
                    chunk_n += 1
                    fut = executor.submit(
                        process_chunk, chunk, src_crs, existing_ids
                    )
                    futures[fut] = chunk_n
                    chunk = []

                # Drain completed futures periodically to free memory
                if len(futures) >= N_WORKERS * 2:
                    done = [f for f in futures if f.done()]
                    for f in done:
                        n = futures.pop(f)
                        result = f.result()
                        append_to_csv(result)
                        print(f"  Chunk {n} complete ({len(result):,} rows written)")

            # Submit final partial chunk
            if chunk:
                chunk_n += 1
                fut = executor.submit(process_chunk, chunk, src_crs, existing_ids)
                futures[fut] = chunk_n

            # Collect remaining futures
            for f in as_completed(futures):
                n = futures[f]
                result = f.result()
                append_to_csv(result)
                print(f"  Chunk {n} complete ({len(result):,} rows written)")

    # -----------------------------------------------------------------------
    # Deduplicate (a lake_id could appear in both existing and remaining
    # if the existing file was loaded from CSV and IDs overlap at chunk edges)
    # -----------------------------------------------------------------------
    print("Deduplicating output...")
    final = pd.read_csv(output_csv)
    n_before = len(final)
    final = final.drop_duplicates('lake_id').reset_index(drop=True)
    print(f"  Rows before dedup: {n_before:,}  after: {len(final):,}")
    final.to_csv(output_csv, index=False)

    print(f"\nDone.")
    print(f"  Total lakes : {len(final):,}")
    print(f"  EXTENT breakdown:\n{final['EXTENT'].value_counts(dropna=False).to_string()}")
    print(f"  on_glacier=1: {(final['on_glacier'] == 1).sum():,}")
    print(f"  Output      : {output_csv}")


if __name__ == "__main__":
    main()
