"""
Northern Permafrost Lake Filter and Attribute Assignment
=========================================================
Filters the SWOT PLD lake dataset to lakes that are:
  1. North of 50 deg N
  2. Intersect permafrost zones (EXTENT = C, D, I, S)
  3. Do NOT intersect the glacier shapefile

Then assigns the following attributes to each qualifying lake (largest
intersection area rule where a lake spans multiple polygons):
  - EXTENT       : Permafrost extent class (C/D/I/S)
  - CONTENT      : Ground ice content (Low/Medium/High)
  - glaciated    : Whether lake was glaciated at LGM (glaciated/unglaciated)
  - yedoma       : Yedoma domain (yedoma / non-yedoma)
  - thermokarst_lakes    : TkThLP thermokarst lake/pond potential
  - thermokarst_wetlands : TKWP thermokarst wetland potential
  - biome        : WWF biome name
  - shield       : Canadian Shield (True/False)

------------------------------
All spatial filtering (permafrost intersection, glacier exclusion) is performed
in EPSG:3575 (North Pole LAEA).

Output
------
  PLD_PF.gpkg  -- filtered + attributed lakes (GeoPackage)
  PLD_PF.csv   -- same, without geometry
"""

import os
import gc
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.validation import make_valid

warnings.filterwarnings('ignore', 'GeoSeries.notna', UserWarning)

# ---------------------------------------------------------------------------
# File paths -- edit these
# ---------------------------------------------------------------------------
lakes_gdb_path      = r'SWOT_PLD_v201_02042025_attributes_updated.gdb' ## Wang et al., 2025
permafrost_fp       = r'permafrost_clean_reprojected.shp' ## Brown et al., 2002
glacial_history_fp  = r'LGM_best_estimate.shp' ## Batchelor et al., 2025
glaciers_fp         = r'glaciers_reprojected.shp' ## GLIMS and NSIDC 2026
yedoma_fp           = r"yedoma_domain_reprojected.shp" ## Strauss et al 2021
thermokarst_fp      = r"Circumpolar_Thermokarst_Landscapes/Thermokarst_reprojected.shp" ## Olefeldt et al 2016
biomes_fp           = r"biomes_reprojected.shp" ## Olson et al 2001
shield_fp           = r"canadian_shield_reprojected.shp" ## Natural Resources Canada 2022

output_dir          = r''
output_gpkg         = os.path.join(output_dir, "PLD_PF.gpkg")
output_csv          = os.path.join(output_dir, "PLD_PF.csv")

# Target CRS for ALL spatial operations -- North Pole LAEA
TARGET_CRS = "EPSG:3575"
CHUNK_SIZE = 100_000   # lakes per chunk for spatial join

# ---------------------------------------------------------------------------
# Biome mapping
# ---------------------------------------------------------------------------
BIOME_MAPPING = {
    1:    "Tropical & Subtropical Moist Broadleaf Forests",
    2:    "Tropical & Subtropical Dry Broadleaf Forests",
    3:    "Tropical & Subtropical Coniferous Forests",
    4:    "Temperate Broadleaf & Mixed Forests",
    5:    "Temperate Conifer Forests",
    6:    "Boreal Forests/Taiga",
    7:    "Tropical & Subtropical Grasslands, Savannas & Shrublands",
    8:    "Temperate Grasslands, Savannas & Shrublands",
    9:    "Flooded Grasslands & Savannas",
    10:   "Montane Grasslands & Shrublands",
    11:   "Tundra",
    12:   "Mediterranean Forests, Woodlands & Scrub",
    13:   "Deserts & Xeric Shrublands",
    14:   "Mangroves",
    98:   "Water",
    99:   "Ice",
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def fix_geometries(gdf):
    """Fix invalid and null geometries."""
    n_null = gdf.geometry.isna().sum()
    if n_null > 0:
        print(f"  Dropping {n_null} null geometries...")
        gdf = gdf[gdf.geometry.notna()].copy()
    invalid = ~gdf.geometry.is_valid
    if invalid.any():
        print(f"  Fixing {invalid.sum()} invalid geometries...")
        gdf.loc[invalid, 'geometry'] = gdf.loc[invalid, 'geometry'].apply(
            lambda g: make_valid(g) if g is not None else g
        )
    return gdf


def to_target_crs(gdf):
    """Reproject a GeoDataFrame to TARGET_CRS if not already there."""
    if gdf.crs is None:
        raise ValueError("GeoDataFrame has no CRS set.")
    if gdf.crs.to_epsg() != int(TARGET_CRS.split(':')[1]):
        return gdf.to_crs(TARGET_CRS)
    return gdf


def assign_by_largest_intersection(lakes, overlay_gdf, field, label):
    """
    For each lake, assign the value of `field` from the polygon in overlay_gdf
    that has the largest intersection area with the lake.
    Returns lakes GeoDataFrame with the new field added.
    """
    print(f"  Assigning {label} ({field})...")
    lakes = lakes.copy()
    lakes['_orig_idx'] = lakes.index

    overlay_slim = overlay_gdf[[field, 'geometry']].copy()

    intersections = gpd.overlay(
        lakes[['_orig_idx', 'geometry']],
        overlay_slim,
        how='intersection'
    )
    intersections['_area'] = intersections.geometry.area

    if intersections.empty:
        lakes[field] = np.nan
    else:
        max_idx = intersections.groupby('_orig_idx')['_area'].idxmax()
        largest = intersections.loc[max_idx, ['_orig_idx', field]].set_index('_orig_idx')
        lakes[field] = lakes['_orig_idx'].map(largest[field])

    lakes.drop(columns=['_orig_idx'], inplace=True)
    return lakes


def assign_glaciated(lakes, glacial_gdf, label="glaciated status"):
    """
    Mark each lake as glaciated (1) or unglaciated (0) based on whether
    its geometry intersects the LGM glacial extent.
    """
    print(f"  Assigning {label}...")
    lakes = lakes.copy()
    lakes['_orig_idx'] = lakes.index

    intersections = gpd.overlay(
        lakes[['_orig_idx', 'geometry']],
        glacial_gdf[['geometry']],
        how='intersection'
    )
    intersections['_area'] = intersections.geometry.area
    glaciated_idxs = set(
        intersections.groupby('_orig_idx')['_area'].idxmax().index.tolist()
    ) if not intersections.empty else set()

    lakes['glaciated'] = lakes['_orig_idx'].isin(glaciated_idxs).astype(int)
    lakes.drop(columns=['_orig_idx'], inplace=True)
    return lakes

def assign_shield(lakes, shield_gdf, label="Canadian Shield"):
    """
    Mark each lake as within the Canadian Shield (True) or not (False)
    based on whether its geometry intersects the shield boundary.
    """
    print(f"  Assigning {label}...")
    lakes = lakes.copy()
    lakes['_orig_idx'] = lakes.index

    intersections = gpd.overlay(
        lakes[['_orig_idx', 'geometry']],
        shield_gdf[['geometry']],
        how='intersection'
    )
    shield_idxs = set(intersections['_orig_idx'].unique()) if not intersections.empty else set()

    lakes['shield'] = lakes['_orig_idx'].isin(shield_idxs)
    lakes.drop(columns=['_orig_idx'], inplace=True)
    return lakes


def assign_yedoma(lakes, yedoma_gdf, label="yedoma"):
    """
    Mark each lake as yedoma (1) or not (0) based on whether its geometry
    intersects the yedoma domain. 
    Output values: 'yedoma' / 'non-yedoma' 
    """
    print(f"  Assigning {label}...")
    lakes = lakes.copy()
    lakes['_orig_idx'] = lakes.index

    intersections = gpd.overlay(
        lakes[['_orig_idx', 'geometry']],
        yedoma_gdf[['geometry']],
        how='intersection'
    )
    yedoma_idxs = set(intersections['_orig_idx'].unique()) if not intersections.empty else set()

    lakes['yedoma'] = lakes['_orig_idx'].isin(yedoma_idxs).map({True: 'yedoma', False: 'non-yedoma'})
    lakes.drop(columns=['_orig_idx'], inplace=True)
    return lakes


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    os.makedirs(output_dir, exist_ok=True)

    # -----------------------------------------------------------------------
    # Step 1: Load lakes and reproject to TARGET_CRS
    # -----------------------------------------------------------------------
    print("=" * 60)
    print("Step 1: Loading lakes")
    print("=" * 60)

    print("Loading SWOT PLD lakes...")
    lakes_all = gpd.read_file(lakes_gdb_path, layer='SWOT_PLD_v201_02042025_attributes_updated')
    print(f"  Total lakes: {len(lakes_all):,}")

    # Quick lat filter before the expensive reproject
    lakes_north = lakes_all[lakes_all['lat'] > 50].copy()
    del lakes_all
    gc.collect()
    print(f"  Lakes north of 50 deg N: {len(lakes_north):,}")

    # Reproject to TARGET_CRS before any spatial operations.
    # EPSG:3575 is a planar projection centred on the North Pole, so there is
    # no antimeridian and polygon edges near 180 deg lon project correctly.
    print(f"  Reprojecting lakes to {TARGET_CRS}...")
    lakes_north = to_target_crs(lakes_north)
    lakes_north = fix_geometries(lakes_north)

    # -----------------------------------------------------------------------
    # Step 2: Filter to lakes that intersect permafrost zones
    # -----------------------------------------------------------------------
    print("=" * 60)
    print("Step 2: Filtering to permafrost zones")
    print("=" * 60)

    print("  Loading permafrost shapefile...")
    permafrost = gpd.read_file(permafrost_fp)
    permafrost = permafrost[permafrost['EXTENT'].isin(['C', 'D', 'I', 'S'])].copy()
    permafrost = to_target_crs(permafrost)
    permafrost = fix_geometries(permafrost)
    print(f"  Permafrost features: {len(permafrost)}")

    print("  Filtering lakes (chunked spatial join in EPSG:3575)...")
    results = []
    for i in range(0, len(lakes_north), CHUNK_SIZE):
        chunk = lakes_north.iloc[i:i + CHUNK_SIZE].copy()
        joined = gpd.sjoin(
            chunk,
            permafrost[['EXTENT', 'geometry']],
            how='inner',
            predicate='intersects'
        )
        joined = joined.drop(columns=['index_right'])
        results.append(joined)
        print(f"  Chunk {i // CHUNK_SIZE + 1}: {len(joined)} lakes intersect permafrost")
        gc.collect()

    lakes_pf = pd.concat(results, ignore_index=True)
    lakes_pf = lakes_pf.drop_duplicates(subset='lake_id')
    lakes_pf = gpd.GeoDataFrame(lakes_pf, geometry='geometry', crs=TARGET_CRS)
    del results, lakes_north
    gc.collect()
    print(f"  Lakes intersecting permafrost zones: {len(lakes_pf):,}")

    # -----------------------------------------------------------------------
    # Step 3: Exclude lakes that intersect the glacier shapefile
    # -----------------------------------------------------------------------
    print("=" * 60)
    print("Step 3: Excluding lakes that intersect glaciers")
    print("=" * 60)

    print("  Loading glacier shapefile...")
    glaciers = gpd.read_file(glaciers_fp)
    glaciers = to_target_crs(glaciers)
    glaciers = fix_geometries(glaciers)
    print(f"  Glacier features: {len(glaciers):,}")

    glacier_join = gpd.sjoin(
        lakes_pf[['lake_id', 'geometry']],
        glaciers[['geometry']],
        how='inner',
        predicate='intersects'
    )
    glacier_lake_ids = set(glacier_join['lake_id'].unique())
    print(f"  Lakes intersecting glaciers (excluded): {len(glacier_lake_ids):,}")

    lakes_pf = lakes_pf[~lakes_pf['lake_id'].isin(glacier_lake_ids)].copy()
    print(f"  Lakes remaining: {len(lakes_pf):,}")
    del glaciers, glacier_join, glacier_lake_ids
    gc.collect()

    # -----------------------------------------------------------------------
    # Step 4: Load attribute shapefiles
    # -----------------------------------------------------------------------
    print("=" * 60)
    print("Step 4: Loading attribute shapefiles")
    print("=" * 60)

    print("  Permafrost (full attributes)...")
    pf_attrs = gpd.read_file(permafrost_fp)
    pf_attrs = pf_attrs[pf_attrs['EXTENT'].isin(['C', 'D', 'I', 'S'])].copy()
    pf_attrs = to_target_crs(pf_attrs)
    pf_attrs = fix_geometries(pf_attrs)
    pf_attrs['CONTENT'] = pf_attrs['CONTENT'].str.strip().str.lower().map(
        {'l': 'Low', 'm': 'Medium', 'h': 'High'}
    )

    print("  Glacial history...")
    glacial = gpd.read_file(glacial_history_fp)
    glacial = to_target_crs(glacial)
    glacial = fix_geometries(glacial)

    print("  Yedoma...")
    yedoma = gpd.read_file(yedoma_fp)
    yedoma = to_target_crs(yedoma)
    yedoma = fix_geometries(yedoma)
    # No confidence column — polygon presence = yedoma (binary, like glaciated)

    print("  Thermokarst...")
    thermokarst = gpd.read_file(thermokarst_fp)
    thermokarst = to_target_crs(thermokarst)
    thermokarst = fix_geometries(thermokarst)
    tk_rename = {}
    for col in thermokarst.columns:
        if col.lower() in ('tkthlp', 'tkth_lp', 'tkthlp '):
            tk_rename[col] = 'TkThLP'
        elif col.lower() in ('tkwp', 'tk_wp', 'tkwp '):
            tk_rename[col] = 'TKWP'
    if tk_rename:
        thermokarst = thermokarst.rename(columns=tk_rename)

    print("  Biomes...")
    biomes = gpd.read_file(biomes_fp)
    biomes = to_target_crs(biomes)
    biomes = fix_geometries(biomes)
    biome_col = next((c for c in biomes.columns if c.upper() == 'BIOME'), 'BIOME')
    biomes[biome_col] = pd.to_numeric(biomes[biome_col], errors='coerce').astype('Int64')
    biomes['biome_name'] = biomes[biome_col].map(BIOME_MAPPING)

    print("  Canadian Shield...")
    shield = gpd.read_file(shield_fp)
    shield = to_target_crs(shield)
    shield = fix_geometries(shield)

    # -----------------------------------------------------------------------
    # Step 5: Assign attributes
    # -----------------------------------------------------------------------
    print("=" * 60)
    print("Step 5: Assigning attributes")
    print("=" * 60)

    lakes_out = lakes_pf.copy()

    lakes_out = assign_by_largest_intersection(lakes_out, pf_attrs, 'EXTENT',  'permafrost extent')
    lakes_out = assign_by_largest_intersection(lakes_out, pf_attrs, 'CONTENT', 'ground ice content')

    lakes_out = assign_glaciated(lakes_out, glacial)
    lakes_out['glaciated'] = lakes_out['glaciated'].map({1: 'glaciated', 0: 'unglaciated'})

    lakes_out = assign_yedoma(lakes_out, yedoma)

    if 'TkThLP' in thermokarst.columns:
        lakes_out = assign_by_largest_intersection(
            lakes_out, thermokarst, 'TkThLP', 'thermokarst lakes')
        lakes_out = lakes_out.rename(columns={'TkThLP': 'thermokarst_lakes'})
    if 'TKWP' in thermokarst.columns:
        lakes_out = assign_by_largest_intersection(
            lakes_out, thermokarst, 'TKWP', 'thermokarst wetlands')
        lakes_out = lakes_out.rename(columns={'TKWP': 'thermokarst_wetlands'})

    lakes_out = assign_by_largest_intersection(lakes_out, biomes, 'biome_name', 'biome')

    lakes_out = assign_shield(lakes_out, shield)

    # -----------------------------------------------------------------------
    # Step 6: Save outputs
    # -----------------------------------------------------------------------
    print("=" * 60)
    print("Step 6: Saving outputs")
    print("=" * 60)

    print(f"  Writing GeoPackage to:\n    {output_gpkg}")
    lakes_out.to_file(output_gpkg, layer='PLD_PF', driver='GPKG')

    print(f"  Writing CSV to:\n    {output_csv}")
    lakes_out.drop(columns='geometry').to_csv(output_csv, index=False)

    print("\nDone.")
    print(f"  Final lake count: {len(lakes_out):,}")
    print(f"  Columns: {lakes_out.columns.tolist()}")


if __name__ == "__main__":
    main()
