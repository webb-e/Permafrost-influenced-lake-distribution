#!/usr/bin/env python3
"""
Computes lake density and lake fraction for each permafrost extent class (C, D, I, S),
subdivided by glacial history (previously glaciated vs. unglaciated) and Canadian Shield,
across a 50 km equal-area grid.

For each permafrost-glacial category, the script calculates per grid cell:
  - Category area (m²): area of the permafrost-glacial category within the grid cell,
    with current glacier extent subtracted
  - Lake area: total area of lake polygons intersecting the category within the grid cell
  - Lake fraction: lake area / category area
  - Lake count: number of lakes with centroids (lon/lat attributes) falling within the category
  - Lake density: lake count / category area

Lake fraction is computed from polygon intersections; lake density is computed from lake centroids
using the lon/lat attribute fields on the lakes file, reprojected to EPSG:3575.

Categories:
  glac_{C/D/I/S}    glaciated × permafrost extent
  ungl_{C/D/I/S}    unglaciated × permafrost extent
  glac_sh_{C/D/I/S}  glaciated × Canadian Shield × permafrost extent
  glac_nsh_{C/D/I/S} glaciated × NOT Canadian Shield × permafrost extent
                     (shield is entirely glaciated, so no unglaciated shield categories)

Output: a GeoPackage (PF_x_Glaciation.gpkg) containing the grid with per-category
area, lake area, lake count, fraction, and density fields.
"""

import os
import geopandas as gpd
import fiona
from shapely.ops import unary_union
from shapely.geometry import GeometryCollection
from shapely.validation import make_valid
import numpy as np
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

# -------------------------
# File paths
# -------------------------
_BASE = ()
glacial_history_fp = os.path.join(_BASE, "LGM_best_estimate.shp") ## Batchelor et al 2025
permafrost_fp      = os.path.join(_BASE, "permafrost_clean_reprojected.shp") ## Brown et al 2002
grid_fp            = os.path.join(_BASE, "Northern_grid_50km.shp")
glaciers_fp        = os.path.join(_BASE, "glaciers_reprojected.shp") ## GLIMS and NSIDC 2026
shield_fp          = (r"canadian_shield_reprojected.shp") ## Natural Resources Canada 2022

lakes_fp = (PLD_PF.gpkg)

out_dir  = ()
out_path = os.path.join(out_dir, "PF_x_Glaciation.gpkg")

# Categories: glaciated and unglaciated × each permafrost extent
EXTENTS    = {"C", "D", "I", "S"}
CATEGORIES = []
for e in sorted(EXTENTS):
    CATEGORIES.append(("glac", e))
    CATEGORIES.append(("ungl", e))
    CATEGORIES.append(("glac_sh", e))    # glaciated × Canadian Shield
    CATEGORIES.append(("glac_nsh", e))  # glaciated × NOT Canadian Shield

target_crs = "EPSG:3575"

# -------------------------
# Geometry helpers (same pattern as grid_area_attributes.py)
# -------------------------

def ensure_crs(gdf, target=target_crs):
    from pyproj import CRS
    from pyproj.exceptions import CRSError
    proj4 = "+proj=laea +lat_0=90 +lon_0=10 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    try:    tcrs = CRS.from_user_input(target)
    except: tcrs = CRS.from_user_input(proj4)
    if gdf.crs is None:
        return gdf.set_crs(tcrs, allow_override=True)
    try:
        if CRS.from_user_input(gdf.crs) == tcrs:
            return gdf
    except Exception:
        pass
    try:    return gdf.to_crs(tcrs)
    except CRSError:
        try:    return gdf.to_crs(proj4)
        except: return gdf.set_crs(proj4, allow_override=True)
    except Exception:
        try:    return gdf.to_crs(proj4)
        except: return gdf.set_crs(proj4, allow_override=True)


def safe_union(gdf):
    if gdf is None or len(gdf) == 0:
        return GeometryCollection()
    geoms = [
        g.buffer(0) if (g is not None and not g.is_valid) else g
        for g in gdf.geometry
        if g is not None and not g.is_empty
    ]
    return unary_union(geoms) if geoms else GeometryCollection()


def safe_intersect(a, b):
    if a is None or a.is_empty or b is None or b.is_empty:
        return GeometryCollection()
    try:
        return a.intersection(b)
    except Exception:
        try:    return make_valid(a).intersection(make_valid(b))
        except: return GeometryCollection()


def subtract(geom, mask):
    """Subtract mask from geom; returns geom unchanged if mask is empty."""
    if mask is None or mask.is_empty:
        return geom
    try:    return geom.difference(mask)
    except: return make_valid(geom).difference(make_valid(mask))


# -------------------------
# 1. Load grid
# -------------------------
print("Loading grid...")
try:
    grid = gpd.read_file(grid_fp, engine='fiona')
except Exception as e:
    print(f"  fiona load failed ({e}); falling back.")
    with fiona.open(grid_fp) as src:
        grid = gpd.GeoDataFrame.from_features(list(src), crs=None)

if grid.crs is None:
    grid = grid.set_crs(target_crs, allow_override=True)
else:
    grid = ensure_crs(grid, target_crs)

if 'grid_id' not in grid.columns:
    grid = grid.reset_index(drop=True)
    grid['grid_id'] = grid.index.astype(int)

print(f"  {len(grid)} grid cells loaded.")

# -------------------------
# 2. Load permafrost, glacial history, and current glaciers
# -------------------------
print("Loading permafrost, glacial history, glacier, and shield layers...")
permafrost = ensure_crs(gpd.read_file(permafrost_fp), target_crs)
glac_hist  = ensure_crs(gpd.read_file(glacial_history_fp), target_crs)
glaciers   = ensure_crs(gpd.read_file(glaciers_fp), target_crs)
shield_gdf = ensure_crs(gpd.read_file(shield_fp), target_crs)

# Build current glacier union for subtraction
print("Building current glacier union...")
glacier_union = safe_union(glaciers)
print(f"  Glacier area: {glacier_union.area / 1e6:,.0f} km²")

# -------------------------
# 3. Filter permafrost to extents C, D, I, S
# -------------------------
permafrost['EXTENT'] = permafrost['EXTENT'].astype(str).str.strip()
pf_sel = permafrost[permafrost['EXTENT'].isin(EXTENTS)].copy()
if pf_sel.empty:
    raise ValueError("No permafrost polygons found with EXTENT in C, D, I, S.")
pf_sel['geometry'] = pf_sel.geometry.buffer(0)

# -------------------------
# 4. Build LGM glaciated union and Canadian Shield union
# -------------------------
print("Building LGM glaciated union...")
glac_hist = glac_hist[~glac_hist.geometry.is_empty & glac_hist.geometry.notna()].copy()
lgm_union = safe_union(glac_hist)
lgm_union = subtract(lgm_union, glacier_union)

print("Building Canadian Shield union...")
shield_union = safe_union(shield_gdf)
# Pre-compute lgm ∩ shield and lgm ∩ (not shield) to reuse across all extents
lgm_shield_union  = safe_intersect(lgm_union, shield_union)
lgm_nshield_union = subtract(lgm_union, shield_union)
print(f"  LGM ∩ Shield area:     {lgm_shield_union.area / 1e6:,.0f} km²")
print(f"  LGM ∩ not-Shield area: {lgm_nshield_union.area / 1e6:,.0f} km²")

# -------------------------
# 5. Build category geometries: permafrost × (glaciated / unglaciated),
#    with current glaciers subtracted from all areas
# -------------------------
print("Constructing category geometries...")
category_gdfs = {}
for extent in sorted(EXTENTS):
    sub = pf_sel[pf_sel['EXTENT'] == extent].copy()
    if sub.empty:
        for pref in ('glac', 'ungl'):
            category_gdfs[f'{pref}_{extent}'] = gpd.GeoDataFrame(
                geometry=[], crs=target_crs)
        continue

    inter_geoms = []
    diff_geoms  = []
    for geom in sub.geometry:
        if geom is None or geom.is_empty:
            inter_geoms.append(None)
            diff_geoms.append(None)
            continue
        # intersect with LGM extent, subtract current glaciers
        glac_part = subtract(safe_intersect(geom, lgm_union), glacier_union)
        ungl_part = subtract(geom.difference(lgm_union) if not lgm_union.is_empty else geom,
                             glacier_union)
        inter_geoms.append(glac_part if (glac_part is not None and not glac_part.is_empty) else None)
        diff_geoms.append(ungl_part  if (ungl_part  is not None and not ungl_part.is_empty)  else None)

    sub_glac = sub.copy(); sub_glac['geometry'] = inter_geoms
    sub_ungl = sub.copy(); sub_ungl['geometry'] = diff_geoms
    sub_glac = sub_glac[sub_glac.geometry.notna() & ~sub_glac.geometry.is_empty].copy()
    sub_ungl = sub_ungl[sub_ungl.geometry.notna() & ~sub_ungl.geometry.is_empty].copy()

    category_gdfs[f'glac_{extent}'] = sub_glac.set_crs(target_crs, allow_override=True) if not sub_glac.empty else gpd.GeoDataFrame(geometry=[], crs=target_crs)
    category_gdfs[f'ungl_{extent}'] = sub_ungl.set_crs(target_crs, allow_override=True) if not sub_ungl.empty else gpd.GeoDataFrame(geometry=[], crs=target_crs)

    # glaciated × Canadian Shield and glaciated × NOT Canadian Shield
    sh_geoms  = []
    nsh_geoms = []
    for geom in sub.geometry:
        if geom is None or geom.is_empty:
            sh_geoms.append(None);  nsh_geoms.append(None)
            continue
        sh_part  = subtract(safe_intersect(geom, lgm_shield_union),  glacier_union)
        nsh_part = subtract(safe_intersect(geom, lgm_nshield_union), glacier_union)
        sh_geoms.append(sh_part   if (sh_part  is not None and not sh_part.is_empty)  else None)
        nsh_geoms.append(nsh_part if (nsh_part is not None and not nsh_part.is_empty) else None)

    sub_sh  = sub.copy(); sub_sh['geometry']  = sh_geoms
    sub_nsh = sub.copy(); sub_nsh['geometry'] = nsh_geoms
    sub_sh  = sub_sh[sub_sh.geometry.notna()   & ~sub_sh.geometry.is_empty].copy()
    sub_nsh = sub_nsh[sub_nsh.geometry.notna() & ~sub_nsh.geometry.is_empty].copy()
    category_gdfs[f'glac_sh_{extent}']  = sub_sh.set_crs(target_crs,  allow_override=True) if not sub_sh.empty  else gpd.GeoDataFrame(geometry=[], crs=target_crs)
    category_gdfs[f'glac_nsh_{extent}'] = sub_nsh.set_crs(target_crs, allow_override=True) if not sub_nsh.empty else gpd.GeoDataFrame(geometry=[], crs=target_crs)

# -------------------------
# 6. Initialise output columns on grid
# -------------------------
print("Preparing grid output columns...")
grid = ensure_crs(grid.copy(), target_crs)
for pref, ext in CATEGORIES:
    grid[f"area_{pref}_{ext}"]      = 0.0
    grid[f"lakearea_{pref}_{ext}"]  = 0.0
    grid[f"lakecount_{pref}_{ext}"] = 0
    grid[f"dens_{pref}_{ext}"]      = np.nan
    grid[f"frac_{pref}_{ext}"]      = np.nan

# -------------------------
# 7. Load lakes; build centroid GeoDataFrame from lon/lat attributes
# -------------------------
print("Loading lakes...")
lakes = ensure_crs(gpd.read_file(lakes_fp), target_crs)
lakes = lakes[~lakes.geometry.is_empty & lakes.geometry.notna()].copy()
lakes = lakes.reset_index(drop=True).rename_axis('lake_idx').reset_index()
lakes['lake_geom_area'] = lakes.geometry.area

print("Building centroid GeoDataFrame from lon/lat attributes...")
centroids_gdf = gpd.GeoDataFrame(
    lakes[['lake_idx']].copy(),
    geometry=gpd.points_from_xy(lakes['lon'].astype(float), lakes['lat'].astype(float)),
    crs="EPSG:4326"
).to_crs(target_crs)

try:
    centroids_sindex = centroids_gdf.sindex
except Exception:
    centroids_sindex = None

# -------------------------
# 8. Per-category: overlay grid × category, compute area, fraction, density
# -------------------------
print("Processing categories...")
for pref, ext in CATEGORIES:
    key     = f"{pref}_{ext}"
    cat_gdf = category_gdfs.get(key)
    if cat_gdf is None or cat_gdf.empty:
        print(f"  {key}: no geometry, skipping.")
        continue

    print(f"  Processing {key} ...")
    area_col       = f"area_{pref}_{ext}"
    lake_area_col  = f"lakearea_{pref}_{ext}"
    lake_count_col = f"lakecount_{pref}_{ext}"
    frac_col       = f"frac_{pref}_{ext}"
    density_col    = f"dens_{pref}_{ext}"

    # --- Grid × category overlay ---
    try:
        grid_candidates = gpd.sjoin(
            grid[['grid_id', 'geometry']], cat_gdf[['geometry']],
            how='inner', predicate='intersects')
        grid_sub = grid[grid['grid_id'].isin(grid_candidates['grid_id'].unique())].copy()
        if grid_sub.empty:
            print(f"  {key}: no overlapping grid cells.")
            continue
        inter = gpd.overlay(grid_sub[['grid_id', 'geometry']], cat_gdf[['geometry']], how='intersection')
    except Exception:
        inter = gpd.overlay(grid[['grid_id', 'geometry']], cat_gdf[['geometry']], how='intersection')

    if inter.empty:
        print(f"  {key}: overlay empty.")
        continue

    inter = inter.set_crs(target_crs, allow_override=True)
    inter['cat_area'] = inter.geometry.area
    inter = inter.reset_index(drop=True).rename_axis('inter_idx').reset_index()

    # Aggregate category area per grid cell
    area_by_grid = inter.groupby('grid_id')['cat_area'].sum().to_dict()
    for gid, val in area_by_grid.items():
        grid.loc[grid['grid_id'] == gid, area_col] = float(val)

    inter_grid_map  = inter.set_index('inter_idx')['grid_id'].to_dict()
    inter_geom_dict = inter.set_index('inter_idx')['geometry'].to_dict()

    # --- Lake fraction (polygon intersection) ---
    try:
        possible_idx = []
        for bounds in inter.geometry.bounds.values:
            possible_idx.extend(list(lakes.sindex.intersection(tuple(bounds))))
        lakes_sub = lakes.loc[sorted(set(possible_idx))].copy() if possible_idx else lakes.iloc[0:0].copy()
    except Exception:
        lakes_sub = lakes.copy()

    lakearea_by_inter = {}
    if not lakes_sub.empty:
        try:
            joined = gpd.sjoin(lakes_sub.set_geometry('geometry'),
                               inter[['inter_idx', 'geometry']], how='inner', predicate='intersects')
        except Exception:
            joined = gpd.sjoin(lakes_sub.set_geometry('geometry'),
                               inter[['inter_idx', 'geometry']], how='inner', op='intersects')

        for _, row in joined.iterrows():
            inter_geom = inter_geom_dict.get(row['inter_idx'])
            if inter_geom is None or inter_geom.is_empty:
                continue
            try:
                area = row['geometry'].intersection(inter_geom).area
                if area > 0:
                    lakearea_by_inter[row['inter_idx']] = lakearea_by_inter.get(row['inter_idx'], 0.0) + area
            except Exception:
                continue

    lakearea_by_grid = {}
    for inter_idx, area_val in lakearea_by_inter.items():
        gid = inter_grid_map.get(inter_idx)
        if gid is not None:
            lakearea_by_grid[gid] = lakearea_by_grid.get(gid, 0.0) + float(area_val)

    for gid, val in lakearea_by_grid.items():
        grid.loc[grid['grid_id'] == gid, lake_area_col] = float(val)

    # Compute fraction
    for gid in set(area_by_grid) | set(lakearea_by_grid):
        area_val      = grid.loc[grid['grid_id'] == gid, area_col].values[0]
        lake_area_val = grid.loc[grid['grid_id'] == gid, lake_area_col].values[0]
        if area_val and not np.isnan(area_val) and area_val > 0:
            grid.loc[grid['grid_id'] == gid, frac_col] = float(lake_area_val) / float(area_val)

    # --- Lake density (centroid-based) ---
    lakecount_by_inter = {}
    try:
        possible_c_idx = []
        if centroids_sindex is not None:
            for b in inter.geometry.bounds.values:
                possible_c_idx.extend(list(centroids_sindex.intersection(tuple(b))))
            centroids_sub = centroids_gdf.iloc[sorted(set(possible_c_idx))].copy() if possible_c_idx else centroids_gdf.iloc[0:0].copy()
        else:
            centroids_sub = centroids_gdf.copy()

        if not centroids_sub.empty:
            try:
                c_joined = gpd.sjoin(centroids_sub.set_geometry('geometry'),
                                     inter[['inter_idx', 'geometry']], how='inner', predicate='within')
            except Exception:
                c_joined = gpd.sjoin(centroids_sub.set_geometry('geometry'),
                                     inter[['inter_idx', 'geometry']], how='inner', op='within')
            if not c_joined.empty:
                lakecount_by_inter = c_joined.groupby('inter_idx')['lake_idx'].nunique().to_dict()
    except Exception:
        pass

    lakecount_by_grid = {}
    for inter_idx, count_val in lakecount_by_inter.items():
        gid = inter_grid_map.get(inter_idx)
        if gid is not None:
            lakecount_by_grid[gid] = lakecount_by_grid.get(gid, 0) + int(count_val)

    for gid, val in lakecount_by_grid.items():
        grid.loc[grid['grid_id'] == gid, lake_count_col] = int(val)

    # Compute density
    for gid in set(area_by_grid) | set(lakecount_by_grid):
        area_val  = grid.loc[grid['grid_id'] == gid, area_col].values[0]
        count_val = grid.loc[grid['grid_id'] == gid, lake_count_col].values[0]
        if area_val and not np.isnan(area_val) and area_val > 0:
            grid.loc[grid['grid_id'] == gid, density_col] = float(count_val) / float(area_val)

    print(f"    {key}: {len(area_by_grid)} cells with area, "
          f"{len(lakearea_by_grid)} with lake area, "
          f"{len(lakecount_by_grid)} with lake centroids.")

# -------------------------
# 9. Write output
# -------------------------
print("Writing output...")
os.makedirs(out_dir, exist_ok=True)
if os.path.exists(out_path):
    try: os.remove(out_path)
    except Exception: pass

if 'grid_area' not in grid.columns:
    grid['grid_area'] = grid.geometry.area

grid.to_file(out_path, driver="GPKG", layer="grid_summary")
print(f"Output written to: {out_path}")

csv_path = out_path.replace(".gpkg", ".csv")
grid.drop(columns="geometry").to_csv(csv_path, index=False)
print(f"CSV written to: {csv_path}")

added_cols = [c for c in grid.columns if c.startswith(('area_', 'lakearea_', 'lakecount_', 'dens_', 'frac_'))]
print("Added columns:")
for c in sorted(added_cols):
    print(f"  - {c}")

print("Processing complete.")
