#!/usr/bin/env python3
"""
================================
Calculates total land area (km²) within the northern permafrost domain
(>50°N, EXTENT = C/D/I/S) for each landscape attribute category.

Uses the same direct union/intersection approach as pf_glaciated_fraction.py:
  - Clip permafrost to >50°N
  - Build unary_union for each category
  - Compute area(pf_union ∩ category_union) with glacier subtraction

Output
------
  northern_pf_land_area.csv
"""

import warnings
import numpy as np
import geopandas as gpd
import pandas as pd
from shapely.ops import unary_union
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.validation import make_valid

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_BASE = ()
_IN = ()

permafrost_fp  = _BASE + "permafrost_clean_reprojected.shp" ## Brown et al., 2002
glaciers_fp    = _BASE + "glaciers_reprojected.shp" ## GLIMS and NSIDC 2026
glacial_fp     = _BASE + "LGM_best_estimate.shp" ## Batchelor et al., 2025
thermokarst_fp = _BASE + "Thermokarst_reprojected.shp" ##Olefeldt et al., 2016
biomes_fp      = _BASE + "biomes_reprojected.shp" ## Olson et al., 2001
yedoma_fp      = _IN   + "yedoma_domain_reprojected.shp" ## Strauss et al., 2021
shield_fp      = _IN   + "/canadian_shield_reprojected.shp" ##Natural Resources Canada 2022

output_csv = ("lake_density_fraction/northern_pf_land_area.csv")

TARGET_CRS = "EPSG:3575"

# ---------------------------------------------------------------------------
# Geometry helpers — identical to pf_glaciated_fraction.py
# ---------------------------------------------------------------------------
def safe_union(gdf):
    geoms = [
        g.buffer(0) if (g is not None and not g.is_valid) else g
        for g in gdf.geometry if g is not None and not g.is_empty
    ]
    return unary_union(geoms) if geoms else None

def safe_intersect(a, b):
    if a is None or b is None or a.is_empty or b.is_empty:
        return None
    try:
        return a.intersection(b)
    except Exception:
        try:    return make_valid(a).intersection(make_valid(b))
        except: return None

def subtract(geom, mask):
    if geom is None or mask is None or mask.is_empty:
        return geom
    try:    return geom.difference(mask)
    except: return make_valid(geom).difference(make_valid(mask))

def area_km2(geom):
    if geom is None or geom.is_empty:
        return 0.0
    return geom.area / 1e6

def compute(pf_geom, cat_union):
    """area(pf ∩ cat) with glaciers already subtracted from pf_geom."""
    g = safe_intersect(pf_geom, cat_union)
    return round(area_km2(g))

def to_crs(gdf):
    if gdf.crs is None:
        return gdf.set_crs(TARGET_CRS, allow_override=True)
    if gdf.crs.to_epsg() != 3575:
        return gdf.to_crs(TARGET_CRS)
    return gdf

# ---------------------------------------------------------------------------
# 1. Load permafrost, clip to 50°N
# ---------------------------------------------------------------------------
print("Loading permafrost and clipping to north of 50°N...")
pf = to_crs(gpd.read_file(permafrost_fp))
pf['EXTENT']  = pf['EXTENT'].astype(str).str.strip()
pf['CONTENT'] = pf['CONTENT'].astype(str).str.strip()
pf = pf[pf['EXTENT'].isin(['C','D','I','S']) & pf.geometry.notna() & ~pf.geometry.is_empty].copy()

lons      = np.linspace(-180, 179.9, 720)
ring      = [(lon, 50.0) for lon in lons] + [(lon, 89.9) for lon in lons[::-1]]
north_gdf = gpd.GeoDataFrame(geometry=[ShapelyPolygon(ring)], crs="EPSG:4326").to_crs(TARGET_CRS)
pf        = gpd.clip(pf, north_gdf)
pf        = pf[pf.geometry.notna() & ~pf.geometry.is_empty].copy()
print(f"  {len(pf)} permafrost polygons after 50°N clip")

# ---------------------------------------------------------------------------
# 2. Build glacier union; subtract from pf to get pf_union (net of glaciers)
# ---------------------------------------------------------------------------
print("Building glacier union...")
glaciers      = to_crs(gpd.read_file(glaciers_fp))
glacier_union = safe_union(glaciers)
print(f"  Glacier area: {area_km2(glacier_union):,.0f} km²")

print("Building permafrost union (may take a few minutes)...")
pf_union_raw = safe_union(pf)
pf_union     = subtract(pf_union_raw, glacier_union)
print(f"  Total pf area: {area_km2(pf_union):,.0f} km²")

# Per-extent unions (also glacier-subtracted)
print("Building per-extent unions...")
ext_unions = {}
for ext in ['C','D','I','S']:
    sub = pf[pf['EXTENT']==ext]
    ext_unions[ext] = subtract(safe_union(sub), glacier_union)
    print(f"  Extent {ext}: {area_km2(ext_unions[ext]):,.0f} km²")

# Ground ice unions (from CONTENT column, no extra shapefile needed)
print("Building ground ice unions...")
ICE_MAP = {'l':'Low','m':'Medium','h':'High','low':'Low','medium':'Medium','high':'High'}
pf['_ice'] = pf['CONTENT'].str.lower().map(ICE_MAP)
ice_unions = {}
for level in ['Low','Medium','High']:
    sub = pf[pf['_ice']==level]
    ice_unions[level] = subtract(safe_union(sub), glacier_union)
    print(f"  Ice {level}: {area_km2(ice_unions[level]):,.0f} km²")

# ---------------------------------------------------------------------------
# 3. Load and union each category shapefile
# ---------------------------------------------------------------------------
print("Loading thermokarst...")
tk         = to_crs(gpd.read_file(thermokarst_fp))
tkthlp_col = next((c for c in tk.columns if c.lower() in ('tkthlp','tkth_lp')), None)
tkwp_col   = next((c for c in tk.columns if c.lower() in ('tkwp','tk_wp')),     None)
HIGH_VALS  = ['very high','very_high','high']
MOD_VALS   = ['low','moderate','low/moderate','moderate/low']

def tk_union(col, vals):
    mask = tk[col].astype(str).str.strip().str.lower().isin(vals)
    return subtract(safe_union(tk[mask]), glacier_union)

def tk_none_union(col, h_vals, m_vals):
    mask = ~tk[col].astype(str).str.strip().str.lower().isin(h_vals + m_vals)
    return subtract(safe_union(tk[mask]), glacier_union)

print("  Building thermokarst unions...")
tklp_high_u = tk_union(tkthlp_col, HIGH_VALS)
tklp_mod_u  = tk_union(tkthlp_col, MOD_VALS)
tklp_none_u = tk_none_union(tkthlp_col, HIGH_VALS, MOD_VALS)
tkwp_high_u = tk_union(tkwp_col, HIGH_VALS)
tkwp_mod_u  = tk_union(tkwp_col, MOD_VALS)
tkwp_none_u = tk_none_union(tkwp_col, HIGH_VALS, MOD_VALS)
print("  Thermokarst unions done.")

print("Loading yedoma...")
yedoma_u = subtract(safe_union(to_crs(gpd.read_file(yedoma_fp))), glacier_union)
print(f"  Yedoma: {area_km2(yedoma_u):,.0f} km²")

print("Loading biomes...")
biomes    = to_crs(gpd.read_file(biomes_fp))
biome_col = next((c for c in biomes.columns if c.upper() == 'BIOME'), None)
biomes[biome_col] = pd.to_numeric(biomes[biome_col], errors='coerce')
biome6_u  = subtract(safe_union(biomes[biomes[biome_col]==6]),  glacier_union)
biome11_u = subtract(safe_union(biomes[biomes[biome_col]==11]), glacier_union)
print(f"  Biome6: {area_km2(biome6_u):,.0f}  Biome11: {area_km2(biome11_u):,.0f} km²")

print("Loading glacial history...")
lgm       = to_crs(gpd.read_file(glacial_fp))
lgm       = lgm[lgm.geometry.notna() & ~lgm.geometry.is_empty].copy()
lgm_union = subtract(safe_union(lgm), glacier_union)
print(f"  LGM: {area_km2(lgm_union):,.0f} km²")

print("Loading Canadian Shield...")
shield_gdf = to_crs(gpd.read_file(shield_fp))
shield_u   = subtract(safe_union(shield_gdf), glacier_union)
print(f"  Shield: {area_km2(shield_u):,.0f} km²")

# ---------------------------------------------------------------------------
# 4. Compute all areas
# ---------------------------------------------------------------------------
print("\nComputing areas...")
total_pf    = round(area_km2(pf_union))
C_area      = round(area_km2(ext_unions['C']))
D_area      = round(area_km2(ext_unions['D']))
S_area      = round(area_km2(ext_unions['S']))
I_area      = round(area_km2(ext_unions['I']))
ice_high    = round(area_km2(ice_unions['High']))
ice_med     = round(area_km2(ice_unions['Medium']))
ice_low     = round(area_km2(ice_unions['Low']))

tklp_high   = compute(pf_union, tklp_high_u)
tklp_mod    = compute(pf_union, tklp_mod_u)
tklp_none   = compute(pf_union, tklp_none_u)
tkwp_high   = compute(pf_union, tkwp_high_u)
tkwp_mod    = compute(pf_union, tkwp_mod_u)
tkwp_none   = compute(pf_union, tkwp_none_u)
yedoma      = compute(pf_union, yedoma_u)
biome6      = compute(pf_union, biome6_u)
biome11     = compute(pf_union, biome11_u)
glaciated   = compute(pf_union, lgm_union)
unglaciated = round(total_pf - glaciated)
shield      = compute(pf_union, shield_u)

# ---------------------------------------------------------------------------
# 5. Build and output table
# ---------------------------------------------------------------------------
rows = [
    ("Permafrost conditions",             ""),
    ("Continuous",                        C_area),
    ("Discontinuous",                     D_area),
    ("Sporadic",                          S_area),
    ("Isolated",                          I_area),
    ("Thermokarst conditions",            ""),
    ("Thermokarst Lakes",                 ""),
    ("High/very high",                    tklp_high),
    ("Low/moderate",                      tklp_mod),
    ("None",                              tklp_none),
    ("Thermokarst Wetlands",              ""),
    ("High/very high",                    tkwp_high),
    ("Low/moderate",                      tkwp_mod),
    ("None",                              tkwp_none),
    ("Ground Ice Content",                ""),
    ("High",                              ice_high),
    ("Medium",                            ice_med),
    ("Low",                               ice_low),
    ("Yedoma",                            yedoma),
    ("Biomes",                            ""),
    ("Boreal Forests/Taiga",              biome6),
    ("Tundra",                            biome11),
    ("Glacial history",                   ""),
    ("Glaciated",                         glaciated),
    ("Not Glaciated",                     unglaciated),
    ("Canadian Shield",                   shield),
    ("Entire northern permafrost domain", total_pf),
]

print()
print(f"{'Category':<40} {'Land area (km²)':>16}")
print("-" * 58)
for label, val in rows:
    val_str = f"{val:>16,.0f}" if val != "" else ""
    print(f"{label:<40} {val_str}")
print()

df = pd.DataFrame(rows, columns=["Category", "Land area (km2)"])
df.to_csv(output_csv, index=False)
print(f"CSV written to:\n  {output_csv}")
