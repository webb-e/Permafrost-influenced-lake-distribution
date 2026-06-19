`global_PF_extent.py` Determine whether lakes intersect with glaciers and are within the permafrost region. Requires `PLD_PF.gpkg`, which is an output of `create_PF_PLD.py`; landscape characteristic shapefiles; and global lakes inventory; see main page for references.

`global_lake_density_fraction.py` Calculates global lake density and fraction (coverage). Relies on `PLD_PF.gpkg`, which is an output of `create_PF_PLD.py`;  `Grid_global_50km_non_glacier_area.gpkg`, which is an output of `global_area_grid.py`; and global lakes inventory, see main page for references

`global_grid_area.py` Calculates the area of land that is not currently occupied by glaciers using the global grid. Relies on `Grid_global_50km.shp`, which an output of `make_global_grid`; and GLIMS and NSIDC 2026; see main page for references

`create_PF_PLD.py` Determine which lakes are within the northern permafrost region and intersect lakes with geospatial datasets to assign landscape characteristic attributes. Relies on landscape characteristic shapefiles and global lakes inventory; see main page for references.

`glacial_pf_intersection_sheild.py` Stratify permafrost extent by glacial history and the Canadian Shield; determine lake density and coverage. Relies on `PLD_PF.gpkg`, which is an output of `create_PF_PLD.py`; `Northern_grid_50km.shp`, which is an output of `make_northern_grid`; and landscape characterstic shapefiles; see main paing for references.

`future_permafrost_analysis.py` Intersect PLD (lake inventory) with future warming scenarios from Cadburn et al., 2017. Relies on data from Chadburn et al., 2017 and Brown et al., 2002; see main page for references.

`make_global_grid` Run in Google Earth Engine code editor; creates global equal area grid

`make_northern_grid` Run in Google Earth Engine code editor; creates equal area grid for land north of 50 DegN

Code used to generate lake density and coverage across the northern permafrost region is archived through the Arctic Data Center (LINK), along with the GeoPackage of lake density and coverage.
