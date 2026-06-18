library(data.table)


# ────────────────────────────────────────────────────────────────
# ── load & filter data ──────────────────────────────────────────
# ────────────────────────────────────────────────────────────────

data <- freadPF_lake_density_coverage.csv')
globaldata <- fread('global_lake_density_fraction.csv')

# Remove grid cells with insufficient land coverage
# 50x50 km grid cells = 2500 km² total area; 10% threshold = 250 km²
AREA_THRESHOLD_KM2 <- 2500 * 0.10  # 250 km²

# ────────────────────────────────────────────────────────────────
# ── helper functions ────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
# Compute mean, SD, median, IQR — only for cells where the zone has area >= 10% of grid cell
summarise_grid <- function(x, area) {
  keep <- !is.na(area) & area >= AREA_THRESHOLD_KM2
  v    <- x[keep]
  list(
    avg = mean(v,   na.rm = TRUE),
    sd  = sd(v,     na.rm = TRUE),
    med = median(v, na.rm = TRUE),
    iqr = IQR(v,    na.rm = TRUE)
  )
}

make_row <- function(subcategory, frac_col, dens_col, area_col) {
  sf <- summarise_grid(frac_col, area_col)
  sd <- summarise_grid(dens_col, area_col)
  data.table(
    Category                                 = NA_character_,
    Subcategory                              = subcategory,
    `Average lake coverage (%)`                = round(sf$avg * 100, 1),
    `SD lake coverage (%)`                     = round(sf$sd  * 100, 1),
    `Median lake coverage (%)`                 = round(sf$med * 100, 1),
    `IQR lake coverage (%)`                    = round(sf$iqr * 100, 1),
    `Average lake density (# lakes/100km2)`  = round(sd$avg, 1),
    `SD lake density`                        = round(sd$sd,  1),
    `Median lake density (# lakes/100km2)`   = round(sd$med, 1),
    `IQR lake density (# lakes/100km2)`      = round(sd$iqr, 1)
  )
}

header_row <- function(category) {
  data.table(
    Category                                 = category,
    Subcategory                              = NA_character_,
    `Average lake coverage (%)`                = NA_real_,
    `SD lake coverage (%)`                     = NA_real_,
    `Median lake coverage (%)`                 = NA_real_,
    `IQR lake coverage (%)`                    = NA_real_,
    `Average lake density (# lakes/100km2)`  = NA_real_,
    `SD lake density`                        = NA_real_,
    `Median lake density (# lakes/100km2)`   = NA_real_,
    `IQR lake density (# lakes/100km2)`      = NA_real_
  )
}

# ──────────────────────────────────────────────────────────────
# ── 1. Permafrost conditions ──────────────────────────────────
# ──────────────────────────────────────────────────────────────

rows_pf <- rbindlist(list(
  header_row("Permafrost conditions"),
  make_row("Isolated",      data$lake_fraction_pf_I, data$lake_density_pf_I, data$pf_I_km2),
  make_row("Sporadic",      data$lake_fraction_pf_S, data$lake_density_pf_S, data$pf_S_km2),
  make_row("Discontinuous", data$lake_fraction_pf_D, data$lake_density_pf_D, data$pf_D_km2),
  make_row("Continuous",    data$lake_fraction_pf_C, data$lake_density_pf_C, data$pf_C_km2)
))

# ──────────────────────────────────────────────────────────────
# ── 2. Glacial history ────────────────────────────────────────
# ──────────────────────────────────────────────────────────────

rows_glacial <- rbindlist(list(
  header_row("Glacial history"),
  make_row("Not Glaciated", data$lake_fraction_pf_unglaciated, data$lake_density_pf_unglaciated, data$pf_unglaciated_km2),
  make_row("Glaciated",     data$lake_fraction_pf_glaciated,   data$lake_density_pf_glaciated,   data$pf_glaciated_km2)
))

# ──────────────────────────────────────────────────────────────
# ── 3. Thermokarst conditions ─────────────────────────────────
# ──────────────────────────────────────────────────────────────

rows_tk <- rbindlist(list(
  header_row("Thermokarst conditions"),
  header_row("Thermokarst Lakes"),
  make_row("None",           data$lake_fraction_pf_TkThLP_None,     data$lake_density_pf_TkThLP_None,     data$pf_TkThLP_None_km2),
  make_row("Low/moderate",   data$lake_fraction_pf_TkThLP_Moderate, data$lake_density_pf_TkThLP_Moderate, data$pf_TkThLP_Moderate_km2),
  make_row("High/very high", data$lake_fraction_pf_TkThLP_High,     data$lake_density_pf_TkThLP_High,     data$pf_TkThLP_High_km2),
  header_row("Thermokarst Wetlands"),
  make_row("None",           data$lake_fraction_pf_TKWP_None,       data$lake_density_pf_TKWP_None,       data$pf_TKWP_None_km2),
  make_row("Low/moderate",   data$lake_fraction_pf_TKWP_Moderate,   data$lake_density_pf_TKWP_Moderate,   data$pf_TKWP_Moderate_km2),
  make_row("High/very high", data$lake_fraction_pf_TKWP_High,       data$lake_density_pf_TKWP_High,       data$pf_TKWP_High_km2)
))

# ──────────────────────────────────────────────────────────────
# ── 4. Ground Ice Content ─────────────────────────────────────
# ──────────────────────────────────────────────────────────────

rows_ground_ice <- rbindlist(list(
  header_row("Ground Ice Content"),
  make_row("Low",    data$lake_fraction_pf_ice_Low,    data$lake_density_pf_ice_Low,    data$pf_ice_Low_km2),
  make_row("Medium", data$lake_fraction_pf_ice_Medium, data$lake_density_pf_ice_Medium, data$pf_ice_Medium_km2),
  make_row("High",   data$lake_fraction_pf_ice_High,   data$lake_density_pf_ice_High,   data$pf_ice_High_km2)
))

# ──────────────────────────────────────────────────────────────
# ── 5. Biomes ─────────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────

rows_biomes <- rbindlist(list(
  header_row("Biomes"),
  make_row("Boreal Forests/Taiga", data$lake_fraction_pf_biome6,  data$lake_density_pf_biome6,  data$pf_biome6_km2),
  make_row("Tundra",               data$lake_fraction_pf_biome11, data$lake_density_pf_biome11, data$pf_biome11_km2)
))

# ──────────────────────────────────────────────────────────────
# ── 6. Yedoma ─────────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────

rows_yedoma <- rbindlist(list(
  header_row("Yedoma"),
  make_row(" ", data$lake_fraction_pf_yedoma, data$lake_density_pf_yedoma, data$pf_yedoma_km2)
))

# ──────────────────────────────────────────────────────────────
# ── 7. Canadian Shield ────────────────────────────────────────
# ──────────────────────────────────────────────────────────────

row_shield <- make_row("Canadian Shield",
                       data$lake_fraction_pf_Shield,
                       data$lake_density_pf_Shield,
                       data$pf_Shield_km2)
row_shield[, Category    := "Canadian Shield"]
row_shield[, Subcategory := NA_character_]

# ──────────────────────────────────────────────────────────────
# ── 8. Totals rows ────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────

pf_stats_frac <- summarise_grid(data$lake_fraction_pf_area, data$pf_area_km2)
pf_stats_dens <- summarise_grid(data$lake_density_pf_area,  data$pf_area_km2)

row_pf_total <- data.table(
  Category                                 = "Entire northern permafrost domain",
  Subcategory                              = NA_character_,
  `Average lake coverage (%)`                = round(pf_stats_frac$avg * 100, 1),
  `SD lake coverage (%)`                     = round(pf_stats_frac$sd  * 100, 1),
  `Median lake coverage (%)`                 = round(pf_stats_frac$med * 100, 1),
  `IQR lake coverage (%)`                    = round(pf_stats_frac$iqr * 100, 1),
  `Average lake density (# lakes/100km2)`  = round(pf_stats_dens$avg, 1),
  `SD lake density`                        = round(pf_stats_dens$sd,  1),
  `Median lake density (# lakes/100km2)`   = round(pf_stats_dens$med, 1),
  `IQR lake density (# lakes/100km2)`      = round(pf_stats_dens$iqr, 1)
)

# Global row — calculated from globaldata using same area threshold
global_frac <- summarise_grid(globaldata$lake_fraction,          globaldata$non_glacier_land_area_km2)
global_dens <- summarise_grid(globaldata$lake_density_per100km2, globaldata$non_glacier_land_area_km2)

row_global <- data.table(
  Category                                 = "Globally",
  Subcategory                              = NA_character_,
  `Average lake coverage (%)`                = round(global_frac$avg * 100, 1),
  `SD lake coverage (%)`                     = round(global_frac$sd  * 100, 1),
  `Median lake coverage (%)`                 = round(global_frac$med * 100, 1),
  `IQR lake coverage (%)`                    = round(global_frac$iqr * 100, 1),
  `Average lake density (# lakes/100km2)`  = round(global_dens$avg, 1),
  `SD lake density`                        = round(global_dens$sd,  1),
  `Median lake density (# lakes/100km2)`   = round(global_dens$med, 1),
  `IQR lake density (# lakes/100km2)`      = round(global_dens$iqr, 1)
)

# ──────────────────────────────────────────────────────────────
# ── assemble & write ──────────────────────────────────────────
# ──────────────────────────────────────────────────────────────

table2 <- rbindlist(list(
  rows_pf,
  rows_glacial,
  rows_tk,
  rows_ground_ice,
  rows_biomes,
  rows_yedoma,
  row_shield,
  row_pf_total,
  row_global
), fill = TRUE)

fwrite(table2, 'lake_dist_table2_descending.csv')
