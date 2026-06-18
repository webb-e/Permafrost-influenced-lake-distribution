library(data.table)

# ────────────────────────────────────────────────────────────────
# ── load data ───────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────

pld <- fread('PF_PLD.csv')
land <- fread('northern_pf_land_area.csv')
setnames(land, c("Category", "land_area_km2"))

# ────────────────────────────────────────────────────────────────
# ── helper functions ────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────
summarise_lakes <- function(poly_area) {
  list(total_area = sum(poly_area, na.rm = TRUE),
       n_lakes    = length(poly_area),
       med_area   = median(poly_area, na.rm = TRUE),
       iqr_area   = IQR(poly_area,    na.rm = TRUE))
}

make_row <- function(category, subcategory, poly_area,
                     total_pf_area, total_pf_n,
                     global_area, global_n) {
  s <- summarise_lakes(poly_area)
  data.table(
    Category                                   = NA_character_,
    Subcategory                                = subcategory,
    `Total Lake area (km2)`                    = round(s$total_area, 0),
    `% of lake area within northern PF region` = round(s$total_area / total_pf_area * 100, 1),
    `% of global lake area`                    = round(s$total_area / global_area   * 100, 1),
    `Median lake area (km2)`                   = round(s$med_area, 3),
    `Lake area IQR (km2)`                      = round(s$iqr_area, 3),
    `Total number of lakes`                    = s$n_lakes,
    `% of lakes within northern PF region`     = round(s$n_lakes / total_pf_n * 100, 1),
    `% of global lakes`                        = round(s$n_lakes / global_n   * 100, 1))
}

# Blank header row with just the category name in col 0
header_row <- function(category) {
  data.table(
    Category                                   = category,
    Subcategory                                = NA_character_,
    `Total Lake area (km2)`                    = NA_real_,
    `% of lake area within northern PF region` = NA_real_,
    `% of global lake area`                    = NA_real_,
    `Median lake area (km2)`                   = NA_real_,
    `Lake area IQR (km2)`                      = NA_real_,
    `Total number of lakes`                    = NA_integer_,
    `% of lakes within northern PF region`     = NA_real_,
    `% of global lakes`                        = NA_real_
  )
}

# ──────────────────────────────────────────────────────────────
# ── lake area/lake number denominators ────────────────────────
# ──────────────────────────────────────────────────────────────

total_pf_area <- sum(pld$poly_area, na.rm = TRUE)
total_pf_n    <- nrow(pld)

global_area <- 2623118   # km²  (from global PLD)
global_n    <- 5897941   # lakes (from global PLD)

# ──────────────────────────────────────────────────────────────
# ── 1. Permafrost conditions ──────────────────────────────────
# ──────────────────────────────────────────────────────────────
# EXTENT: "C" = Continuous, "D" = Discontinuous, "S" = Sporadic, "I" = Isolated

pf_zones <- list(Isolated = "I", Sporadic = "S", Discontinuous = "D", Continuous = "C")

rows_pf <- rbindlist(c(
  list(header_row("Permafrost conditions")),
  lapply(names(pf_zones), function(label) {
    make_row("Permafrost conditions", label,
             pld[EXTENT == pf_zones[[label]], poly_area],
             total_pf_area, total_pf_n, global_area, global_n)
  })))

# ──────────────────────────────────────────────────────────────
# ── 2. Glacial history ────────────────────────────────────────
# ──────────────────────────────────────────────────────────────

pld[, glaciated_class := fifelse(
  tolower(trimws(glaciated)) == "glaciated", "Glaciated", "Not Glaciated")]

rows_glacial <- rbindlist(c(
  list(header_row("Glacial history")),
  lapply(c("Not Glaciated", "Glaciated"), function(g) {
    make_row("Glacial history", g,
             pld[glaciated_class == g, poly_area],
             total_pf_area, total_pf_n, global_area, global_n)
  })))

# ──────────────────────────────────────────────────────────────
# ── 3. Thermokarst conditions ─────────────────────────────────
# ──────────────────────────────────────────────────────────────
tk_levels <- c("None", "Low/moderate", "High/very high")

classify_tk <- function(col_vals) {
  out <- as.character(col_vals)
  out[out %in% c("High", "Very High")] <- "High/very high"
  out[out %in% c("Low", "Moderate")]   <- "Low/moderate"
  out[is.na(out) | out == ""]          <- "None"
  out
}

pld[, tk_lakes_class    := classify_tk(thermokarst_lakes)]
pld[, tk_wetlands_class := classify_tk(thermokarst_wetlands)]

rows_tk <- rbindlist(c(
  list(header_row("Thermokarst conditions")),
  list(header_row("Thermokarst Lakes")),
  lapply(tk_levels, function(lv) {
    make_row("Thermokarst conditions", lv,
             pld[tk_lakes_class == lv, poly_area],
             total_pf_area, total_pf_n, global_area, global_n)
  }),
  list(header_row("Thermokarst Wetlands")),
  lapply(tk_levels, function(lv) {
    make_row("Thermokarst conditions", lv,
             pld[tk_wetlands_class == lv, poly_area],
             total_pf_area, total_pf_n, global_area, global_n)
  })))

# ──────────────────────────────────────────────────────────────
# ── 4. Ground Ice Content ─────────────────────────────────────
# ──────────────────────────────────────────────────────────────
# CONTENT: "High", "Medium", "Low"

rows_ground_ice <- rbindlist(c(
  list(header_row("Ground Ice Content")),
  lapply(c("Low", "Medium", "High"), function(g) {
    make_row("Ground Ice Content", g,
             pld[CONTENT == g, poly_area],
             total_pf_area, total_pf_n, global_area, global_n)
  })))

# ──────────────────────────────────────────────────────────────
# ── 5. Biomes ─────────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────

rows_biomes <- rbindlist(list(
  header_row("Biomes"),
  make_row("Biomes", "Boreal Forests/Taiga",
           pld[biome_name == "Boreal Forests/Taiga", poly_area],
           total_pf_area, total_pf_n, global_area, global_n),
  make_row("Biomes", "Tundra",
           pld[biome_name == "Tundra", poly_area],
           total_pf_area, total_pf_n, global_area, global_n)))

# ──────────────────────────────────────────────────────────────
# ── 6. Yedoma ─────────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────

pld[, yedoma_class := fifelse(
  tolower(trimws(yedoma)) == "yedoma",
  "Yedoma", "Non-yedoma")]

rows_yedoma <- rbindlist(list(
  header_row("Yedoma"),
  make_row("Yedoma", "Yedoma",
           pld[yedoma_class == "Yedoma", poly_area],
           total_pf_area, total_pf_n, global_area, global_n)))

# ──────────────────────────────────────────────────────────────
# ── 7. Canadian Shield ────────────────────────────────────────
# ──────────────────────────────────────────────────────────────

row_shield <- make_row("Canadian Shield", "Canadian Shield",
                       pld[shield == TRUE, poly_area],
                       total_pf_area, total_pf_n, global_area, global_n)
row_shield[, Category    := "Canadian Shield"]
row_shield[, Subcategory := NA_character_]

# ──────────────────────────────────────────────────────────────
# ── 8. Totals rows ────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────
row_pf_total <- data.table(
  Category                                   = "Entire northern permafrost domain",
  Subcategory                                = NA_character_,
  `Total Lake area (km2)`                    = round(total_pf_area, 0),
  `% of lake area within northern PF region` = 100.0,
  `% of global lake area`                    = round(total_pf_area / global_area * 100, 1),
  `Median lake area (km2)`                   = round(median(pld$poly_area, na.rm = TRUE), 3),
  `Lake area IQR (km2)`                      = round(IQR(pld$poly_area,    na.rm = TRUE), 3),
  `Total number of lakes`                    = total_pf_n,
  `% of lakes within northern PF region`     = NA_real_,
  `% of global lakes`                        = round(total_pf_n / global_n * 100, 1))

row_global <- data.table(
  Category                                   = "Globally",
  Subcategory                                = NA_character_,
  `Total Lake area (km2)`                    = global_area,
  `% of lake area within northern PF region` = NA_real_,
  `% of global lake area`                    = NA_real_,
  `Median lake area (km2)`                   = 0.036,
  `Lake area IQR (km2)`                      = 0.076,
  `Total number of lakes`                    = global_n,
  `% of lakes within northern PF region`     = NA_real_,
  `% of global lakes`                        = NA_real_)

# ──────────────────────────────────────────────────────────────
# ── assemble ──────────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────
table1 <- rbindlist(list(
  rows_pf,
  rows_glacial,
  rows_tk,
  rows_ground_ice,
  rows_biomes,
  rows_yedoma,
  row_shield,
  row_pf_total,
  row_global), fill = TRUE)

# ──────────────────────────────────────────────────────────────
# ── join land area ────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────
# Match on Subcategory first; fall back to Category for header/total rows
lookup_key <- fifelse(!is.na(table1$Subcategory), table1$Subcategory, table1$Category)
table1[, `Land area (km2)` := land[match(lookup_key, Category), land_area_km2]]
table1[Category == "Globally", `Land area (km2)` := 134573595]

# Reorder so Land area appears after Subcategory
setcolorder(table1, c("Category", "Subcategory", "Land area (km2)"))

# ──────────────────────────────────────────────────────────────
# ── write ─────────────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────
fwrite(table1, lake_dist_table1_descending.csv')
