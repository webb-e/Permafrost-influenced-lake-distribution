library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)

# ── File paths ───────────────────────────────────────────────────────────────
future_pf_path <- 'future_PLD.gpkg'
domain_path    <- "latitude_50N_EPSG3575.shp"
out_path       <- "futurePFmap_2deg.png"

# ── Parameters ───────────────────────────────────────────────────────────────
rotation_angle <- 10
crs_use_epsg   <- 3575
XLIM           <- c(-4600000, 4600000)
YLIM           <- c(-4400000, 4400000)
BBOX_CROP      <- c(xmin = -5500000, ymin = -5500000, xmax = 5500000, ymax = 5500000)
LAND_GREY      <- "grey40"

# ── Helpers ──────────────────────────────────────────────────────────────────
rotate_geom <- function(sf_obj, angle_deg) {
  angle_rad  <- angle_deg * pi / 180
  rot_matrix <- matrix(c(cos(angle_rad), -sin(angle_rad),
                         sin(angle_rad),  cos(angle_rad)), nrow = 2)
  st_geometry(sf_obj) <- st_geometry(sf_obj) * rot_matrix
  sf_obj
}

sf::sf_use_s2(FALSE)

# ════════════════════════════════════════════════════════════════════════════
# PREPARE BASEMAPS
# ════════════════════════════════════════════════════════════════════════════
land   <- ne_download(scale = "medium", type = "land", category = "physical", returnclass = "sf")
land_t <- st_transform(land, crs_use_epsg)
land_t <- st_crop(land_t, st_bbox(BBOX_CROP, crs = st_crs(crs_use_epsg)))
land_r <- rotate_geom(land_t, rotation_angle)

domain_raw <- st_read(domain_path, quiet = TRUE)
if (st_crs(domain_raw)$epsg != crs_use_epsg)
  domain_raw <- st_transform(domain_raw, crs_use_epsg)
domain_r <- rotate_geom(domain_raw, rotation_angle)

# ════════════════════════════════════════════════════════════════════════════
# LOAD & PROCESS FUTURE PERMAFROST LAKES (2°C only)
# ════════════════════════════════════════════════════════════════════════════
pf_lakes_r <- st_read(future_pf_path, quiet = FALSE,
                      query = "SELECT pf_category_2deg, geom FROM future_PLD_rotated") %>%
  st_set_crs(NA) %>%
  mutate(pf_category = factor(pf_category_2deg, levels = c(
    "Stable (no change at 2°C)",
    "Reduction of continuous permafrost",
    "Reduction of discontinuous permafrost",
    "Reduction of sporadic permafrost",
    "Reduction of isolated permafrost")))

# ════════════════════════════════════════════════════════════════════════════
# COLOR PALETTE
# ════════════════════════════════════════════════════════════════════════════
pf_future_cols <- c(
  "Reduction of continuous permafrost"    = "#f0e442",  
  "Reduction of discontinuous permafrost" = "#e69500",   
  "Reduction of sporadic permafrost"      = "#cc3399",   
  "Reduction of isolated permafrost"      = "#8b0000",   
  "Stable (no change at 2°C)"        = "#56b4e9"    )

# ════════════════════════════════════════════════════════════════════════════
# BUILD PLOT
# ════════════════════════════════════════════════════════════════════════════
map_p <- ggplot() +
  geom_sf(data = land_r, fill = LAND_GREY, color = NA, inherit.aes = FALSE) +
  geom_sf(data = pf_lakes_r, aes(color = pf_category),
          fill = NA, show.legend = FALSE, linewidth = 0.03, inherit.aes = FALSE) +
  # Invisible point layer to drive the legend
  geom_point(
    data = data.frame(
      x = rep(0, 5), y = rep(0, 5),
      pf_category = factor(levels(pf_lakes_r$pf_category),
                           levels = levels(pf_lakes_r$pf_category))),
    aes(x = x, y = y, color = pf_category), size = 0, inherit.aes = FALSE) +
  geom_sf(data = domain_r, fill = NA, color = "white",
          linewidth = 0.6, linetype = "dashed", inherit.aes = FALSE) +
  scale_color_manual(values = pf_future_cols,
    name   = "Permafrost loss\nat 2°C warming",
    na.translate = FALSE,
    guide  = guide_legend(
      title.position = "top", title.hjust = 0,
      override.aes   = list(shape = 16, size = 6, linewidth = 0))) +
  theme_void() +
  theme(panel.background  = element_rect(fill = "grey18", color = NA),
    panel.border      = element_rect(color = "white", fill = NA, linewidth = 0.3),
    plot.background   = element_rect(fill = "grey18", color = NA),
    legend.position      = c(0.03, 0.05),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = alpha("grey18", 0.7), color = NA),
    legend.title      = element_text(color = "white", size = 20, face = "bold"),
    legend.text       = element_text(color = "white", size = 16),
    legend.key        = element_rect(fill = "grey18", color = NA),
    plot.margin       = margin(10, 10, 10, 10) ) +
  coord_sf(datum = NA, xlim = XLIM, ylim = YLIM, expand = FALSE)

# ════════════════════════════════════════════════════════════════════════════
# SAVE
# ════════════════════════════════════════════════════════════════════════════
ggsave(out_path, plot = map_p,
       width = 12, height = 12, units = "in", dpi = 300, bg = "grey18")
