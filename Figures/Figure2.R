library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(viridis)

# ── File paths ───────────────────────────────────────────────────────────────
gpkg_lakes  <- 'PF_lake_density_coverage.gpkg'
domain_path <- 'latitude_50N_EPSG3575.shp"
out_path    <- "Lake_density_Fig2.png"

# ── Parameters ───────────────────────────────────────────────────────────────

rotation_angle <- 10
crs_use_epsg   <- 3575
XLIM           <- c(-5000000, 5000000)
YLIM           <- c(-4800000, 4800000)
BBOX_CROP      <- c(xmin = -5500000, ymin = -5500000, xmax = 5500000, ymax = 5500000)
LAND_GREY <- "grey40"   

# ── Helpers ──────────────────────────────────────────────────────────────────
rotate_geom <- function(sf_obj, angle_deg) {
  angle_rad  <- angle_deg * pi / 180
  rot_matrix <- matrix(
    c(cos(angle_rad), -sin(angle_rad),
      sin(angle_rad),  cos(angle_rad)),nrow = 2)
  st_geometry(sf_obj) <- st_geometry(sf_obj) * rot_matrix
  sf_obj
}

sf::sf_use_s2(FALSE)

# ── Basemap ──────────────────────────────────────────────────────────────────
land <- ne_download(scale = "medium", type = "land",
                    category = "physical", returnclass = "sf")
land_t <- st_transform(land, crs_use_epsg)
land_t <- st_crop(land_t, st_bbox(BBOX_CROP, crs = st_crs(crs_use_epsg)))
land_r <- rotate_geom(land_t, rotation_angle)

# ── Study domain (50°N line) ─────────────────────────────────────────────────
domain_raw <- st_read(domain_path, quiet = TRUE)
if (st_crs(domain_raw)$epsg != crs_use_epsg)
  domain_raw <- st_transform(domain_raw, crs_use_epsg)
domain_r <- rotate_geom(domain_raw, rotation_angle)

# ── Lake data ────────────────────────────────────────────────────────────────
lakes_raw <- st_read(gpkg_lakes, quiet = TRUE)
if (st_crs(lakes_raw)$epsg != crs_use_epsg)
  lakes_raw <- st_transform(lakes_raw, crs_use_epsg)

lakes_frac_r <- lakes_raw %>%
  filter(!is.na(lake_fraction_pf_area)) %>%
  rotate_geom(rotation_angle)

lakes_dens_r <- lakes_raw %>%
  filter(!is.na(lake_density_pf_area)) %>%
  rotate_geom(rotation_angle)

# ── Shared theme & coord ─────────────────────────────────────────────────────
base_theme <- function() {
  theme_minimal(base_size = 20) +
    theme(
      panel.background     = element_rect(fill = "grey18", color = NA),
      panel.border         = element_rect(color = "white", fill = NA, linewidth = 0.3),
      axis.text            = element_blank(),
      axis.ticks           = element_blank(),
      panel.grid           = element_blank(),
      plot.title           = element_blank(),
      plot.margin          = margin(2, 2, 2, 2),
      legend.position      = c(0.01, 0.01),
      legend.justification = c(0, 0),
      legend.direction     = "vertical",
      legend.background = element_rect(fill = alpha("grey18", 0.7), color = NA),
      legend.title      = element_text(color = "white", size = 20, face = "bold"),
      legend.text       = element_text(color = "white", size = 18),
      legend.margin        = margin(4, 7, 4, 7),
      legend.key.size      = unit(0.75, "cm"))
  }
     # legend.title         = element_text(face = "bold", size = 22),
     # legend.text          = element_text(size = 20))
     

shared_coord <- function() {
  coord_sf(datum = NA, xlim = XLIM, ylim = YLIM, expand = FALSE)}

land_layer <- function() {
  geom_sf(data = land_r, fill = LAND_GREY, color = NA, inherit.aes = FALSE)}

domain_layer <- function() {
  geom_sf(data = domain_r, fill = NA, color = "white",
          linewidth = 0.5, linetype = "dashed", inherit.aes = FALSE)}

# ── Panel A – Lake fraction ──────────────────────────────────────────────────
pA <- ggplot() +
  land_layer() +
  geom_sf(data = lakes_frac_r,
          aes(fill = lake_fraction_pf_area * 100),
          color = NA, linewidth = 0, inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "mako", direction = -1,
                       na.value = "transparent",
                       name = "Lake cover (%)",
                       trans = scales::pseudo_log_trans(base = 10),
                       breaks = c(0, 5, 25),
                       guide = guide_colorbar(title.position = "top",
                                              title.hjust = 0.5,
                                              direction = "horizontal",
                                              barwidth = 6, barheight = 1.22)) +
  domain_layer() +
  shared_coord() +
  base_theme()

# ── Panel B – Lake density ───────────────────────────────────────────────────
pB <- ggplot() +
  land_layer() +
  geom_sf(data = lakes_dens_r,
          aes(fill = lake_density_pf_area),
          color = NA, linewidth = 0, inherit.aes = FALSE) +
  scale_fill_distiller(palette = "GnBu", direction = 1,
                       na.value = "transparent",
                       breaks = c(0, 10, 100),
                       name = "Lake density\n(lakes/100 km²)",
                       trans = scales::pseudo_log_trans(base = 10),
                       guide = guide_colorbar(title.position = "top",
                                              title.hjust = 0.5,
                                              direction = "horizontal",
                                              barwidth = 6, barheight = 1.2)) +
  domain_layer() +
  shared_coord() +
  base_theme()

# ── Assemble & save ──────────────────────────────────────────────────────────
fig1 <- wrap_plots(pA, pB, ncol = 2) + # plot_annotation(tag_levels = "A")
        plot_layout(guides = "keep") &  
  theme( plot.background = element_rect(fill = "grey18", color = NA),
    plot.margin = margin(5, 10, 5, 10))
  #  plot.tag.position = c(0.03, 0.97),  # top-left
 #   plot.tag = element_text(size = 24, face = "bold", color = "white"))

ggsave(out_path, plot = fig1,       bg = "grey18",
       width = 14, height = 7, units = "in", dpi = 300)
