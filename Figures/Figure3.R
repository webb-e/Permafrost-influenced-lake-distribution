library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rnaturalearth)

# ── File paths ───────────────────────────────────────────────────────────────
pf_path     <- ## Brown et al., 2002 (reprojected) 
glac_path   <- ## Batchelor et al., 2025 (reprojected)
biome_path  <- ## Olson et al., 2001 (reprojected)
tk_path     <-  ## Olefeldt et al., 2016 (reprojected)
yedoma_path <- ## Strauss et al., 2021 (reprojected)
shield_path <- ## Natural Resources Canada 2022 (reprojected)
domain_path <- "/latitude_50N_EPSG3575.shp"
out_path    <- "Figure3.png"

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
      sin(angle_rad),  cos(angle_rad)),
    nrow = 2)
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

# ── Permafrost (shared between extent + ground ice panels) ───────────────────
pf_raw <- st_read(pf_path, quiet = TRUE)
if (st_crs(pf_raw)$epsg != crs_use_epsg)
  pf_raw <- st_transform(pf_raw, crs_use_epsg)

# ── Pre-process datasets ─────────────────────────────────────────────────────

### Permafrost extent
pf_extent <- pf_raw %>%
  mutate(EXTENT = dplyr::case_when(
    EXTENT == "C" ~ "Continuous",
    EXTENT == "D" ~ "Discontinuous",
    EXTENT == "S" ~ "Sporadic",
    EXTENT == "I" ~ "Isolated",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(EXTENT)) %>%
  mutate(EXTENT = factor(EXTENT,
                         levels = c("Continuous","Discontinuous","Sporadic","Isolated"))) %>%
  group_by(EXTENT) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

pf_extent_r <- rotate_geom(pf_extent, rotation_angle)

### Glacial history
glac_raw <- st_read(glac_path, quiet = TRUE)
if (st_crs(glac_raw)$epsg != crs_use_epsg)
  glac_raw <- st_transform(glac_raw, crs_use_epsg)

glac_union <- st_union(glac_raw) %>% st_as_sf()
glac_union$glaciation <- "Glaciated"

land_union <- st_union(land_t) %>% st_as_sf()
unglac     <- st_difference(land_union, glac_union) %>% st_as_sf()
unglac$glaciation <- "Unglaciated"

glac_combined <- bind_rows(
  glac_union %>% rename(geometry = x),
  unglac     %>% rename(geometry = x)) %>%
  mutate(glaciation = factor(glaciation, levels = c("Glaciated","Unglaciated")))

glac_r <- rotate_geom(glac_combined, rotation_angle)

### Biomes
biome_raw <- st_read(biome_path, quiet = TRUE)
if (st_crs(biome_raw)$epsg != crs_use_epsg)
  biome_raw <- st_transform(biome_raw, crs_use_epsg)

biome_filt <- biome_raw %>%
  st_make_valid() %>%
  filter(BIOME %in% c(6, 11)) %>%
  mutate(BIOME_label = dplyr::case_when(
    BIOME == 6  ~ "Boreal Forest",
    BIOME == 11 ~ "Tundra")) %>%
  mutate(BIOME_label = factor(BIOME_label, levels = c("Boreal Forest","Tundra"))) %>%
  group_by(BIOME_label) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

biome_r <- rotate_geom(biome_filt, rotation_angle)

### Ground ice
pf_ice <- pf_raw %>%
  mutate(CONTENT = dplyr::case_when(
    CONTENT == "h" ~ "High",
    CONTENT == "m" ~ "Medium",
    CONTENT == "l" ~ "Low",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(CONTENT)) %>%
  mutate(CONTENT = factor(CONTENT, levels = c("High","Medium","Low"))) %>%
  group_by(CONTENT) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

pf_ice_r <- rotate_geom(pf_ice, rotation_angle)

### Thermokarst lakes
tk_raw <- st_read(tk_path, quiet = TRUE)
if (st_crs(tk_raw)$epsg != crs_use_epsg)
  tk_raw <- st_transform(tk_raw, crs_use_epsg)

tk_lakes <- tk_raw %>%
  st_make_valid() %>%
  mutate(TK_lakes = dplyr::case_when(
    TKThLP %in% c("Very High","High") ~ "High/Very High",
    TKThLP %in% c("Moderate","Low")   ~ "Low/Moderate",
    TKThLP == "None"                   ~ "None",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(TK_lakes)) %>%
  mutate(TK_lakes = factor(TK_lakes,
                           levels = c("High/Very High","Low/Moderate","None"))) %>%
  group_by(TK_lakes) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

tk_lakes_r <- rotate_geom(tk_lakes, rotation_angle)

### Thermokarst wetlands
tk_wetlands <- tk_raw %>%
  st_make_valid() %>%
  mutate(TK_wet = dplyr::case_when(
    TKHP %in% c("Very High","High") ~ "High/Very High",
    TKHP %in% c("Moderate","Low")   ~ "Low/Moderate",
    TKHP == "None"                   ~ "None",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(TK_wet)) %>%
  mutate(TK_wet = factor(TK_wet,
                         levels = c("High/Very High","Low/Moderate","None"))) %>%
  group_by(TK_wet) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

tk_wet_r <- rotate_geom(tk_wetlands, rotation_angle)

### Yedoma
yedoma_raw <- st_read(yedoma_path, quiet = TRUE)
if (st_crs(yedoma_raw)$epsg != crs_use_epsg)
  yedoma_raw <- st_transform(yedoma_raw, crs_use_epsg)

yedoma_r       <- rotate_geom(yedoma_raw, rotation_angle)
yedoma_r$label <- "Yedoma"

### Canadian Shield
shield_raw <- st_read(shield_path, quiet = TRUE)
if (st_crs(shield_raw)$epsg != crs_use_epsg)
  shield_raw <- st_transform(shield_raw, crs_use_epsg)

shield_r       <- rotate_geom(shield_raw, rotation_angle)
shield_r$label <- "Canadian shield"

# ── Shared theme & coord ─────────────────────────────────────────────────────
base_theme <- function() {
  theme_void(base_size = 20) +
    theme(panel.background  = element_rect(fill = "grey18", color = NA),
          plot.background   = element_rect(fill = "grey18", color = NA),
          legend.position   = c(0.03, 0.05),
          legend.justification = c(0, 0),
          panel.border = element_rect(color = "white", fill = NA, linewidth = 0.5),
          legend.background = element_rect(fill = alpha("grey18", 0.7), color = NA),
          legend.title      = element_text(color = "white", size = 25, face = "bold"),
          legend.text       = element_text(color = "white", size = 22),
          legend.key        = element_rect(fill = "grey18", color = NA),
          plot.margin       = margin(10, 10, 10, 10)) 
}

shared_coord <- function() {
  coord_sf(datum = NA, xlim = XLIM, ylim = YLIM, expand = FALSE)
}

land_layer <- function() {
  geom_sf(data = land_r, fill = LAND_GREY, color = NA, inherit.aes = FALSE)
}

# ── Panels ───────────────────────────────────────────────────────────────────

### C – Permafrost extent
pf_cols <- c("Continuous"    = "#6a3d9a",
             "Discontinuous" = "#9e7cc5",
             "Sporadic"      = "#c8b0e0",
             "Isolated"      = "#ecdff5")

pC <- ggplot() +
  land_layer() +
  geom_sf(data = pf_extent_r, aes(fill = EXTENT),
          color = NA, alpha = 0.88, inherit.aes = FALSE) +
  scale_fill_manual(values = pf_cols, name = "Permafrost extent",
                    na.translate = FALSE,
                    guide = guide_legend(title.position = "top", ncol = 2,
                                         title.hjust = 0.5)) +
  shared_coord() +
  base_theme()

### D – Glacial history
glac_cols <- c("Glaciated"   = "#2166ac",
               "Unglaciated" = LAND_GREY)

 # scale_fill_manual(values = c("Yedoma" = "#99AFD7"), name = "Yedoma domain",

pD <- 
  ggplot() +
  land_layer() +
  geom_sf(data = glac_r, aes(fill = glaciation),
          color = NA, alpha = 0.85, inherit.aes = FALSE) +
  scale_fill_manual(values = glac_cols, name = "Glacial history",
                    na.translate = FALSE,
                    guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
  shared_coord() +
  base_theme()

### E – Biomes
biome_cols <- c("Boreal Forest" = "#5c7a3e",
                "Tundra"        = "#d4b483")

pE <- ggplot() +
  land_layer() +
  geom_sf(data = biome_r, aes(fill = BIOME_label),
          color = NA, alpha = 0.85, inherit.aes = FALSE) +
  scale_fill_manual(values = biome_cols, name = "Biome",
                    na.translate = FALSE,
                    guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
  shared_coord() +
  base_theme()

### F – Ground ice content
ice_cols <- c("High"   = "#c0392b",
              "Medium" = "#e8927c",
              "Low"    = "#f9d4c8")

pF <- ggplot() +
  land_layer() +
  geom_sf(data = pf_ice_r, aes(fill = CONTENT),
          color = NA, alpha = 0.85, inherit.aes = FALSE) +
  scale_fill_manual(values = ice_cols, name = "Ground ice content",
                    na.translate = FALSE,
                    guide = guide_legend(title.position = "top", ncol = 2,
                                         title.hjust = 0.5)) +
  shared_coord() +
  base_theme()

### G – Thermokarst lakes
tk_lake_cols <- c("High/Very High" = "#d73027",
                  "Low/Moderate"   = "#fc8d59",
                  "None"           = "#fee08b")

pG <- ggplot() +
  land_layer() +
  geom_sf(data = tk_lakes_r, aes(fill = TK_lakes),
          color = NA, alpha = 0.85, inherit.aes = FALSE) +
  scale_fill_manual(values = tk_lake_cols, name = "Thermokarst lakes",
                    na.translate = FALSE,
                    guide = guide_legend(title.position = "top", ncol = 2,
                                         title.hjust = 0.5)) +
  shared_coord() +
  base_theme()

### H – Thermokarst wetlands
tk_wet_cols <- c("High/Very High" = "#810f7c",
                 "Low/Moderate"   = "#df65b0",
                 "None"           = "#f1eef6")

pH <- ggplot() +
  land_layer() +
  geom_sf(data = tk_wet_r, aes(fill = TK_wet),
          color = NA, alpha = 0.85, inherit.aes = FALSE) +
  scale_fill_manual(values = tk_wet_cols, name = "Thermokarst wetlands",
                    na.translate = FALSE,
                    guide = guide_legend(title.position = "top", ncol = 2,
                                         title.hjust = 0.5)) +
  shared_coord() +
  base_theme()

### I – Yedoma domain
pI <- ggplot() +
  land_layer() +
  geom_sf(data = yedoma_r,
          aes(fill = label),
          color = NA, alpha = 0.85, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Yedoma" = "#9ecae1"), name = NULL,
                    guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
  shared_coord() +
  base_theme() +
  theme(legend.position = c(0.01, 0.10))

### J – Canadian Shield
pJ <- ggplot() +
  land_layer() +
  geom_sf(data = shield_r,
          aes(fill = label),
          color = NA, alpha = 0.85, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Canadian shield" = "#35978f"), name = NULL,
                    guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
  shared_coord() +
  base_theme() +
  theme(legend.position = c(0.01, 0.10))

# ── Assemble  layout & save ─────────────────────────────────────────────

blank <- plot_spacer() & theme(plot.background = element_rect(fill = "grey18", color = NA))

fig2 <- wrap_plots(pC, pD, pE, pF, pG, pH, pI, pJ, blank, ncol = 3) +
  plot_layout(guides = "keep") &
  theme(plot.margin = margin(5, 10, 5, 10))
ggsave(out_path, plot = fig2, bg = "grey18",
       height = 22, width = 23.45, units = "in", dpi = 300)
