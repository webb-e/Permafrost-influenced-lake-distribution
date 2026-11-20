library(dplyr)
library(purrr)
library(data.table)
library(ggplot2)
library(forcats)

df <- fread('future_PLD.csv')
## ────────────────────────────────────────────────────────────────
## ── DATA WRANGLING ──────────────────────────────────────────────
## ────────────────────────────────────────────────────────────────

# Group thermokarst classes
group_tk <- function(x) {
  x <- as.character(x)
  case_when(
    x %in% c("High", "Very High") ~ "High/very high",
    x %in% c("Low", "Moderate")   ~ "Low/moderate",
    TRUE                          ~ "None")
}

## Recode PF based on future projections (which lakes experience a reduction?)
pf <- df %>%
  filter(EXTENT %in% c("C", "D", "S", "I")) %>%     
  mutate(modern_class = case_when(
      EXTENT == "C" ~ 4L,
      EXTENT == "D" ~ 3L,
      EXTENT == "S" ~ 2L,
      EXTENT == "I" ~ 1L),
    pf_loss_2deg_min  = pf_2degC_min  < modern_class,
    pf_loss_2deg_mean = pf_2degC_mean < modern_class,
    pf_loss_2deg_max  = pf_2degC_max  < modern_class,
    thermokarst_lakes_grp    = group_tk(thermokarst_lakes),
    thermokarst_wetlands_grp = group_tk(thermokarst_wetlands),
    thermokarst_combined = case_when(
      thermokarst_lakes_grp == "High/very high" |
        thermokarst_wetlands_grp == "High/very high" ~ "High/very high",
      thermokarst_lakes_grp == "Low/moderate" |
        thermokarst_wetlands_grp == "Low/moderate"   ~ "Low/moderate",
      TRUE                                           ~ "None"))

# Check for missing projections
colSums(is.na(pf[, c("pf_2degC_min", "pf_2degC_mean", "pf_2degC_max")]))

## ────────────────────────────────────────────────────────────────
## ── DATA SUMMARIES ────────────────────────────────
## ────────────────────────────────────────────────────────────────
summarise_group <- function(data, var, category) {
  data %>% filter(!is.na(.data[[var]])) %>%
    group_by(Subcategory = as.character(.data[[var]])) %>%
    summarise(
      n_total         = n(),
      # percent of lakes with a valid projection that lose permafrost
      pct_loss_minmap = 100 * mean(pf_loss_2deg_min,  na.rm = TRUE),
      pct_loss        = 100 * mean(pf_loss_2deg_mean, na.rm = TRUE),
      pct_loss_maxmap = 100 * mean(pf_loss_2deg_max,  na.rm = TRUE),
      .groups = "drop") %>%
    mutate(
      pct_lo   = pmin(pct_loss_minmap, pct_loss_maxmap, pct_loss),
      pct_hi   = pmax(pct_loss_minmap, pct_loss_maxmap, pct_loss),
      Category = category)
}

plot_data <- bind_rows(
  summarise_group(pf, "glaciated", "Glacial history") %>%
    mutate(Subcategory = recode(Subcategory,
                                "glaciated"   = "Glaciated",
                                "unglaciated" = "Unglaciated")),
  summarise_group(pf, "CONTENT", "Ground ice"),
  summarise_group(pf %>% filter(biome_name %in% c("Tundra", "Boreal Forests/Taiga")),
                  "biome_name", "Biome"),
  summarise_group(pf, "thermokarst_combined",     "Thermokarst combined"),
  summarise_group(pf, "thermokarst_lakes_grp",    "Thermokarst lakes"),
  summarise_group(pf, "thermokarst_wetlands_grp", "Thermokarst wetlands")) %>%
  filter(Subcategory %in% c("Glaciated", "Unglaciated", "High", "Medium", "Low",
                            "Boreal Forests/Taiga", "Tundra",
                            "High/very high", "Low/moderate", "None")) %>%
  mutate( Category = recode(Category,
                      "Glacial history"      = "Glacial\nhistory",
                      "Thermokarst combined" = "Thermokarst\ncombined",
                      "Thermokarst lakes"    = "Thermokarst\nlakes",
                      "Thermokarst wetlands" = "Thermokarst\nwetlands"),
    Category = factor(Category, levels = c(
      "Glacial\nhistory", "Ground ice", "Biome",
      "Thermokarst\ncombined", "Thermokarst\nlakes", "Thermokarst\nwetlands")),
    FillGroup = case_when(
      Category == "Thermokarst\ncombined" ~ paste("Combined", Subcategory),
      Category == "Thermokarst\nlakes"    ~ paste("Lakes",    Subcategory),
      Category == "Thermokarst\nwetlands" ~ paste("Wetlands", Subcategory),
      TRUE ~ Subcategory),
    Subcat_f = factor(Subcategory, levels = c(
      "Glaciated", "Unglaciated",
      "High", "Medium", "Low",
      "Boreal Forests/Taiga", "Tundra",
      "High/very high", "Low/moderate", "None")))

## ────────────────────────────────────────────────────────────────
## ── FIGURE ──────────────────────────────────────────────────────
## ────────────────────────────────────────────────────────────────

map_labels <- c("2°C min"  = "Minimum projected permafrost extent",
                "2°C mean" = "Mean projected permafrost extent",
                "2°C max"  = "Maximum projected permafrost extent")
fill_cols <- c(
  "Glaciated"   = "#2166ac",
  "Unglaciated" = "grey39",
  "High"   = "#c0392b",
  "Medium" = "#e8927c",
  "Low"    = "#f9d4c8",
  "Boreal Forests/Taiga" = "#5c7a3e",
  "Tundra"               = "#d4b483",
  "Combined High/very high" = "#01665e",
  "Combined Low/moderate"   = "#5ab4ac",
  "Combined None"           = "#c7eae5",
  "Lakes High/very high" = "#d73027",
  "Lakes Low/moderate"   = "#fc8d59",
  "Lakes None"           = "#fee08b",
  "Wetlands High/very high" = "#810f7c",
  "Wetlands Low/moderate"   = "#df65b0",
  "Wetlands None"           = "#dadaeb")

pts <- bind_rows(
  plot_data %>% transmute(Category, Subcat_f, FillGroup, pct = pct_loss_minmap, Map = "2°C min"),
  plot_data %>% transmute(Category, Subcat_f, FillGroup, pct = pct_loss,        Map = "2°C mean"),
  plot_data %>% transmute(Category, Subcat_f, FillGroup, pct = pct_loss_maxmap, Map = "2°C max")) %>%
  mutate(Map = factor(Map, levels = c("2°C max", "2°C mean", "2°C min")))

overall_mean <- 100 * mean(pf$pf_loss_2deg_mean, na.rm = TRUE)

plot <- ggplot() +
  geom_vline(xintercept = overall_mean, linetype = "solid", color = "grey40", linewidth = 0.4) +
  geom_segment(data = plot_data,
               aes(x = pct_lo, xend = pct_hi, y = Subcat_f, yend = Subcat_f, color = FillGroup),
               linewidth = 0.7) +
  geom_point(data = pts,
             aes(x = pct, y = Subcat_f, fill = FillGroup, shape = Map, size = Map),
             color = "black", stroke = 0.4) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = fill_cols, guide = "none") +
  scale_color_manual(values = fill_cols, guide = "none") +
  scale_shape_manual(values = c("2°C min" = 24, "2°C mean" = 21, "2°C max" = 22),
                     labels = map_labels, name = NULL) +
  scale_size_manual(values = c("2°C min" = 5, "2°C mean" = 7, "2°C max" = 5),
                    labels = map_labels, name = NULL) +
  guides(shape = guide_legend(override.aes = list(fill = "grey70"))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02)),
                     labels = function(x) paste0(x, "%"), breaks = seq(0, 100, by = 20)) +
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = "Percentage of lakes", y = NULL) +
  theme_bw(base_size = 25) +
  theme(
    panel.border       = element_rect(color = "grey50", fill = NA, linewidth = 0.5),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.spacing      = unit(0.3, "lines"),
    strip.background   = element_blank(),
    strip.text.y       = element_text(angle = 0, hjust = 0, face = "bold", size = 14),
    legend.position    = "bottom",
    legend.direction = "vertical",
    plot.margin        = margin(t = 10, r = 10, b = 10, l = 10))

## ────────────────────────────────────────────────────────────────
## ── SAVE FIGURE ─────────────────────────────────────────────────
## ────────────────────────────────────────────────────────────────
ggsave("Fig7.png",plot, width = 10, height = 9, dpi = 500)

## ────────────────────────────────────────────────────────────────
## ── VALUES FOR TEXT ─────────────────────────────────────────────
## ────────────────────────────────────────────────────────────────
loss_summary <- data.frame(
  map = c("2°C min", "2°C mean", "2°C max"),
  n_lakes_loss = c(sum(pf$pf_loss_2deg_min,  na.rm = TRUE),
                   sum(pf$pf_loss_2deg_mean, na.rm = TRUE),
                   sum(pf$pf_loss_2deg_max,  na.rm = TRUE)),
  area_loss_km2 = c(sum(pf$poly_area[pf$pf_loss_2deg_min  %in% TRUE]),
                    sum(pf$poly_area[pf$pf_loss_2deg_mean %in% TRUE]),
                    sum(pf$poly_area[pf$pf_loss_2deg_max  %in% TRUE]))) %>%
  mutate(pct_lakes = 100 * n_lakes_loss  / nrow(pf),
         pct_area  = 100 * area_loss_km2 / sum(pf$poly_area, na.rm = TRUE))

loss_summary
