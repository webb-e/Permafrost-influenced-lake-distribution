library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(patchwork)
library(scales)

df<-fread('PF_x_Glaciation.csv')

#### helpers for later
order_pf <- c("Continuous", "Discontinuous", "Sporadic", "Isolated")

pf_map <- c(C = "Continuous", D = "Discontinuous", S = "Sporadic", I = "Isolated")

region_map <- c(glac_sh="Glaciated · Shield", glac_nsh="Glaciated · Non-shield",ungl="Unglaciated")

## ── 1. PIVOT wide → long ──────────────────────────────────────────────────────

all_long <- df %>% as_tibble() %>%
  select(grid_id, matches("^(dens|frac)_(glac_sh|glac_nsh|ungl)_[CDSI]$")) %>%
  pivot_longer(cols  = -grid_id,names_to  = c("metric", "region_code", "pf_code"),
    names_pattern = "^(dens|frac)_(glac_sh|glac_nsh|ungl)_([CDSI])$") %>%
  filter(!is.na(value)) %>%
  mutate(
    # dens_* stored as lakes/m²      → convert to lakes/100 km² (* 1e8)
    # frac_* stored as fractions 0-1 → convert to percent (* 100)
    value      = if_else(metric == "dens", value * 1e8, value * 100),
    variable   = if_else(metric == "dens", "Lake density (lakes/100 km²)", "Lake coverage (%)"),
    permafrost = factor(pf_map[pf_code], levels = order_pf),
    region     = factor(region_map[region_code],
                        levels = c("Glaciated · Shield", "Glaciated · Non-shield", "Unglaciated")),
    variable   = factor(variable, levels = c("Lake density (lakes/100 km²)", "Lake coverage (%)")) ) %>%
  select(grid_id, variable, region, permafrost, value)

## ── 2. SUMMARY STATS ─────────────────────────────────────────────────────────

sum_stats <- all_long %>%
  group_by(variable, region, permafrost) %>%
  summarise(median = median(value, na.rm = TRUE),
    q1     = quantile(value, 0.25, na.rm = TRUE),
    q3     = quantile(value, 0.75, na.rm = TRUE),.groups = "drop")

sum_dens <- sum_stats %>% filter(variable == "Lake density (lakes/100 km²)")
sum_frac <- sum_stats %>% filter(variable == "Lake coverage (%)")

## ── 3. FIGURE ────────────────────────────────────────────────────────────────

region_colors <- c("Glaciated · Shield"= "#2166ac","Glaciated · Non-shield" = "#2166ac", "Unglaciated"  = "black")
region_lty <- c("Glaciated · Shield"= "dashed","Glaciated · Non-shield" = "dotted","Unglaciated" = "solid")
region_lwd <- c( "Glaciated · Shield" = 0.7, "Glaciated · Non-shield" = 0.7, "Unglaciated" = 1.1)

base_theme <- theme_bw(base_size = 15) +
  theme(legend.position = "none",
    axis.text.y     = element_text(margin = margin(r = 10)),
    axis.text.x     = element_text(angle = -30, hjust = 0, vjust = 1),
    plot.margin     = margin(5, 75, 5, 5))

build_panel <- function(df, y_breaks, y_labeller, y_title, y_limits = c(NA, NA)) {
  labels_df <- df %>%
    filter(permafrost == tail(levels(permafrost), 1)) %>%
    mutate(label = recode(as.character(region),
                          "Glaciated · Shield"     = "Glaciated\nShield",
                          "Glaciated · Non-shield" = "Glaciated\nNon-shield",
                          "Unglaciated"            = "Unglaciated"))
  ggplot(df, aes(x = permafrost, y = median,
                 color = region, fill = region, group = region)) +
    geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.18, color = NA, na.rm = TRUE) +
    geom_line(aes(linetype = region, linewidth = region)) +
    geom_point(size = 2.2) +
    geom_text( data  = labels_df, aes(label  = label),
               hjust = 0, nudge_x = 0.3,size= 4.5, lineheight = 0.9, show.legend = FALSE) +
    scale_x_discrete(expand = expansion(add = c(0.1, 0))) +
    scale_y_log10(breaks = y_breaks, labels = y_labeller) +
    scale_color_manual(values     = region_colors, guide = "none") +
    scale_fill_manual(values      = region_colors, guide = "none") +
    scale_linetype_manual(values  = region_lty,    guide = "none") +
    scale_linewidth_manual(values = region_lwd,    guide = "none") +
    coord_cartesian(xlim = c(0.9, 4.2), ylim = y_limits, clip = "off") +
    labs(x = NULL, y = y_title) +
    base_theme
}

p_dens <- build_panel(sum_dens,y_breaks   = c(0.1, 1, 10, 100), 
                      y_labeller = function(x) formatC(x, format = "g"),
                      y_title = "Lake density (lakes / 100 km²)", y_limits   = c(0.01, 100))

p_cov <- build_panel(sum_frac,y_breaks   = c(0.1, 1, 10),y_labeller = function(x) paste0(x, "%"),
               y_title = "Lake coverage (%)",y_limits   = c(0.01, NA))

fig4 <- p_dens | p_cov

ggsave('Fig4.png',
       fig4, width = 9.5, height = 4.8, dpi = 300)
