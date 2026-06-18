library(data.table)
library(ggplot2)
library(scales)
library(patchwork)

#--------------------------
#---------- read in data
#--------------------------
## created from global_PF_extent.py
dt<-fread('PLD_PF_global.csv')

#--------------------------
#---------- wrangle data into usable format
#--------------------------
dt[, permafrost := fifelse(EXTENT %in% c("C","c"), "Continuous PF",
                           fifelse(EXTENT %in% c("D","d"), "Discontinuous PF",
                                   fifelse(EXTENT %in% c("S","s"), "Sporadic PF",
                                           fifelse(EXTENT %in% c("I","i"), "Isolated PF",
                                                   fifelse(EXTENT %in% c("none","None","NONE", ""), "No permafrost", as.character(EXTENT))))))]

# latitude bins 
bin <- 1
breaks <- seq(-50, 81, by = bin)
dt[, lat_bin := bin * floor(lat / bin)]
dt <- dt[lat_bin >= min(breaks) & lat_bin <= max(breaks)]

# grouping
dt[, lake_group := fifelse(on_glacier == 1, "Glacial",
                           fifelse(permafrost == "No permafrost",
                                   "Permafrost-free",
                                   "Permafrost-influenced"))]

# Exclude glacial lakes
dt <- dt[on_glacier != 1]

# Desired order (legend)
dt[, lake_group := factor(lake_group,
                          levels = c("Permafrost-free",
                                     "Permafrost-influenced"))]

## aggregate
agg <- dt[, .(n_lakes  = .N, 
              area_km2 = sum(poly_area, na.rm = TRUE)),
          by = .(lat_bin, lake_group)]
agg[, lat_f := factor(lat_bin, levels = sort(unique(lat_bin)))]

### colors
cols <- c("Permafrost-free"       = "#212e53",
          "Permafrost-influenced" = "#bed3c3")

#--------------------------
#---------- plots
#--------------------------
bs <- 16
### area plot
p_area <- ggplot(agg, aes(x = area_km2, y = lat_f, fill = lake_group, color = lake_group)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = cols, breaks = names(cols)) +
  guides(color = "none") +
  scale_color_manual(values = cols) +
  scale_x_continuous(labels = comma) +
  scale_y_discrete(breaks = as.character(seq(-50, 80, by = 10)),
                   labels = seq(-50, 80, by = 10)) +
  labs(x = "Lake area (km²)", y = "Latitude", fill = NULL) +
  theme_bw(base_size = bs) +
  theme(panel.grid = element_blank(), legend.position = "none")

## count plot
p_count <- ggplot(agg, aes(x = n_lakes, y = lat_f, fill = lake_group, color = lake_group)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  guides(color = "none") +
  scale_fill_manual(values = cols, breaks = names(cols)) +
  scale_color_manual(values = cols) +
  scale_x_continuous(labels = comma) +
  scale_y_discrete(breaks = as.character(seq(-50, 80, by = 10)),
                   labels = seq(-50, 80, by = 10)) +
  labs(x = "Number of lakes", y = NULL, fill = NULL) +
  theme_bw(base_size = bs) +
  theme(panel.grid = element_blank(), legend.position = "none")

finalplot <- (p_area | p_count) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top", legend.title = element_blank())

finalplot

#--------------------------
#---------- save
#--------------------------
ggsave(path = 'Lake distribution/figures',
      filename = "LakeDist_Fig1.png",plot = finalplot,
      width = 7, height = 6, units = "in", dpi = 500)

#--------------------------
#---------- percentages permafrost influenced/not
#--------------------------

# Total area and count by lake_group
summary_by_group <- dt[, .(n_lakes = .N,
                           total_area = sum(poly_area, na.rm = TRUE)),
                          by = lake_group]

# Global totals
total_lakes <- dt[, .N]
total_area  <- dt[, sum(poly_area, na.rm = TRUE)]

# Add percentages
summary_by_group[, `:=`(
  pct_lakes = round(n_lakes    / total_lakes * 100, 1),
  pct_area  = round(total_area / sum(total_area) * 100, 1))]

# Total area and count by lake_group
summary_by_group <- dt[, .(
  n_lakes    = .N,
  total_area = sum(poly_area, na.rm = TRUE)), by = lake_group]

# Global totals
total_lakes <- dt[, .N]
total_area  <- dt[, sum(poly_area, na.rm = TRUE)]

# Add percentages
summary_by_group[, `:=`(
  pct_lakes = round(n_lakes    / total_lakes * 100, 1),
  pct_area  = round(total_area / sum(total_area) * 100, 1))]

# Pull the permafrost-influenced numbers
pf <- summary_by_group[lake_group == "Permafrost-influenced"]
cat(sprintf("%.1f%% of global lake area and %.1f%% of lakes are in regions affected by permafrost\n",
            pf$pct_area, pf$pct_lakes))
