library(dplyr)
library(purrr)
library(data.table)
library(ggplot2)
library(forcats)

df <- fread('future_PLD.csv')

pf <- df %>%
  filter(EXTENT %in% c("C", "D", "S", "I"),
         biome_name %in% c("Tundra", "Boreal Forests/Taiga")) %>%
  mutate(
    # Current permafrost class from EXTENT (matching map code)
    modern_class = case_when(
      EXTENT == "C" ~ 4L,
      EXTENT == "D" ~ 3L,
      EXTENT == "S" ~ 2L,
      EXTENT == "I" ~ 1L ),
    future_class_2deg = round(pf_2degC),
    # Loss = permafrost extent decreases at 2°C
    pf_loss_2deg = future_class_2deg < modern_class,
    thermokarst_lakes_grp = case_when(
      thermokarst_lakes %in% c("High", "Very High") ~ "High/very high",
      thermokarst_lakes %in% c("Low", "Moderate")   ~ "Low/moderate",
      thermokarst_lakes == "None"                   ~ "None",
      TRUE ~ NA_character_ ),
    thermokarst_wetlands_grp = case_when(
      thermokarst_wetlands %in% c("High", "Very High") ~ "High/very high",
      thermokarst_wetlands %in% c("Low", "Moderate")   ~ "Low/moderate",
      thermokarst_wetlands == "None"                  ~ "None",
      TRUE ~ NA_character_ ),
    thermokarst_combined = case_when(
      thermokarst_lakes %in% c("High", "Very High") |
        thermokarst_wetlands %in% c("High", "Very High") ~ "High/very high",
      thermokarst_lakes %in% c("Low", "Moderate") |
        thermokarst_wetlands %in% c("Low", "Moderate")   ~ "Low/moderate",
      thermokarst_lakes == "None" & thermokarst_wetlands == "None" ~ "None",
      TRUE ~ NA_character_))

summarise_group <- function(data, var, order, category) {
  totals <- data %>% filter(!is.na(.data[[var]])) %>%
    count(Subcategory = .data[[var]], name = "n_total")
  
  data %>% filter(!is.na(.data[[var]])) %>%
    group_by(Subcategory = .data[[var]]) %>%
    summarise(n_loss = sum(pf_loss_2deg, na.rm = TRUE), .groups = "drop") %>%
    left_join(totals, by = "Subcategory") %>%
    mutate(
      pct_loss = 100 * n_loss / n_total,
      Category = category ) %>%
    mutate(Subcategory = factor(Subcategory, levels = order))
}

plot_data <- bind_rows(
  summarise_group(pf, "glaciated", c("glaciated","unglaciated"), "Glacial history") %>%
    mutate(Subcategory = fct_recode(Subcategory,
                                    "Glaciated"   = "glaciated",
                                    "Unglaciated" = "unglaciated")),
  summarise_group(pf, "CONTENT",    c("High","Medium","Low"),           "Ground ice"),
  summarise_group(pf, "biome_name", c("Boreal Forests/Taiga","Tundra"), "Biome"),
  summarise_group(pf, "thermokarst_combined",
                  c("High/very high","Low/moderate","None"), "Thermokarst combined"),
  summarise_group(pf, "thermokarst_lakes_grp",
                  c("High/very high","Low/moderate","None"), "Thermokarst lakes"),
  summarise_group(pf, "thermokarst_wetlands_grp",
                  c("High/very high","Low/moderate","None"), "Thermokarst wetlands")) %>%
  mutate(
    Category = case_when(
      Category == "Glacial history"      ~ "Glacial\nhistory",
      Category == "Ground ice"           ~ "Ground ice",
      Category == "Thermokarst combined" ~ "Thermokarst\ncombined",
      Category == "Thermokarst lakes"    ~ "Thermokarst\nlakes",
      Category == "Thermokarst wetlands" ~ "Thermokarst\nwetlands",
      TRUE ~ Category ),
    Category = factor(Category, levels = c(
      "Glacial\nhistory", "Ground ice", "Biome",
      "Thermokarst\ncombined", "Thermokarst\nlakes", "Thermokarst\nwetlands")),
    FillGroup = case_when(
      Category == "Thermokarst\ncombined"  ~ paste("Combined",  Subcategory),
      Category == "Thermokarst\nlakes"     ~ paste("Lakes",     Subcategory),
      Category == "Thermokarst\nwetlands"  ~ paste("Wetlands",  Subcategory),
      TRUE ~ as.character(Subcategory) ))

category_order <- levels(plot_data$Category)

subcat_order <- list(
  "Glacial\nhistory"      = c("Glaciated", "Unglaciated"),
  "Ground ice"            = c("High", "Medium", "Low"),
  "Biome"                 = c("Boreal Forests/Taiga", "Tundra"),
  "Thermokarst\ncombined" = c("High/very high", "Low/moderate", "None"),
  "Thermokarst\nlakes"    = c("High/very high", "Low/moderate", "None"),
  "Thermokarst\nwetlands" = c("High/very high", "Low/moderate", "None"))

within_group_gap  <- 0.75
between_group_gap <- 1.45

position_df <- { y_cursor <- 1
bind_rows(lapply(rev(category_order), function(cat) {
  subcats <- subcat_order[[cat]]
  out <- data.frame(
    Category        = cat,
    Subcategory_chr = subcats,
    y_pos           = y_cursor + (seq_along(subcats) - 1) * within_group_gap)
  y_cursor <<- max(out$y_pos) + between_group_gap
  out }))}

plot_data <- plot_data %>%
  mutate(Subcategory_chr = as.character(Subcategory)) %>%
  left_join(position_df, by = c("Category", "Subcategory_chr"))

separator_positions <- position_df %>%
  group_by(Category) %>%
  summarise(y_min = min(y_pos), y_max = max(y_pos), .groups = "drop") %>%
  arrange(y_min) %>%
  mutate(separator = (y_max + lead(y_min)) / 2) %>%
  pull(separator) %>%
  na.omit()

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

x_max <- max(plot_data$pct_loss, na.rm = TRUE)

category_label_df <- plot_data %>%
  group_by(Category) %>%
  summarise(y = max(y_pos), .groups = "drop") %>%
  mutate(y = case_when(Category %in% c("Biome", "Ground ice") ~ y + 0.35, TRUE ~ y))

category_x <- 75
subcat_x   <- 1.5
bar_height  <- 0.64

plot <- ggplot(plot_data, aes(fill = FillGroup, color = FillGroup)) +
  geom_hline(yintercept = separator_positions, color = "grey60", linewidth = 0.4, linetype = "dashed") +
  geom_rect(aes(xmin = 0, xmax = pct_loss,
                ymin = y_pos - bar_height / 2, ymax = y_pos + bar_height / 2),
            alpha = 0.8, linewidth = 0.6) +
  geom_text(data = plot_data, aes(x = subcat_x, y = y_pos, label = Subcategory_chr),
            hjust = 0, vjust = 0.5, size = 5, color = "black", inherit.aes = FALSE) +
  geom_text(data = category_label_df, aes(x = category_x, y = y, label = Category),
            hjust = 0, vjust = 0.5, fontface = "bold", size = 5, color = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = fill_cols, guide = "none") +
  scale_color_manual(values = fill_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.03)),
                     labels = function(x) paste0(x, "%"), breaks = seq(0, 100, by = 20)) +
  labs(x = "Percentage of lakes", y = NULL) +
  theme_bw(base_size = 15) +
  theme(
    panel.border    = element_rect(color = "grey50", fill = NA, linewidth = 0.5),
    panel.grid      = element_blank(),
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.title.x    = element_text(size = 16),
    axis.text.x     = element_text(size = 14),
    plot.margin     = margin(t = 10, r = 15, b = 10, l = 11))

ggsave(
  "/Users/elizabethwebb/Library/CloudStorage/GoogleDrive-webb.elizabeth.e@gmail.com/My Drive/PostDoc/Lake distribution/figures/barplot_fig6_2deg.png",
  plot, width = 6, height = 8, dpi = 300)
