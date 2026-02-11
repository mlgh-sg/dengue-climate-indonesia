################################################################################
# Figure 5: Climate-Dengue Phase Lag Heatmap
#
# This script generates Figure 5 for the manuscript, showing:
#   - 5A: Phase lag heatmap (province x climate variable)
#   - 5B: Phase coherence heatmap (optional/supplementary)
#
# Dependencies: Requires data-processing.R to be run first
################################################################################

#### Load dependencies (all packages loaded in data-processing.R)
# source("script/data-processing.R")

################################################################################
# SECTION 1: Prepare Data for Heatmap
################################################################################

# Prepare data for heatmap
phase_heatmap_data <- climate_dengue_phase_df %>%
  mutate(
    climate_var_label = case_when(
      climate_var == "precipitation_era" ~ "Precipitation",
      climate_var == "temperature_era" ~ "Temperature",
      climate_var == "rel_humidity_era" ~ "Humidity",
      TRUE ~ climate_var
    ),
    climate_var_label = factor(climate_var_label,
                                levels = c("Precipitation", "Temperature", "Humidity"))
  ) %>%
  left_join(admin1_EN, by = "idadmin1") %>%
  arrange(idadmin1)

phase_heatmap_data <- phase_heatmap_data %>%
  mutate(shapeName = factor(shapeName, levels = levels(admin1_national_EN$shapeName)[-1]))

################################################################################
# SECTION 2: Figure 5A - Phase Lag Heatmap
################################################################################

fig_5A_supp <- phase_heatmap_data %>%
  ggplot(aes(x = climate_var_label, y = shapeName, fill = phase_lag_months)) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-6, 6),
    oob = scales::squish,
    name = "Phase lag\n(months)"
  )  +
  # Add region separators
  geom_hline(yintercept = as.numeric(admin1_national_EN$shapeName[c(11, 18, 20, 25, 31, 33)])-0.5,
             color = "black", size = 0.4) +
  theme_Publication(base_size = 10) +
  geom_text(aes(label = round(phase_lag_months,1)), color = "black", size = 3, fontface = "bold") +
  labs(
    x = NULL,
    y = NULL,
    tag = "A"
  ) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

fig_5A <- phase_heatmap_data %>% filter(climate_var != "rel_humidity_era") %>% 
  ggplot(aes(x = climate_var_label, y = shapeName, fill = phase_lag_months)) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-6, 6),
    oob = scales::squish,
    name = "Phase lag\n(months)"
  )  +
  # Add region separators
  geom_hline(yintercept = as.numeric(admin1_national_EN$shapeName[c(11, 18, 20, 25, 31, 33)])-0.5,
             color = "black", size = 0.4) +
  theme_Publication(base_size = 10) +
  geom_text(aes(label = round(phase_lag_months,1)), color = "black", size = 3, fontface = "bold") +
  labs(
    x = NULL,
    y = NULL,
    tag = "A"
  ) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

################################################################################
# SECTION 3: Figure 5B - Phase Coherence Heatmap
################################################################################

fig_5B_supp <- phase_heatmap_data %>%
  ggplot(aes(x = climate_var_label, y = shapeName, fill = phase_coherence*100)) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_viridis_c(
    option = "plasma",
    limits = c(0, 100),
    name = "Phase\ncoherence (%)"
  ) +
  geom_hline(yintercept = as.numeric(admin1_national_EN$shapeName[c(11, 18, 20, 25, 31, 33)])-0.5,
             color = "black", size = 0.4) +
  theme_Publication(base_size = 10) +
  geom_text(aes(label = round(phase_coherence*100,0)), color = "black", size = 3, fontface = "bold") +
  labs(
    x = NULL,
    y = NULL,
    tag = "B"
  ) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

fig_5B <- phase_heatmap_data %>% filter(climate_var != "rel_humidity_era") %>% 
  ggplot(aes(x = climate_var_label, y = shapeName, fill = phase_coherence*100)) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_viridis_c(
    option = "plasma",
    limits = c(0, 100),
    name = "Phase\ncoherence (%)"
  ) +
  geom_hline(yintercept = as.numeric(admin1_national_EN$shapeName[c(11, 18, 20, 25, 31, 33)])-0.5,
             color = "black", size = 0.4) +
  theme_Publication(base_size = 10) +
  geom_text(aes(label = round(phase_coherence*100,0)), color = "black", size = 3, fontface = "bold") +
  labs(
    x = NULL,
    y = NULL,
    tag = "B"
  ) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

################################################################################
# SECTION 4: Combined Figure and Save
################################################################################

# Combined figure (A + B side by side)
fig_5_combined <- fig_5A + fig_5B +
  plot_layout(widths = c(1, 1))

fig_5_supp_combined <- fig_5A_supp + fig_5B_supp +
  plot_layout(widths = c(1, 1))

# Save figures
ggsave("output/fig_5.jpg", fig_5_combined, height = 19, width = 17, unit = "cm", dpi = 600)
ggsave("output/fig_5_supp.jpg", fig_5_supp_combined, height = 19, width = 17, unit = "cm", dpi = 600)