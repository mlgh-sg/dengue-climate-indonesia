################################################################################
# Figure 4: Climate Summary and Climate-Dengue Relationships
#
# This script generates Figure 4 for the manuscript, showing:
#   - 4A: ONI/DMI time series (2009-2024)
#   - 4B: Precipitation heatmap by province (z-scored)
#   - 4C: Temperature heatmap by province (z-scored)
#   - 4D: Relative humidity heatmap by province (z-scored)
#
# Dependencies: Requires data-processing.R to be run first
################################################################################

#### Load dependencies (all packages loaded in data-processing.R)
# source("script/data-processing.R")

################################################################################
# SECTION 1: Figure 4A - ONI and DMI Time Series
################################################################################

# Create time series plot showing large-scale climate indices (2009-2024)
# ONI: Oceanic Niño Index (ENSO indicator)
# DMI: Dipole Mode Index (Indian Ocean Dipole indicator)
# Horizontal lines indicate different intensity thresholds
fig_oni_dmi_A <- oni_iod_2009_2024 %>% 
  filter(year < 2025) %>% 
  dplyr::select(date, oni, dmi) %>% 
  pivot_longer(-date, names_to = "var", values_to = "val") %>% 
  mutate(var = factor(var, levels = c("dmi", "oni"),
                      labels = c("DMI", "ONI"))) %>% 
  ggplot(aes(x = date, y = val, colour = var)) +
  theme_Publication(base_size = 10) +
  labs(x = NULL, y = "ONI/DMI", colour = NULL, tag = "A") +
  # scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_x_date(
    breaks = as.Date(paste0(2009:2024, "-07-01")),
    labels = 2009:2024,
    expand = c(0, 0),
    name = NULL
  ) +
  geom_vline(
    xintercept = as.numeric(as.Date(paste0(2009:2024, "-01-01"))),
    color = "gray",
    linewidth = 0.5, linetype = 2
  ) +
  scale_colour_colorblind() +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.25, linetype = 2) +
  geom_hline(yintercept = -0.5, colour = "#4662D7FF", linewidth = 0.25, linetype = 2) +
  geom_hline(yintercept = -1, colour = "#4662D7FF", linewidth = 0.25) +
  geom_hline(yintercept = -1.5, colour = "#4662D7FF", linewidth = 1) +
  geom_hline(yintercept = -2, colour = "#4662D7FF", linewidth = 1.5) +
  geom_hline(yintercept = 0.5, colour = "#CB2A04FF", linewidth = 0.25, linetype = 2) +
  geom_hline(yintercept = 1, colour = "#CB2A04FF", linewidth = 0.25) +
  geom_hline(yintercept = 1.5, colour = "#CB2A04FF", linewidth = 1) +
  geom_hline(yintercept = 2, colour = "#CB2A04FF", linewidth = 1.5) +
  geom_line(linewidth = 1.75) +
  theme(
    # axis.text.y = element_text(size = 7),
    # axis.text.x = element_text(size = 8),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid = element_blank(),
    # legend.position = "none",
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.spacing = unit(0, "lines")
  )

################################################################################
# SECTION 2: Figure 4B - Precipitation Heatmap
################################################################################

# Create heatmap showing precipitation patterns across provinces and time
# Values are z-scored within each province to show relative deviations
fig_climate_monthly_B <- climate_era_2009_2024_nolag %>% 
  left_join(admin1_national_EN) %>% 
  # filter(!(shapeName %in% c("Southwest Papua","South Papua","Central Papua","Highland Papua"))) %>% 
  mutate(
    date_start = floor_date(date, "month"),
    date_end = ceiling_date(date, "month") - days(1)
  ) %>%
  ggplot() +
  geom_rect(
    aes(
      xmin = date_start, 
      xmax = date_end + days(1),
      ymin = as.numeric(factor(shapeName)) - 0.5,
      ymax = as.numeric(factor(shapeName)) + 0.5,
      fill = precipitation_era_scaled
    )
  ) +
  geom_vline(
    xintercept = as.numeric(as.Date(paste0(2009:2024, "-01-01"))),
    color = "black",
    linewidth = 0.2, linetype = 2
  ) +
  scale_y_discrete(
    limits = rev(unique(admin1_national_EN$shapeName))
  ) +
  scale_fill_viridis_c(option="turbo",limits=c(-3, 3),oob=scales::squish) +
  theme_Publication(base_size = 10) +
  labs(x = NULL, y = "Province", tag = "B", title = "Precipitation") +
  scale_x_date(
    breaks = as.Date(paste0(2009:2024, "-07-01")),
    labels = 2009:2024,
    expand = c(0, 0),
    name = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 5.5),
    # axis.text.x = element_text(size = 8),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.spacing = unit(0, "lines")
  ) +
  # Add horizontal lines to separate geographic regions
  geom_hline(yintercept = as.numeric(admin1_national_EN$shapeName[c(11, 18, 20, 25, 31, 33, 35)]) + 0.5)

################################################################################
# SECTION 3: Figure 4C - Temperature Heatmap
################################################################################

# Create heatmap showing temperature patterns across provinces and time
# Values are z-scored within each province
fig_climate_monthly_C <- climate_era_2009_2024_nolag %>% 
  left_join(admin1_national_EN) %>%
  # filter(!(shapeName %in% c("Southwest Papua","South Papua","Central Papua","Highland Papua"))) %>% 
  mutate(
    date_start = floor_date(date, "month"),
    date_end = ceiling_date(date, "month") - days(1)
  ) %>%
  ggplot() +
  geom_rect(
    aes(
      xmin = date_start, 
      xmax = date_end + days(1),
      ymin = as.numeric(factor(shapeName)) - 0.5,
      ymax = as.numeric(factor(shapeName)) + 0.5,
      fill = temperature_era_scaled
    )
  ) +
  geom_vline(
    xintercept = as.numeric(as.Date(paste0(2009:2024, "-01-01"))),
    color = "black",
    linewidth = 0.2, linetype = 2
  ) +
  scale_y_discrete(
    limits = rev(unique(admin1_national_EN$shapeName))
  ) +
  scale_fill_viridis_c(option="turbo",limits=c(-3, 3),oob=scales::squish) +
  theme_Publication(base_size = 10) +
  labs(x = NULL, y = "Province", tag = "C", title = "Temperature") +
  scale_x_date(
    breaks = as.Date(paste0(2009:2024, "-07-01")),
    labels = 2009:2024,
    expand = c(0, 0),
    name = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 5.5),
    # axis.text.x = element_text(size = 8),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.spacing = unit(0, "lines")
  ) +
  # Add horizontal lines to separate geographic regions
  geom_hline(yintercept = as.numeric(admin1_national_EN$shapeName[c(11, 18, 20, 25, 31, 33, 35)]) + 0.5)

################################################################################
# SECTION 4: Figure 4D - Relative Humidity Heatmap
################################################################################

# Create heatmap showing relative humidity patterns across provinces and time
# Values are z-scored within each province
fig_climate_monthly_D <- climate_era_2009_2024_nolag %>% 
  left_join(admin1_national_EN) %>%
  # filter(!(shapeName %in% c("Southwest Papua","South Papua","Central Papua","Highland Papua"))) %>% 
  mutate(
    date_start = floor_date(date, "month"),
    date_end = ceiling_date(date, "month") - days(1)
  ) %>%
  ggplot() +
  geom_rect(
    aes(
      xmin = date_start, 
      xmax = date_end + days(1),
      ymin = as.numeric(factor(shapeName)) - 0.5,
      ymax = as.numeric(factor(shapeName)) + 0.5,
      fill = rel_humidity_era_scaled
    )
  ) +
  geom_vline(
    xintercept = as.numeric(as.Date(paste0(2009:2024, "-01-01"))),
    color = "black",
    linewidth = 0.2, linetype = 2
  ) +
  scale_y_discrete(
    limits = rev(unique(admin1_national_EN$shapeName))
  ) +
  scale_fill_viridis_c(option="turbo",limits=c(-3, 3),oob=scales::squish) +
  theme_Publication(base_size = 10) +
  labs(x = NULL, y = "Province", tag = "D", title = "Relative humidity",
       fill = "z-score\n(Climate variables)") +
  scale_x_date(
    breaks = as.Date(paste0(2009:2024, "-07-01")),
    labels = 2009:2024,
    expand = c(0, 0),
    name = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 5.5),
    # axis.text.x = element_text(size = 8),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    # legend.position = "none",
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.spacing = unit(0, "lines")
  ) +
  # Add horizontal lines to separate geographic regions
  geom_hline(yintercept = as.numeric(admin1_national_EN$shapeName[c(11, 18, 20, 25, 31, 33, 35)]) + 0.5)

################################################################################
# SECTION 5: Combine and Save Figure 4
################################################################################

# Combine all panels vertically: ONI/DMI time series + 3 climate heatmaps
fig_4 <- fig_oni_dmi_A / fig_climate_monthly_B / fig_climate_monthly_C / fig_climate_monthly_D

# Save combined figure in multiple formats
ggsave("output/fig_4.jpg", fig_4, height = 28, width = 17, unit = "cm", dpi = 600)
ggsave("output/fig_4.svg", fig_4, height = 28, width = 17, unit = "cm", dpi = 600)