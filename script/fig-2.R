################################################################################
# Figure 2: Wavelet Phase Analysis and Spatial Gradients
#
# This script generates Figure 2 for the manuscript, showing:
#   - 2A: Choropleth map of average phase lags (diverging scale)
#   - 2B: Choropleth map of peak outbreak month
#   - 2C: West-to-east gradient analysis by epidemic year
#
# Methods:
#   - Continuous wavelet transform (Morlet) at 12-month period
#   - July-June epidemic year for epidemic period definition
#
# Dependencies: Requires data-processing.R to be run first
################################################################################

#### Load dependencies (all packages loaded in data-processing.R)
# source("script/data-processing.R")

################################################################################
# SECTION 1: Supplementary Figure - Phase Lags with Error Bars
################################################################################

# Create detailed plot showing phase lags with confidence intervals for all provinces
# This is used as a supplementary figure
fig_phase_stats_supplementary <- province_phase_stats %>% 
  dplyr::select(idadmin1, region, admin1 = province, median_phase_lag, lower_ci, upper_ci) %>% 
  left_join(admin1_EN) %>% 
  mutate(admin1_name = factor(shapeName,
                              levels = as.character(admin1_EN$shapeName)[match(levels(province_phase_stats$province),
                                                                                admin1_indonesia$admin1)])) %>%
  ggplot(aes(x = admin1_name, y = median_phase_lag, col = region)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci, col = region), width = 0.3) +
  theme_minimal() +
  labs(
    x = NULL,
    y = "Phase lag from other provinces\n(in months)",
    col = NULL
  ) +
  scale_colour_colorblind() +
  coord_cartesian(ylim = c(-4.5, 4.5)) +
  scale_y_continuous(breaks = seq(-5, 5, 1)) +
  geom_hline(yintercept = 0, colour = "red", linetype = 2) +
  theme_Publication(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black")
  ) + 
  guides(color = guide_legend(nrow = 1))

ggsave("output/fig_phase_stats_supplementary.jpg", fig_phase_stats_supplementary, height = 10, width = 17, unit = "cm", dpi = 600)

################################################################################
# SECTION 2: Figure 2A - Average Phase Lag Map
################################################################################

# Create choropleth showing median phase lag across all provinces
# Uses diverging color scale: red = positive lag (late), blue = negative lag (early)
fig_2A <- ggplot() +
  geom_sf(data = other_shp, fill = "gray80", size=0.2) +  
  geom_sf(data = admin1_shp_with_lags, aes(fill = median_phase_lag), color = "black", size = 0.2) +
  scale_fill_gradient2(
    low = "#4662D7FF", 
    mid = "white", 
    high = "#CB2A04FF", 
    midpoint = 0,
    limits = c(-3, 3),
    name = "Average phase lag\n(in months)",
  ) +
  theme_Publication(base_size = 10) +
  labs(tag = "A", x = NULL, y = NULL) +
  geom_shadowtext(
    data = admin1_ir_centroids,
    aes(x = X, y = Y, label = admin1_no),
    size = 2.0,  
    color = "white",
    bg.color = 'grey10',
    fontface = "bold",  
    check_overlap = FALSE
  ) +
  coord_sf(xlim = c(95.01098, 141.0194),
           ylim = c(-11.00759, 5.906897))

################################################################################
# SECTION 3: Figure 2B - Peak Outbreak Month Map
################################################################################

# Create choropleth showing median peak month (derived from wavelet phase analysis)
# Peak month represents the typical time of maximum incidence within July-June epidemic year
fig_2B <- ggplot() +
  geom_sf(data = other_shp, fill = "gray80", size=0.2) +  
  geom_sf(data = admin1_shp_with_peak, aes(fill = median_peak), color = "black", size = 0.2) +
  theme_Publication(base_size = 10) +
  labs(tag = "B", x = NULL, y = NULL, fill="Peak month") +
  geom_shadowtext(
    data = admin1_ir_centroids,
    aes(x = X, y = Y, label = admin1_no),
    size = 2.0,  
    color = "white",
    bg.color = 'grey10',
    fontface = "bold",  
    check_overlap = FALSE
  ) +
  coord_sf(xlim = c(95.01098, 141.0194),
           ylim = c(-11.00759, 5.906897)) +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

################################################################################
# SECTION 4: Combine and Save Figure 2
################################################################################

# Combine panels A and B vertically
fig_2 <- (fig_2A / fig_2B) + plot_layout(heights=c(1,1))

# Save the combined figure
ggsave("output/fig_2.jpg", fig_2, height = 20, width = 17, unit = "cm", dpi = 600)
