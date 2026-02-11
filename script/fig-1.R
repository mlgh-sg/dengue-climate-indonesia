################################################################################
# Figure 1: Spatiotemporal Patterns of Dengue Incidence
#
# This script generates Figure 1 for the manuscript, showing:
#   - 1A: Choropleth map of average provincial incidence (2010-2024)
#   - 1B: Time series of regional and national incidence (log scale)
#   - 1C: Heatmap of provincial incidence over time (z-scored)
#
# Additional outputs:
#   - Monthly average heatmap (supplementary)
#
# Dependencies: Requires data-processing.R to be run first
################################################################################

#### Load dependencies (all packages loaded in data-processing.R)
# source("script/data-processing.R")

################################################################################
# SECTION 0: Create Supplementary Regional Map
################################################################################

# Calculate centroids for region labels on map
region_centroids <- region_shp %>%
  st_centroid() %>%
  cbind(st_coordinates(.)) %>%
  mutate(region = factor(region_shp$region,levels=c("Sumatra","Java & Bali","Kalimantan","Nusa Tenggara",
                                                    "Sulawesi","Maluku","Papua"))) %>%
  st_drop_geometry()

# Set region order for consistent display
region_shp$region <- factor(region_shp$region,levels=c("Sumatra","Java & Bali","Kalimantan","Nusa Tenggara",
                                                       "Sulawesi","Maluku","Papua"))

# Create regional map with labeled regions (for supplementary materials)
fig_region <- region_shp %>%
  ggplot() +
  # Add neighboring countries as background
  geom_sf(data = other_shp, fill = "gray80", linewidth = 0.2) +
  # Color regions
  geom_sf(aes(fill = region), color = 'black', linewidth = 0.2,
          show.legend = c(fill = TRUE, color = FALSE)) +
  # Add region labels
  geom_shadowtext(
    data = region_centroids,
    aes(x = X, y = Y, label = region),
    size = 2.0,
    color = "white",
    bg.color = 'grey10',
    fontface = "bold",
    check_overlap = FALSE
  ) +
  theme_Publication(base_size = 10) +
  scale_fill_colorblind() +
  coord_sf(xlim = c(95.01098, 141.0194),
           ylim = c(-11.00759, 5.906897)) +
  labs(x=NULL,y=NULL,fill=NULL)

# Save regional map
ggsave("output/fig_region.jpg", fig_region, height = 8, width = 17, unit = "cm", dpi = 600)

################################################################################
# SECTION 1: Figure 1A - Average Provincial Incidence Map
################################################################################

# Adjust centroid position for province 31 (overlapping label issue)
admin1_ir_centroids <- admin1_ir_centroids %>%
  mutate(Y=ifelse(idadmin1==31,-5.8,Y))

# Create choropleth map of average incidence rates (2010-2024)
fig_1A_plot <- admin1_ir_shp %>% 
  ggplot() +
  geom_sf(data = other_shp, fill = "gray80", linewidth = 0.2) +  
  geom_sf(aes(fill = incidence_rate), color = 'black', linewidth = 0.2,
          show.legend = c(fill = TRUE, color = FALSE)) +
  geom_shadowtext(
    data = admin1_ir_centroids,
    aes(x = X, y = Y, label = admin1_no),
    size = 2.0,  
    color = "white",
    bg.color = 'grey10',
    fontface = "bold",  
    check_overlap = FALSE
  ) +
  theme_Publication(base_size = 10) +
  scale_fill_viridis_c(
    option = "inferno",
    na.value = "gray80",
    breaks = c(0, 50, 100, 150, 200),
    limits = c(0, 200)
  ) +
  guides(fill = guide_colorbar(
    barheight = 0.5,
    label = "none"
  )) +
  coord_sf(xlim = c(95.01098, 141.0194),
           ylim = c(-11.00759, 5.906897)) +
  labs(fill = "Average dengue IR\n(2010-2024)", x = NULL, y = NULL, tag = "A")

# Extract legend for custom positioning
fig_1A_legend <- ggpubr::get_legend(fig_1A_plot)
fig_1A_plot2 <- fig_1A_plot + theme(legend.position = "none")

# Combine map and legend with custom positioning
fig_1A <- ggdraw() +
  draw_plot(fig_1A_plot2, 0, 0, 1, 1) +
  draw_plot(fig_1A_legend, 0, 0.075, 0.35, 0.25)

################################################################################
# SECTION 2: Figure 1B - Regional and National Time Series
################################################################################

# Create time series plot showing regional patterns and national aggregate
fig_1B <- dengue_data_regnat %>% 
  ggplot() +
  # Add vertical lines to mark year boundaries
  geom_vline(
    xintercept = as.numeric(as.Date(paste0(2010:2024, "-01-01"))),
    linetype = "dashed",
    color = "gray70",
    size = 0.3
  ) +
  # Regional time series (dashed lines)
  geom_line(aes(x = date, y = incidence_rate, col = region),
            size = 0.4,
            linetype = 5) +
  # National aggregate (thick solid line)
  geom_line(data=dengue_data_admin0,
            aes(x = date, y = incidence_rate),
            color = "#CC79A7",
            size = 1.5) +
  scale_y_continuous(
    # trans = "log10",
    labels = function(x) sprintf("%g", x),
    name = "Monthly dengue incidence rate per 100,000"
  ) +
  scale_x_date(
    breaks = as.Date(paste0(2010:2024, "-07-01")),
    labels = 2010:2024,
    expand = c(0.01, 0.01),
    name = NULL
  ) +
  scale_colour_colorblind() +
  theme_Publication(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.y = element_text(margin = margin(r = 0.1)),
    axis.ticks.x = element_blank()
  ) +
  # Customize legend to match line styles
  guides(color = guide_legend(override.aes = list(linewidth = c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,1.5),
                                                  linetype = c(5,5,5,5,5,5,5,1)))) +
  labs(tag = "B",col=NULL)

################################################################################
# SECTION 3: Figure 1C - Provincial Incidence Heatmap
################################################################################

# Create heatmap showing provincial incidence over time (z-scored within province)
fig_1C <- dengue_data_admin1 %>% 
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
      fill = incidence_rate_scaled
    )
  ) +
  geom_vline(
    xintercept = as.numeric(as.Date(paste0(2010:2024, "-01-01"))),
    color = "black",
    linewidth = 0.2, linetype = 2
  ) +
  scale_y_discrete(
    limits = rev(unique(admin1_national_EN$shapeName))
  ) +
  scale_fill_viridis_c(na.value = "#CCCCCC",option="turbo",limits=c(-3, 3),oob=scales::squish) +
  theme_Publication(base_size = 10) +
  labs(x = NULL, y = "Province", tag = "C", fill = "z-score of\nlog(Incidence Rate)") +
  scale_x_date(
    breaks = as.Date(paste0(2010:2024, "-07-01")),
    labels = 2010:2024,
    expand = c(0, 0),
    name = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    # legend.position = "none",
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    panel.spacing = unit(0, "lines")
  ) +
  # Add horizontal lines to separate geographic regions
  geom_hline(yintercept = as.numeric(admin1_national_EN$shapeName[c(11, 18, 20, 25, 31, 33, 35)]) + 0.5)

################################################################################
# SECTION 4: Combine All Panels and Save
################################################################################

# Combine panels A, B, C vertically
fig_1 <- (free(fig_1A_plot) / fig_1B / fig_1C) +
  plot_layout(heights = c(1, 1, 1))

# Save combined figure in multiple formats
ggsave("output/fig_1.jpg", fig_1, height = 28, width = 17, unit = "cm", dpi = 600)
ggsave("output/fig_1.svg", fig_1, height = 28, width = 17, unit = "cm", dpi = 600)