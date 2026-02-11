################################################################################
# Figure 6: Cumulative Relative Risk Forest Plots
#
# This script generates Figure 6 for the manuscript, showing:
#   - 6A: Precipitation cumulative RR at wavelet-derived lag
#   - 6B: Temperature cumulative RR at wavelet-derived lag
#
# Dependencies: Requires data-processing.R to be run first
################################################################################

#### Load dependencies (all packages loaded in data-processing.R)
# source("script/data-processing.R")

################################################################################
# SECTION 1: Forest Plot for Precipitation (90th Percentile)
################################################################################

# Create forest plot showing cumulative RR for provinces with high coherence
plot_pre_RR_pct90 <- phase_heatmap_data_pct90_over85_RR %>%
  filter(clim_var == "Precipitation") %>%
  left_join(admin1_national_EN) %>%
  filter(idadmin1 != 99) %>%
  ggplot(aes(x = cumulRR, y = shapeName,
             xmin = cumulRR_lo, xmax = cumulRR_hi, col=as.character(lag))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbarh(width = 0.3, linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~clim_var, scales = "free_x", ncol = 3) +
  # scale_x_log10(labels = scales::label_number()) +
  scale_colour_manual(breaks=as.character(0:4),values=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442")) +
  theme_Publication(base_size = 10) +
  labs(
    x = "Cumulative Relative Risk (95% CI)",
    y = NULL,
    # title = "Cumulative effect of high exposure (90th percentile)",
    color = "Effect"
  ) +
  theme(
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  )

################################################################################
# SECTION 2: Forest Plot for Temperature (90th Percentile)
################################################################################

# Create forest plot showing cumulative RR for provinces with high coherence
plot_tem_RR_pct90 <- phase_heatmap_data_pct90_over85_RR %>%
  filter(clim_var == "Temperature") %>%
  left_join(admin1_national_EN) %>%
  filter(idadmin1 != 99) %>%
  ggplot(aes(x = cumulRR, y = shapeName,
             xmin = cumulRR_lo, xmax = cumulRR_hi, col=as.character(lag))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbarh(width = 0.3, linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~clim_var, scales = "free_x", ncol = 3) +
  # scale_x_log10(labels = scales::label_number()) +
  scale_colour_manual(breaks=as.character(0:4),values=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442")) +
  theme_Publication(base_size = 10) +
  labs(
    x = "Cumulative Relative Risk (95% CI)",
    y = NULL,
    # title = "Cumulative effect of high exposure (90th percentile)",
    color = "Effect"
  ) +
  theme(
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  )

################################################################################
# SECTION 3: Create Shared Legend
################################################################################

# Create dummy data with all lags 0-4 to generate a complete legend
dummy_data <- tibble(
  lag = 0:4,
  cumulRR = c(1, 1, 1, 1, 1),
  cumulRR_lo = c(0.9, 0.9, 0.9, 0.9, 0.9),
  cumulRR_hi = c(1.1, 1.1, 1.1, 1.1, 1.1),
  clim_var = "Dummy",
  shapeName = "Dummy"
)

# Create the legend plot
legend_plot <- dummy_data %>%
  ggplot(aes(x = cumulRR, y = shapeName,
             xmin = cumulRR_lo, xmax = cumulRR_hi, col=as.character(lag))) +
  geom_errorbarh(width = 0.3, linewidth = 0.5) +
  geom_point(size = 2) +
  scale_colour_manual(
    breaks = as.character(0:4),
    values = c("#000000","#E69F00","#56B4E9","#009E73","#F0E442")
  ) +
  theme_Publication(base_size = 10) +
  labs(color = "Lag") +
  theme(
    legend.position = "bottom"
  )

# Extract just the legend
legend <- cowplot::get_legend(legend_plot)

################################################################################
# SECTION 4: Combine Panels and Save
################################################################################

# Remove legends from individual plots
plot_pre_RR_pct90 <- plot_pre_RR_pct90 + theme(legend.position = "none")
plot_tem_RR_pct90 <- plot_tem_RR_pct90 + theme(legend.position = "none")

# Combine panels: precipitation / temperature / shared legend
fig_6 <- (plot_pre_RR_pct90 / plot_tem_RR_pct90) / legend +
  plot_layout(heights = c(1, 0.4, 0.1))  # Adjust the 0.1 to control legend size

# Save combined figure
ggsave("output/fig_6.jpg", fig_6, height = 15, width = 17, unit = "cm", dpi = 600)