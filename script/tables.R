################################################################################
# Main Tables Generation
#
# This script generates the main text tables:
#   - Table 1: Within-region cohesion (mean DTW distance)
#   - Table 2: Climate-dengue phase relationships by region
#   - Table 3: West-east gradient correlations by province subset
#
# Dependencies: Requires data-processing.R to be run first
################################################################################

#### Load dependencies (all packages loaded in data-processing.R)
# source("script/data-processing.R")

################################################################################
# TABLE 1: West-East Gradient Correlations (Main Text)
################################################################################

# Combine yearly gradient results from all three analyses
table_1_gradient <- bind_rows(
  yearly_gradient %>%
    mutate(subset = "All provinces (n=34)"),
  yearly_gradient_west %>%
    mutate(subset = "Western + Kalimantan (n=22)"),
  yearly_gradient_west2 %>%
    mutate(subset = "Western only (n=17)")
) %>%
  mutate(
    cor_longitude = round(cor_longitude, 3),
    p_value = round(p_value, 4),
    subset = factor(subset, levels = c("All provinces (n=34)",
                                       "Western + Kalimantan (n=22)",
                                       "Western only (n=17)"))
  ) %>%
  dplyr::select(
    `Epidemic year` = Period,
    Subset = subset,
    `Spearman rho` = cor_longitude,
    `p-value` = p_value,
    Significance = significance
  )

# Pivot to wide format for cleaner presentation
table_1_gradient_wide <- table_1_gradient %>%
  mutate(value = paste0(`Spearman rho`, " (", Significance, ")")) %>%
  dplyr::select(`Epidemic year`, Subset, value) %>%
  pivot_wider(names_from = Subset, values_from = value)

# Summary statistics
gradient_summary <- bind_rows(
  yearly_gradient %>%
    mutate(subset = "All provinces") %>%
    summarise(
      subset = first(subset),
      n_significant = sum(p_value < 0.05),
      n_total = n(),
      mean_rho = mean(cor_longitude),
      median_rho = median(cor_longitude)
    ),
  yearly_gradient_west %>%
    mutate(subset = "Western + Kalimantan") %>%
    summarise(
      subset = first(subset),
      n_significant = sum(p_value < 0.05),
      n_total = n(),
      mean_rho = mean(cor_longitude),
      median_rho = median(cor_longitude)
    ),
  yearly_gradient_west2 %>%
    mutate(subset = "Western only") %>%
    summarise(
      subset = first(subset),
      n_significant = sum(p_value < 0.05),
      n_total = n(),
      mean_rho = mean(cor_longitude),
      median_rho = median(cor_longitude)
    )
) %>%
  mutate(
    pct_significant = round(n_significant / n_total * 100, 1),
    mean_rho = round(mean_rho, 3),
    median_rho = round(median_rho, 3)
  )

# Create gt table
table_1_gt <- table_1_gradient_wide %>%
  gt() %>%
  tab_header(
    title = "Table 1. West-to-east gradient in peak timing by province subset",
    subtitle = "Spearman correlation between longitude and peak month (July-June scale)"
  ) %>%
  tab_footnote(
    footnote = "Positive correlation indicates westward provinces peak earlier. Significance: *** p<0.01, ** p<0.05, * p<0.10, n.s. = not significant.",
    locations = cells_column_labels(columns = `All provinces (n=34)`)
  ) %>%
  cols_align(align = "center", columns = everything()) %>%
  cols_align(align = "left", columns = `Epidemic year`)

# Save tables
gtsave(table_1_gt, "output/table_1_gradient.html")
write_xlsx(table_1_gradient_wide, "output/table_1_gradient.xlsx")
write_xlsx(gradient_summary, "output/table_1_gradient_summary.xlsx")

print("Table 1 saved to output/")
print(table_1_gradient_wide)
print("\nGradient summary:")
print(gradient_summary)

################################################################################
# TABLE 2: Within-Region Cohesion (Main Text)
################################################################################

table_2_region_cohesion <- region_cohesion %>%
  arrange(mean_dtw_distance) %>%
  mutate(
    across(where(is.numeric), ~round(., 2)),
    rank = row_number()
  ) %>%
  dplyr::select(
    Region = region,
    `N provinces` = n_provinces,
    `Mean DTW distance` = mean_dtw_distance,
    `SD` = sd_dtw_distance,
    `CV` = cv_dtw_distance,
    Rank = rank
  )

# Create gt table
table_2_gt <- table_2_region_cohesion %>%
  gt() %>%
  tab_header(
    title = "Table 2. Within-region temporal pattern cohesion",
    subtitle = "Based on pairwise DTW distances between provinces"
  ) %>%
  tab_footnote(
    footnote = "Lower mean DTW distance indicates more similar temporal dynamics within the region. CV = coefficient of variation.",
    locations = cells_column_labels(columns = `Mean DTW distance`)
  ) %>%
  cols_align(align = "center", columns = everything()) %>%
  cols_align(align = "left", columns = Region)

# Save table
gtsave(table_2_gt, "output/table_2_region_cohesion.html")
write_xlsx(table_2_region_cohesion, "output/table_2_region_cohesion.xlsx")

print("Table 2 saved to output/")
print(table_2_region_cohesion)

################################################################################
# Summary for Manuscript
################################################################################

cat("\n")
cat("================================================================================\n")
cat("SUMMARY STATISTICS FOR MANUSCRIPT\n")
cat("================================================================================\n")

cat("\n--- REGION COHESION ---\n")
cat("Most cohesive region:",
    region_cohesion$region[which.min(region_cohesion$mean_dtw_distance)],
    "(mean DTW =", round(min(region_cohesion$mean_dtw_distance, na.rm = TRUE), 2), ")\n")
cat("Least cohesive region:",
    region_cohesion$region[which.max(region_cohesion$mean_dtw_distance)],
    "(mean DTW =", round(max(region_cohesion$mean_dtw_distance, na.rm = TRUE), 2), ")\n")

cat("\n--- PHASE RELATIONSHIPS ---\n")
phase_summary <- climate_dengue_phase_df %>%
  group_by(climate_var) %>%
  summarise(
    mean_lag = mean(phase_lag_months, na.rm = TRUE),
    sd_lag = sd(phase_lag_months, na.rm = TRUE),
    mean_coherence = mean(phase_coherence, na.rm = TRUE)
  )
print(phase_summary)

cat("\n================================================================================\n")
cat("All main tables saved to output/ directory\n")
cat("================================================================================\n")
