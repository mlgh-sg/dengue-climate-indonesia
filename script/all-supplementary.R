################################################################################
#### Supplementary Tables and Figures
####
#### This script generates all supplementary materials for the manuscript.
#### It requires data-processing.R to be run first (sourced below).
####
#### Outputs:
#### - Tables S1-S11 (as data frames and gt objects)
#### - Figures S1-S5 (as ggplot objects)
####
#### For Quarto rendering: All objects are created in the global environment
#### and can be directly referenced in the supplementary-materials.qmd document.
################################################################################

# Load dependencies (all packages loaded in data-processing.R)
source("script/data-processing.R")

cat("\n================================================================================\n")
cat("GENERATING SUPPLEMENTARY MATERIALS\n")
cat("================================================================================\n\n")

#### Tables
#### S1: Province names and their indexes (for reference to all maps in the main figure).
table_s1 <- admin1_EN %>% 
  mutate(index=1:34) %>% 
  dplyr::select(index,province=shapeName,idadmin1)

#### S2: Incidence by July-June annual period and their respective z-score to find major outbreak.
# Add epidemic year to dengue data
# Epidemic year "2015-16" = July 2015 to June 2016
dengue_with_epi_year <- dengue_data_admin1 %>%
  filter(idadmin1 != 99) %>%  # exclude national aggregate
  
  mutate(
    month = month(date),
    year = year(date),
    # Epidemic year: July starts new year, so Jul-Dec belong to current year label
    # Jan-Jun belong to previous year's epidemic year
    epi_year_start = if_else(month >= 7, year, year - 1),
    epi_year = paste0(epi_year_start, "-", substr(epi_year_start + 1, 3, 4))
  )

# Calculate annual national incidence by epidemic year
annual_national_epi <- dengue_with_epi_year %>%
  group_by(epi_year, epi_year_start) %>%
  summarise(
    total_cases = sum(cases, na.rm = TRUE),
    total_pop = sum(pop, na.rm = TRUE) / 12,  # average population
    incidence = total_cases / total_pop * 100000,
    n_months = n_distinct(date),  # should be 12 for complete years
    .groups = "drop"
  ) %>%
  # Keep only complete epidemic years (2010-11 to 2023-24)
  filter(epi_year_start >= 2010 & epi_year_start <= 2023)

# Calculate mean and SD for outbreak threshold
mean_incidence <- mean(annual_national_epi$incidence)
sd_incidence <- sd(annual_national_epi$incidence)

# Define outbreak years (> mean + 1 SD)
annual_national_epi <- annual_national_epi %>%
  mutate(
    z_score = (incidence - mean_incidence) / sd_incidence,
    outbreak_year = z_score > 1,
    year_type = case_when(
      z_score > 1 ~ "Major outbreak",
      z_score > 0.5 ~ "Above average",
      z_score < -0.5 ~ "Below average",
      TRUE ~ "Average"
    )
  )

cat("Epidemic year definition: July Year1 to June Year2\n")
cat("Example: '2015-16' = July 2015 - June 2016\n\n")

cat("Annual dengue incidence by epidemic year:\n")
table_s2 <- annual_national_epi %>%
  dplyr::select(epi_year, total_cases, incidence, z_score, year_type) %>%
  mutate(across(where(is.numeric), ~round(., 2)))

cat(sprintf("\nMean incidence: %.1f per 100,000\n", mean_incidence))
cat(sprintf("SD incidence: %.1f per 100,000\n", sd_incidence))
cat(sprintf("Outbreak threshold (mean + 1 SD): %.1f per 100,000\n", mean_incidence + sd_incidence))

#### S3: Silhouette statistics for optimal number of clusters (all provinces). Highlighted in yellow are the number of clusters selected for each scenario, accounting for the silhouette statistics and the balance in cluster sizes.
# highlight yellow no of `No. of clusters` = 12
table_s3 <- res_cvi %>%
  mutate(k = as.numeric(str_remove(k, "^k_"))) %>% 
  dplyr::select(`No. of clusters`=k,Silhouette=Sil)

#### S4: Silhouette statistics for optimal number of clusters (western provinces with Kalimantan). Highlighted in yellow are the number of clusters selected for each scenario, accounting for the silhouette statistics and the balance in cluster sizes.
# highlight yellow no of `No. of clusters` = 10
table_s4 <- res_cvi_west %>%
  mutate(k = as.numeric(str_remove(k, "^k_"))) %>% 
  dplyr::select(`No. of clusters`=k,Silhouette=Sil)

#### S5: Silhouette statistics for optimal number of clusters (western provinces only). Highlighted in yellow are the number of clusters selected for each scenario, accounting for the silhouette statistics and the balance in cluster sizes.
# highlight yellow no of `No. of clusters` = 4
table_s5 <- res_cvi_west2 %>%
  mutate(k = as.numeric(str_remove(k, "^k_"))) %>% 
  dplyr::select(`No. of clusters`=k,Silhouette=Sil)

#### S6: Province dissimilarity within regions. Mean DTW distance from each province to other provinces in the same region, with z-scores calculated within each region. Provinces with |z| > 1.645 are flagged as temporal outliers whose dynamics deviate substantially from their regional pattern.
table_s6_outliers <- region_outliers %>%
  arrange(region, desc(z_score)) %>%
  mutate(
    across(where(is.numeric) & !matches("idadmin1"), ~round(., 2)),
    outlier_flag = ifelse(is_outlier, "Yes", "No")
  ) %>%
  dplyr::select(
    Region = region,
    Province = admin1,
    `Mean DTW dist to region` = mean_dist_to_region,
    `Z-score` = z_score,
    `Outlier (|z| > 2)` = outlier_flag
  )

# Create gt table
table_s6 <- table_s6_outliers %>%
  gt() %>%
  tab_header(
    title = "Table S6. Province dissimilarity within regions",
    subtitle = " Mean DTW distance from each province to other provinces in the same region, with z-scores calculated within each region. Provinces with |z| > 1.645 are flagged as temporal outliers whose dynamics deviate substantially from their regional pattern."
  ) %>%
  tab_footnote(
    footnote = "Z-score calculated within each region. Provinces with |z| > 1.65 flagged as outliers.",
    locations = cells_column_labels(columns = `Z-score`)
  ) %>%
  cols_align(align = "center", columns = everything()) %>%
  cols_align(align = "left", columns = c(Region, Province))

#### S7: Regional Spearman correlations between large-scale climate indices and incidence.
# Codes to generate data
################################################################################
# 1. Calculate Incidence by Epidemic Year (July-June)
################################################################################

cat("\n================================================================================\n")
cat("SECTION 1: EPIDEMIC YEAR DEFINITION (JULY-JUNE)\n")
cat("================================================================================\n\n")

# Add epidemic year to dengue data
# Epidemic year "2015-16" = July 2015 to June 2016
dengue_with_epi_year <- dengue_data_admin1 %>%
  filter(idadmin1 != 99) %>%  # exclude national aggregate
  
  mutate(
    month = month(date),
    year = year(date),
    # Epidemic year: July starts new year, so Jul-Dec belong to current year label
    # Jan-Jun belong to previous year's epidemic year
    epi_year_start = if_else(month >= 7, year, year - 1),
    epi_year = paste0(epi_year_start, "-", substr(epi_year_start + 1, 3, 4))
  )

# Calculate annual national incidence by epidemic year
annual_national_epi <- dengue_with_epi_year %>%
  group_by(epi_year, epi_year_start) %>%
  summarise(
    total_cases = sum(cases, na.rm = TRUE),
    total_pop = sum(pop, na.rm = TRUE) / 12,  # average population
    incidence = total_cases / total_pop * 100000,
    n_months = n_distinct(date),  # should be 12 for complete years
    .groups = "drop"
  ) %>%
  # Keep only complete epidemic years (2010-11 to 2023-24)
  filter(epi_year_start >= 2010 & epi_year_start <= 2023)

# Calculate mean and SD for outbreak threshold
mean_incidence <- mean(annual_national_epi$incidence)
sd_incidence <- sd(annual_national_epi$incidence)

# Define outbreak years (> mean + 1 SD)
annual_national_epi <- annual_national_epi %>%
  mutate(
    z_score = (incidence - mean_incidence) / sd_incidence,
    outbreak_year = z_score > 1,
    year_type = case_when(
      z_score > 1 ~ "Major outbreak",
      z_score > 0.5 ~ "Above average",
      z_score < -0.5 ~ "Below average",
      TRUE ~ "Average"
    )
  )

cat("Epidemic year definition: July Year1 to June Year2\n")
cat("Example: '2015-16' = July 2015 - June 2016\n\n")

cat("Annual dengue incidence by epidemic year:\n")
print(annual_national_epi %>%
        dplyr::select(epi_year, total_cases, incidence, z_score, year_type) %>%
        mutate(across(where(is.numeric), ~round(., 2))))

cat(sprintf("\nMean incidence: %.1f per 100,000\n", mean_incidence))
cat(sprintf("SD incidence: %.1f per 100,000\n", sd_incidence))
cat(sprintf("Outbreak threshold (mean + 1 SD): %.1f per 100,000\n", mean_incidence + sd_incidence))

################################################################################
# 2. Calculate ONI and DMI by Epidemic Year (July-June)
################################################################################

cat("\n================================================================================\n")
cat("SECTION 2: CLIMATE INDICES BY EPIDEMIC YEAR\n")
cat("================================================================================\n\n")

# Add epidemic year to ONI data
oni_with_epi_year <- oni_2009_2024 %>%
  mutate(
    month = month(date),
    year = year(date),
    epi_year_start = if_else(month >= 7, year, year - 1),
    epi_year = paste0(epi_year_start, "-", substr(epi_year_start + 1, 3, 4))
  )

# Calculate epidemic-year average ONI
annual_oni_epi <- oni_with_epi_year %>%
  filter(epi_year_start >= 2010 & epi_year_start <= 2023) %>%
  group_by(epi_year, epi_year_start) %>%
  summarise(
    mean_oni = mean(oni, na.rm = TRUE),
    max_oni = max(oni, na.rm = TRUE),
    min_oni = min(oni, na.rm = TRUE),
    # Count months in El Nino (ONI > 0.5) or La Nina (ONI < -0.5)
    months_el_nino = sum(oni > 0.5),
    months_la_nina = sum(oni < -0.5),
    .groups = "drop"
  ) %>%
  mutate(
    enso_phase = case_when(
      months_el_nino >= 5 ~ "El Niño",
      months_la_nina >= 5 ~ "La Niña",
      TRUE ~ "Neutral"
    )
  )

# Add epidemic year to DMI data
dmi_with_epi_year <- iod_2009_2024 %>%
  mutate(
    month = month(date),
    year = year(date),
    epi_year_start = if_else(month >= 7, year, year - 1),
    epi_year = paste0(epi_year_start, "-", substr(epi_year_start + 1, 3, 4))
  )

# Calculate epidemic-year average DMI
annual_dmi_epi <- dmi_with_epi_year %>%
  filter(epi_year_start >= 2010 & epi_year_start <= 2023) %>%
  group_by(epi_year, epi_year_start) %>%
  summarise(
    mean_dmi = mean(dmi, na.rm = TRUE),
    max_dmi = max(dmi, na.rm = TRUE),
    min_dmi = min(dmi, na.rm = TRUE),
    # Count months in positive IOD (DMI > 0.4) or negative IOD (DMI < -0.4)
    months_positive_iod = sum(dmi > 0.4, na.rm = TRUE),
    months_negative_iod = sum(dmi < -0.4, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    iod_phase = case_when(
      months_positive_iod >= 3 ~ "Positive IOD",
      months_negative_iod >= 3 ~ "Negative IOD",
      TRUE ~ "Neutral IOD"
    )
  )

# Merge all data
indices_outbreak_epi <- annual_national_epi %>%
  left_join(annual_oni_epi, by = c("epi_year", "epi_year_start")) %>%
  left_join(annual_dmi_epi, by = c("epi_year", "epi_year_start"))

cat("Combined epidemic year data:\n")
print(indices_outbreak_epi %>%
        dplyr::select(epi_year, incidence, z_score, year_type,
                      mean_oni, enso_phase, mean_dmi, iod_phase) %>%
        mutate(across(where(is.numeric), ~round(., 2))))

################################################################################
# 3. ONI vs DMI Correlation Comparison
################################################################################

cat("\n================================================================================\n")
cat("SECTION 3: ONI vs DMI CORRELATION COMPARISON\n")
cat("================================================================================\n\n")

# Correlation: ONI vs incidence
cor_oni <- cor.test(indices_outbreak_epi$mean_oni, indices_outbreak_epi$incidence,
                    method = "spearman")
# Correlation: DMI vs incidence
cor_dmi <- cor.test(indices_outbreak_epi$mean_dmi, indices_outbreak_epi$incidence,
                    method = "spearman")
# Correlation: ONI vs DMI (redundancy check)
cor_oni_dmi <- cor.test(indices_outbreak_epi$mean_oni, indices_outbreak_epi$mean_dmi,
                        method = "spearman")

cat("Correlation with national incidence (epidemic year basis):\n")
cat(sprintf("  ONI: rho = %.3f, p = %.4f %s\n",
            cor_oni$estimate, cor_oni$p.value,
            ifelse(cor_oni$p.value < 0.05, "*", "")))
cat(sprintf("  DMI: rho = %.3f, p = %.4f %s\n",
            cor_dmi$estimate, cor_dmi$p.value,
            ifelse(cor_dmi$p.value < 0.05, "*", "")))
cat(sprintf("\nONI-DMI correlation: rho = %.3f, p = %.4f\n",
            cor_oni_dmi$estimate, cor_oni_dmi$p.value))
cat("  (Low correlation = indices capture different dynamics = potentially complementary)\n")
cat("  (High correlation = indices are redundant)\n")

################################################################################
# 4. Incidence by ENSO and IOD Phase
################################################################################

cat("\n================================================================================\n")
cat("SECTION 4: INCIDENCE BY CLIMATE PHASE\n")
cat("================================================================================\n\n")

# By ENSO phase
cat("Mean incidence by ENSO phase:\n")
enso_summary <- indices_outbreak_epi %>%
  group_by(enso_phase) %>%
  summarise(
    n_years = n(),
    mean_incidence = round(mean(incidence), 1),
    sd_incidence = round(sd(incidence), 1),
    outbreak_years = sum(outbreak_year),
    years = paste(epi_year, collapse = ", "),
    .groups = "drop"
  )
print(enso_summary)

# By IOD phase
cat("\nMean incidence by IOD phase:\n")
iod_summary <- indices_outbreak_epi %>%
  group_by(iod_phase) %>%
  summarise(
    n_years = n(),
    mean_incidence = round(mean(incidence), 1),
    sd_incidence = round(sd(incidence), 1),
    outbreak_years = sum(outbreak_year),
    years = paste(epi_year, collapse = ", "),
    .groups = "drop"
  )
print(iod_summary)

# Cross-tabulation: ENSO x IOD
cat("\nCross-tabulation of ENSO and IOD phases:\n")
cross_tab <- indices_outbreak_epi %>%
  group_by(enso_phase, iod_phase) %>%
  summarise(
    n_years = n(),
    mean_incidence = round(mean(incidence), 1),
    outbreak_years = sum(outbreak_year),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_incidence))
print(cross_tab)

################################################################################
# 5. Outbreak Years and Climate Status
################################################################################

cat("\n================================================================================\n")
cat("SECTION 5: OUTBREAK YEARS AND CLIMATE STATUS\n")
cat("================================================================================\n\n")

cat("Major outbreak years (z > 1) and their climate status:\n")
outbreak_details <- indices_outbreak_epi %>%
  filter(outbreak_year) %>%
  dplyr::select(epi_year, incidence, z_score, mean_oni, enso_phase, mean_dmi, iod_phase) %>%
  mutate(across(where(is.numeric), ~round(., 2)))
print(outbreak_details)

# Count outbreaks by phase
cat("\nOutbreak year counts by climate phase:\n")
cat(sprintf("  During El Niño: %d of %d outbreak years\n",
            sum(indices_outbreak_epi$outbreak_year & indices_outbreak_epi$enso_phase == "El Niño"),
            sum(indices_outbreak_epi$outbreak_year)))
cat(sprintf("  During Positive IOD: %d of %d outbreak years\n",
            sum(indices_outbreak_epi$outbreak_year & indices_outbreak_epi$iod_phase == "Positive IOD"),
            sum(indices_outbreak_epi$outbreak_year)))
cat(sprintf("  During both El Niño AND Positive IOD: %d of %d outbreak years\n",
            sum(indices_outbreak_epi$outbreak_year &
                  indices_outbreak_epi$enso_phase == "El Niño" &
                  indices_outbreak_epi$iod_phase == "Positive IOD"),
            sum(indices_outbreak_epi$outbreak_year)))

################################################################################
# 6. Model Comparison (Linear Regression)
################################################################################

cat("\n================================================================================\n")
cat("SECTION 6: MODEL COMPARISON\n")
cat("================================================================================\n\n")

# Simple linear models
model_oni <- lm(incidence ~ mean_oni, data = indices_outbreak_epi)
model_dmi <- lm(incidence ~ mean_dmi, data = indices_outbreak_epi)
model_both <- lm(incidence ~ mean_oni + mean_dmi, data = indices_outbreak_epi)

cat("Linear model comparison (R-squared):\n")
cat(sprintf("  ONI only:      R² = %.3f, Adj R² = %.3f\n",
            summary(model_oni)$r.squared, summary(model_oni)$adj.r.squared))
cat(sprintf("  DMI only:      R² = %.3f, Adj R² = %.3f\n",
            summary(model_dmi)$r.squared, summary(model_dmi)$adj.r.squared))
cat(sprintf("  ONI + DMI:     R² = %.3f, Adj R² = %.3f\n",
            summary(model_both)$r.squared, summary(model_both)$adj.r.squared))

# ANOVA comparison
cat("\nANOVA: Does adding DMI improve ONI-only model?\n")
anova_result <- anova(model_oni, model_both)
print(anova_result)

cat("\nCoefficients from combined model:\n")
print(round(summary(model_both)$coefficients, 4))

################################################################################
# 7. Lagged Analysis
################################################################################

cat("\n================================================================================\n")
cat("SECTION 7: LAGGED ANALYSIS\n")
cat("================================================================================\n\n")

# Create lagged versions
indices_lagged <- indices_outbreak_epi %>%
  arrange(epi_year_start) %>%
  mutate(
    oni_lag1 = lag(mean_oni, 1),  # Previous epidemic year's ONI
    dmi_lag1 = lag(mean_dmi, 1),  # Previous epidemic year's DMI
    enso_phase_lag1 = lag(enso_phase, 1),
    iod_phase_lag1 = lag(iod_phase, 1)
  )

# Lagged correlations
cor_oni_lag <- cor.test(indices_lagged$oni_lag1, indices_lagged$incidence,
                        method = "spearman", use = "complete.obs")
cor_dmi_lag <- cor.test(indices_lagged$dmi_lag1, indices_lagged$incidence,
                        method = "spearman", use = "complete.obs")

cat("Correlation comparison: Same year vs Lagged (previous epidemic year):\n\n")
cat("ONI:\n")
cat(sprintf("  Same epidemic year:     rho = %.3f, p = %.4f\n", cor_oni$estimate, cor_oni$p.value))
cat(sprintf("  Previous epidemic year: rho = %.3f, p = %.4f\n", cor_oni_lag$estimate, cor_oni_lag$p.value))

cat("\nDMI:\n")
cat(sprintf("  Same epidemic year:     rho = %.3f, p = %.4f\n", cor_dmi$estimate, cor_dmi$p.value))
cat(sprintf("  Previous epidemic year: rho = %.3f, p = %.4f\n", cor_dmi_lag$estimate, cor_dmi_lag$p.value))

# Calculate regional incidence by epidemic year
regional_incidence_epi <- dengue_with_epi_year %>%
  group_by(epi_year, epi_year_start, region) %>%
  summarise(
    total_cases = sum(cases, na.rm = TRUE),
    total_pop = sum(pop, na.rm = TRUE) / 12,
    incidence = total_cases / total_pop * 100000,
    .groups = "drop"
  ) %>%
  filter(epi_year_start >= 2010 & epi_year_start <= 2023) %>%
  left_join(annual_oni_epi %>% dplyr::select(epi_year, mean_oni), by = "epi_year") %>%
  left_join(annual_dmi_epi %>% dplyr::select(epi_year, mean_dmi), by = "epi_year")

# Correlations by region
cat("Correlation with incidence by region:\n\n")
table_s7 <- regional_incidence_epi %>%
  group_by(region) %>%
  summarise(
    n_years = n(),
    cor_oni = cor(mean_oni, incidence, method = "spearman", use = "complete.obs"),
    cor_dmi = cor(mean_dmi, incidence, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  ) %>%
  mutate(
    better_index = case_when(
      abs(cor_oni) > abs(cor_dmi) + 0.1 ~ "ONI",
      abs(cor_dmi) > abs(cor_oni) + 0.1 ~ "DMI",
      TRUE ~ "Similar"
    ),
    across(where(is.numeric), ~round(., 3))
  )

#### This section is for COVID-19 period sensitivity analysis
#### Processing Codes:
################################################################################
# SECTION 1: Create Pre-Pandemic Data Subsets (2010-2019 only)
################################################################################

cat("\n=== Creating pre-pandemic data subsets (2010-2019) ===\n")

# Subset dengue data to pre-pandemic period
dengue_data_admin1_prepandemic <- dengue_data_admin1 %>%
  filter(year(date) < 2020) %>%
  group_by(admin1) %>%
  mutate(
    incidence_rate_log = log(incidence_rate_NA),
    incidence_rate_scaled = as.numeric(scale(incidence_rate_log))
  ) %>%
  ungroup()

cat("Dengue data subset created: ",
    nrow(dengue_data_admin1_prepandemic), " rows (vs ",
    nrow(dengue_data_admin1), " in full dataset)\n")

# Subset climate data to pre-pandemic period
climate_2010_2019 <- climate_era_2009_2024_nolag %>%
  filter(year(date) >= 2010, year(date) < 2020) %>%
  group_by(admin1) %>%
  mutate(
    precipitation_era_scaled = as.numeric(scale(precipitation_era)),
    temperature_era_scaled = as.numeric(scale(temperature_era)),
    rel_humidity_era_scaled = as.numeric(scale(rel_humidity_era))
  ) %>%
  ungroup()

cat("Climate data subset created: ",
    nrow(climate_2010_2019), " rows\n")

################################################################################
# SECTION 2: Re-run Wavelet Analysis on Pre-Pandemic Data
################################################################################

cat("\n=== Running wavelet analysis on pre-pandemic data ===\n")

# Prepare province-level dataset for wavelet analysis
province_data_prepandemic <- dengue_data_admin1_prepandemic %>%
  filter(idadmin1 != 99)

# Prepare data for wavelet analysis
province_ts_prepandemic <- prepare_for_wavelet(
  province_data_prepandemic,
  "admin1",
  "incidence_rate_scaled"
)

# Run wavelet analysis
province_wavelet_prepandemic <- conduct_wavelet_analysis(
  province_ts_prepandemic$ts_data,
  province_ts_prepandemic$valid_groups
)

cat("Wavelet analysis completed for",
    length(province_wavelet_prepandemic$phase_angles), "provinces\n")

################################################################################
# SECTION 3: Calculate Phase Lags (Pre-Pandemic)
################################################################################

cat("\n=== Calculating phase lags (pre-pandemic) ===\n")

# Calculate phase differences
province_phase_months_prepandemic <- create_phase_matrix(
  province_wavelet_prepandemic$phase_angles,
  province_ts_prepandemic$valid_groups
)

# Calculate average phase lags
province_phase_lags_prepandemic <- calculate_phase_lags(
  province_phase_months_prepandemic
)

# Calculate phase statistics
province_phase_stats_prepandemic <- calculate_phase_stats(
  province_phase_months_prepandemic
) %>%
  left_join(
    admin1_indonesia %>%
      dplyr::select(region, idadmin1, province = admin1),
    by = "province"
  ) %>%
  mutate(
    region = to_title_case(region),
    region = ifelse(region == "Java Bali", "Java & Bali", region),
    region = factor(region, levels = c("Sumatra", "Java & Bali", "Kalimantan",
                                       "Nusa Tenggara", "Sulawesi", "Maluku", "Papua"))
  )

cat("Phase statistics calculated for",
    nrow(province_phase_stats_prepandemic), "provinces\n")

################################################################################
# SECTION 4: Re-run Climate-Dengue Phase Analysis (Pre-Pandemic)
################################################################################

cat("\n=== Running climate-dengue phase analysis (pre-pandemic) ===\n")

# Prepare climate data for wavelet
climate_for_wavelet_prepandemic <- climate_2010_2019 %>%
  filter(idadmin1 != 99) %>%
  arrange(idadmin1, date)

# Store climate phase angles
climate_phase_angles_prepandemic <- list()

climate_vars_wavelet <- c("precipitation_era_scaled",
                          "temperature_era_scaled",
                          "rel_humidity_era_scaled")

for (clim_var in climate_vars_wavelet) {
  
  # Prepare data for wavelet
  climate_ts <- prepare_for_wavelet(
    climate_for_wavelet_prepandemic,
    "admin1",
    clim_var
  )
  
  # Run wavelet analysis
  climate_wavelet <- conduct_wavelet_analysis(
    climate_ts$ts_data,
    climate_ts$valid_groups
  )
  
  climate_phase_angles_prepandemic[[clim_var]] <- climate_wavelet$phase_angles
  
  cat("  Wavelet completed for:", clim_var, "\n")
}

# Calculate phase lag and coherence for each province and climate variable
calculate_phase_lag <- function(phase_climate, phase_dengue) {
  diff <- phase_climate - phase_dengue
  diff <- ((diff + pi) %% (2 * pi)) - pi
  mean_diff <- mean(circular(diff, type = "angles", units = "radians"))
  lag_months <- as.numeric(mean_diff) * 12 / (2 * pi)
  return(lag_months)
}

calculate_phase_coherence <- function(phase_climate, phase_dengue) {
  diff <- phase_climate - phase_dengue
  diff <- ((diff + pi) %% (2 * pi)) - pi
  coherence <- 1 - (sd(diff, na.rm = TRUE) / pi)
  return(coherence)
}

# Calculate for all provinces and climate variables
climate_dengue_phase_results_prepandemic <- list()

for (prov_idx in seq_along(province_ts_prepandemic$valid_groups)) {
  
  province_name <- province_ts_prepandemic$valid_groups[prov_idx]
  
  for (clim_var in climate_vars_wavelet) {
    
    phase_climate <- climate_phase_angles_prepandemic[[clim_var]][[prov_idx]]
    phase_dengue <- province_wavelet_prepandemic$phase_angles[[prov_idx]]
    
    # Calculate phase lag and coherence
    lag <- calculate_phase_lag(phase_climate, phase_dengue)
    coherence <- calculate_phase_coherence(phase_climate, phase_dengue)
    
    climate_dengue_phase_results_prepandemic[[length(climate_dengue_phase_results_prepandemic) + 1]] <-
      tibble(
        province = province_name,
        climate_var = gsub("_era_scaled", "", clim_var),
        phase_lag_months = lag,
        phase_coherence = coherence
      )
  }
}

climate_dengue_phase_df_prepandemic <- bind_rows(climate_dengue_phase_results_prepandemic) %>%
  left_join(
    admin1_indonesia %>%
      dplyr::select(region, idadmin1, province = admin1),
    by = c("province" = "province")
  ) %>%
  mutate(
    region = to_title_case(region),
    region = ifelse(region == "Java Bali", "Java & Bali", region),
    climate_var = case_when(
      climate_var == "precipitation" ~ "Precipitation",
      climate_var == "temperature" ~ "Temperature",
      climate_var == "rel_humidity" ~ "Rel. Humidity"
    )
  )

cat("Phase analysis completed for",
    length(unique(climate_dengue_phase_df_prepandemic$province)),
    "provinces x",
    length(unique(climate_dengue_phase_df_prepandemic$climate_var)),
    "climate variables\n")

################################################################################
# SECTION 5: Identify Qualifying Provinces for DLNM (Pre-Pandemic)
################################################################################

cat("\n=== Identifying qualifying provinces (pre-pandemic) ===\n")

# Apply same criteria: coherence >= 0.85 and positive/zero phase lag
qualifying_provinces_prepandemic <- climate_dengue_phase_df_prepandemic %>%
  filter(phase_coherence >= 0.85, phase_lag_months >= 0) %>%
  group_by(climate_var) %>%
  summarise(
    n_provinces = n(),
    provinces = list(province)
  )

cat("Qualifying provinces by climate variable (coherence >= 0.85, lag >= 0):\n")
print(qualifying_provinces_prepandemic)

################################################################################
# SECTION 6: Re-run DLNM for Qualifying Provinces (Pre-Pandemic)
################################################################################

cat("\n=== Running DLNM for qualifying provinces (pre-pandemic) ===\n")
cat("NOTE: This may take several minutes...\n")

# Prepare DLNM data (pre-pandemic)
# Need to match the structure of dengue_data_dlnm from data-processing.R
# which includes: region, idadmin1, admin1, date, pop, cases, incidence,
# precipitation, temperature, rel_humidity, oni, dmi
dengue_data_dlnm_prepandemic <- dengue_data_admin1_prepandemic %>%
  # filter(idadmin1 != 99) %>%
  dplyr::select(region, idadmin1, admin1, date, pop, cases, incidence = incidence_rate) %>%
  left_join(
    prec_admin1_era_2009_2024 %>%
      filter(year(date) >= 2010, year(date) < 2020) %>%
      dplyr::select(region, idadmin1, admin1, date, precipitation)
  ) %>%
  left_join(
    temp_admin1_era_2009_2024 %>%
      filter(year(date) >= 2010, year(date) < 2020) %>%
      dplyr::select(region, idadmin1, admin1, date, temperature)
  ) %>%
  left_join(
    humi_admin1_era_2009_2024 %>%
      filter(year(date) >= 2010, year(date) < 2020) %>%
      dplyr::select(region, idadmin1, admin1, date, rel_humidity)
  ) %>%
  left_join(
    oni_2009_2024 %>%
      filter(year(date) >= 2010, year(date) < 2020) %>%
      dplyr::select(date, oni)
  ) %>%
  left_join(
    iod_2009_2024 %>%
      filter(year(date) >= 2010, year(date) < 2020) %>%
      dplyr::select(date, dmi)
  ) %>%
  arrange(idadmin1, date) %>%
  mutate(year = year(date), month = month(date))

# Get precipitation qualifying provinces
precip_qualifying_prepandemic <- climate_dengue_phase_df_prepandemic %>%
  filter(climate_var == "Precipitation",
         phase_coherence >= 0.85,
         round(phase_lag_months) >= 0) %>%
  pull(province)

# Get temperature qualifying provinces
temp_qualifying_prepandemic <- climate_dengue_phase_df_prepandemic %>%
  filter(climate_var == "Temperature",
         phase_coherence >= 0.85,
         round(phase_lag_months) >= 0) %>%
  pull(province)

# Run DLNM for all qualifying provinces
# Note: dlnm_run() runs all 5 climate variables at once (ONI, DMI, precip, temp, humidity)
# So we just need to run it once per province that appears in either precip or temp qualifying lists

all_qualifying_provinces <- unique(c(precip_qualifying_prepandemic, temp_qualifying_prepandemic))

dlnm_prepandemic <- list()

if (length(all_qualifying_provinces) > 0) {
  
  cat("\nRunning DLNM for", length(all_qualifying_provinces),
      "qualifying provinces...\n")
  cat("(This will take several minutes...)\n\n")
  
  for (prov in all_qualifying_provinces) {
    
    cat("  Processing:", prov, "\n")
    
    prov_data <- dengue_data_dlnm_prepandemic %>%
      filter(admin1 == prov)
    
    # Run DLNM for all climate variables
    dlnm_result <- dlnm_run(
      data_input = prov_data,
      lag_value_long = 6,   # for ONI/DMI
      lag_value_short = 4   # for local climate
    )
    
    dlnm_prepandemic[[prov]] <- dlnm_result
  }
}

cat("\nDLNM analysis completed!\n")

################################################################################
# SECTION 7: Compare Full vs Pre-Pandemic Results
################################################################################

cat("\n=== Creating comparison tables ===\n")

# Table 1: Phase Lags Comparison
phase_lags_comparison <- province_phase_stats %>%
  dplyr::select(province, region,
                phase_lag_full = median_phase_lag,
                lower_ci_full = lower_ci,
                upper_ci_full = upper_ci) %>%
  left_join(
    province_phase_stats_prepandemic %>%
      dplyr::select(province,
                    phase_lag_prepandemic = median_phase_lag,
                    lower_ci_prepandemic = lower_ci,
                    upper_ci_prepandemic = upper_ci),
    by = "province"
  ) %>%
  mutate(
    phase_lag_diff = phase_lag_full - phase_lag_prepandemic,
    qualitative_agreement = ifelse(
      sign(phase_lag_full) == sign(phase_lag_prepandemic),
      "Same direction",
      "Different direction"
    )
  ) %>%
  arrange(region, province)

cat("\nPhase lag comparison table created:", nrow(phase_lags_comparison), "provinces\n")

# Table 2: Phase Coherence Comparison (for precipitation and temperature)
phase_coherence_comparison <- climate_dengue_phase_df %>%
  filter(climate_var %in% c("precipitation_era", "temperature_era")) %>%
  mutate(climate_var = case_when(
    climate_var == "precipitation_era" ~ "Precipitation",
    climate_var == "temperature_era" ~ "Temperature",
    .default = "Rel. Humidity"
  )) %>% 
  filter(!is.na(idadmin1), idadmin1 != 99) %>%  # Exclude national level
  dplyr::select(admin1, region, climate_var,
                phase_lag_full = phase_lag_months,
                coherence_full = phase_coherence) %>%
  left_join(
    climate_dengue_phase_df_prepandemic %>%
      filter(climate_var %in% c("Precipitation", "Temperature")) %>%
      dplyr::select(admin1=province, climate_var,
                    phase_lag_prepandemic = phase_lag_months,
                    coherence_prepandemic = phase_coherence),
    by = c("admin1", "climate_var")
  ) %>%
  mutate(
    coherence_diff = coherence_full - coherence_prepandemic,
    lag_diff = phase_lag_full - phase_lag_prepandemic,
    qualified_full = (coherence_full >= 0.85 & phase_lag_full >= 0),
    qualified_prepandemic = (coherence_prepandemic >= 0.85 & phase_lag_prepandemic >= 0),
    qualification_agreement = ifelse(
      qualified_full == qualified_prepandemic,
      "Same",
      "Different"
    )
  ) %>%
  arrange(climate_var, region, admin1)

cat("Phase coherence comparison table created:",
    nrow(phase_coherence_comparison), "rows\n")

# Table 3: DLNM Cumulative RR Comparison (for provinces with both results)
# Extract cumulative RR at wavelet-derived lags for both periods

# For full period (2010-2024) - create from base objects
# First, get qualifying provinces from full period
qualifying_full_precip <- climate_dengue_phase_df %>%
  filter(climate_var == "precipitation_era",
         phase_coherence >= 0.85,
         round(phase_lag_months) >= 0) %>%
  mutate(lag = round(phase_lag_months)) %>%
  dplyr::select(admin1, climate_var, lag)

qualifying_full_temp <- climate_dengue_phase_df %>%
  filter(climate_var == "temperature_era",
         phase_coherence >= 0.85,
         round(phase_lag_months) >= 0) %>%
  mutate(lag = round(phase_lag_months)) %>%
  dplyr::select(admin1, climate_var, lag)

# Combine
qualifying_full <- bind_rows(qualifying_full_precip, qualifying_full_temp) %>% 
  mutate(climate_var = case_when(
    climate_var == "precipitation_era" ~ "Precipitation",
    climate_var == "temperature_era" ~ "Temperature",
    .default = "Rel. Humidity"
  ))

# Get DLNM results for these provinces at their lags
cumul_RR_full_for_comparison <- list()

for (i in seq_len(nrow(qualifying_full))) {
  prov <- qualifying_full$admin1[i]
  clim_variable <- qualifying_full$climate_var[i]
  lag_val <- qualifying_full$lag[i]
  
  # Find province index
  prov_idx <- which(admin1_indonesia$admin1 == prov & admin1_indonesia$idadmin1 != 99)[1]
  
  if (!is.na(prov_idx) && length(dlnm_output_list) >= prov_idx) {
    # Get cumulative RR for this province/climate/lag
    cumul_RR <- dlnm_output_list[[prov_idx]]$cumulative_RR %>%
      filter(clim_var == clim_variable, percentile == 90, lag == lag_val)
    
    if (nrow(cumul_RR) > 0) {
      cumul_RR_full_for_comparison[[length(cumul_RR_full_for_comparison) + 1]] <-
        tibble(
          province = prov,
          climate_var = clim_variable,
          lag_full = lag_val,
          RR_full = cumul_RR$cumulRR,
          RRlow_full = cumul_RR$cumulRR_lo,
          RRhigh_full = cumul_RR$cumulRR_hi,
          significant_full = if_else((cumul_RR$cumulRR_lo > 1 & cumul_RR$cumulRR_hi > 1) | (cumul_RR$cumulRR_lo < 1 & cumul_RR$cumulRR_hi < 1), 
                                     "Yes", "No")
        )
    }
  }
}

cumul_RR_full_for_comparison <- bind_rows(cumul_RR_full_for_comparison)

# For pre-pandemic (2010-2019)
cumul_RR_prepandemic_for_comparison <- list()

# Extract for both precipitation and temperature
for (prov in names(dlnm_prepandemic)) {
  
  # Precipitation
  if (prov %in% qualifying_full_precip$admin1) {
    # Get wavelet-derived lag for this province
    phase_lag <- climate_dengue_phase_df_prepandemic %>%
      filter(province == prov, climate_var == "Precipitation") %>%
      pull(phase_lag_months)
    
    if (length(phase_lag) > 0) {
      lag_rounded <- round(phase_lag)
      
      # Get cumulative RR at this lag (90th percentile)
      cumul_RR <- dlnm_prepandemic[[prov]]$cumulative_RR %>%
        filter(clim_var == "Precipitation", percentile == 90, lag == lag_rounded)
      
      if (nrow(cumul_RR) > 0) {
        cumul_RR_prepandemic_for_comparison[[paste0(prov, "_Precipitation")]] <-
          tibble(
            province = prov,
            climate_var = "Precipitation",
            lag_prepandemic = lag_rounded,
            RR_prepandemic = cumul_RR$cumulRR,
            RRlow_prepandemic = cumul_RR$cumulRR_lo,
            RRhigh_prepandemic = cumul_RR$cumulRR_hi,
            significant_prepandemic = if_else((cumul_RR$cumulRR_lo > 1 & cumul_RR$cumulRR_hi > 1) | (cumul_RR$cumulRR_lo < 1 & cumul_RR$cumulRR_hi < 1), 
                                              "Yes", "No")
          )
      }
    }
  }
  
  # Temperature
  if (prov %in% qualifying_full_temp$admin1) {
    # Get wavelet-derived lag for this province
    phase_lag <- climate_dengue_phase_df_prepandemic %>%
      filter(province == prov, climate_var == "Temperature") %>%
      pull(phase_lag_months)
    
    if (length(phase_lag) > 0) {
      lag_rounded <- round(phase_lag)
      
      # Get cumulative RR at this lag (90th percentile)
      cumul_RR <- dlnm_prepandemic[[prov]]$cumulative_RR %>%
        filter(clim_var == "Temperature", percentile == 90, lag == lag_rounded)
      
      if (nrow(cumul_RR) > 0) {
        cumul_RR_prepandemic_for_comparison[[paste0(prov, "_Temperature")]] <-
          tibble(
            province = prov,
            climate_var = "Temperature",
            lag_prepandemic = lag_rounded,
            RR_prepandemic = cumul_RR$cumulRR,
            RRlow_prepandemic = cumul_RR$cumulRR_lo,
            RRhigh_prepandemic = cumul_RR$cumulRR_hi,
            significant_prepandemic = if_else((cumul_RR$cumulRR_lo > 1 & cumul_RR$cumulRR_hi > 1) | (cumul_RR$cumulRR_lo < 1 & cumul_RR$cumulRR_hi < 1), 
                                              "Yes", "No")
          )
      }
    }
  }
}

cumul_RR_prepandemic_df <- bind_rows(cumul_RR_prepandemic_for_comparison)

# temperature 1 province mismatch: pre = Yogyakarta, all = Sulteng
# precipitation 3 province mismatch: pre = Aceh, Sumbar, all = Kalteng

# Combine full and prepandemic DLNM results
dlnm_RR_comparison <- cumul_RR_full_for_comparison %>%
  full_join(
    cumul_RR_prepandemic_df,
    by = c("province", "climate_var")
  ) %>%
  drop_na() %>% 
  mutate(
    RR_diff = RR_full - RR_prepandemic,
    both_significant = (significant_full == significant_prepandemic),
    effect_direction_agreement = ifelse(
      is.na(RR_full) | is.na(RR_prepandemic),
      NA,
      ifelse(
        (RR_full > 1 & RR_prepandemic > 1) | (RR_full < 1 & RR_prepandemic < 1),
        "Same direction",
        "Different direction"
      )
    )
  ) %>%
  arrange(climate_var, province)

cat("DLNM RR comparison table created:", nrow(dlnm_RR_comparison), "rows\n")
#### S8: Comparison of phase lags estimates with full data and excluding period during and after COVID-19 pandemic
table_s8 <- phase_lags_comparison

#### S9: Comparison of phase coherence with full data and excluding period during and after COVID-19 pandemic
table_s9 <- phase_coherence_comparison

#### S10: Comparison of DLNM cumulative RR with full data and excluding period during and after COVID-19 pandemic
table_s10 <- dlnm_RR_comparison

################################################################################
# SECTION 8: Summary Statistics
################################################################################

cat("\n=== SUMMARY OF COVID-19 SENSITIVITY ANALYSIS ===\n\n")

cat("1. PHASE LAGS (dengue timing relative to other provinces):\n")
cat("   - Provinces with same direction:",
    sum(phase_lags_comparison$qualitative_agreement == "Same direction"), "/",
    nrow(phase_lags_comparison), "\n")
cat("   - Mean absolute difference:",
    round(mean(abs(phase_lags_comparison$phase_lag_diff), na.rm = TRUE), 2),
    "months\n\n")

cat("2. PHASE COHERENCE (Precipitation):\n")
precip_coherence <- phase_coherence_comparison %>%
  filter(climate_var == "Precipitation")
cat("   - Provinces qualifying in both periods:",
    sum(precip_coherence$qualification_agreement == "Same" &
          precip_coherence$qualified_full), "\n")
cat("   - Provinces qualifying in full period only:",
    sum(precip_coherence$qualified_full & !precip_coherence$qualified_prepandemic), "\n")
cat("   - Provinces qualifying in prepandemic only:",
    sum(!precip_coherence$qualified_full & precip_coherence$qualified_prepandemic), "\n")
cat("   - Mean coherence difference:",
    round(mean(precip_coherence$coherence_diff, na.rm = TRUE), 3), "\n\n")

cat("3. PHASE COHERENCE (Temperature):\n")
temp_coherence <- phase_coherence_comparison %>%
  filter(climate_var == "Temperature")
cat("   - Provinces qualifying in both periods:",
    sum(temp_coherence$qualification_agreement == "Same" &
          temp_coherence$qualified_full), "\n")
cat("   - Provinces qualifying in full period only:",
    sum(temp_coherence$qualified_full & !temp_coherence$qualified_prepandemic), "\n")
cat("   - Provinces qualifying in prepandemic only:",
    sum(!temp_coherence$qualified_full & temp_coherence$qualified_prepandemic), "\n")
cat("   - Mean coherence difference:",
    round(mean(temp_coherence$coherence_diff, na.rm = TRUE), 3), "\n\n")

cat("4. DLNM CUMULATIVE RR:\n")
cat("   - Provinces with both significant elevated risks:",
    sum(dlnm_RR_comparison$both_significant, na.rm = TRUE), "\n")
cat("   - Provinces with same effect direction:",
    sum(dlnm_RR_comparison$effect_direction_agreement == "Same direction",
        na.rm = TRUE), "/",
    sum(!is.na(dlnm_RR_comparison$effect_direction_agreement)), "\n")
cat("   - Mean RR difference (where both available):",
    round(mean(dlnm_RR_comparison$RR_diff, na.rm = TRUE), 2), "\n\n")

cat("CONCLUSION:\n")
cat("Results are qualitatively ",
    ifelse(
      sum(phase_lags_comparison$qualitative_agreement == "Same direction") > 30 &&
        sum(precip_coherence$qualification_agreement == "Same") > 15,
      "SIMILAR",
      "DIFFERENT"
    ),
    " between full (2010-2024) and pre-pandemic (2010-2019) periods.\n")
cat("This suggests findings are ",
    ifelse(
      sum(phase_lags_comparison$qualitative_agreement == "Same direction") > 30 &&
        sum(precip_coherence$qualification_agreement == "Same") > 15,
      "ROBUST",
      "SENSITIVE"
    ),
    " to the inclusion of pandemic years.\n\n")

#### S11: Sensitivity analysis of threshold for phase coherence
table_s11 <- phase_significance_threshold_sensitivity

# The higher the threshold, fewer provinces selected. Highest threshold of 0.9 gave 16 provinces for precipitation
# and 5 for temperature (18 and 7 for 0.85 threshold and 22 and 9 for 0.7 threshold). Significant provinces stable
# at 0.85 threshold (11 and 3) and it is 12 and 3 for 0.7 threshold. Less significant provinces in 0.9 threshold
# (9 and 2)

#### Figures
#### S1: Geographical regions of Indonesia.
#### Create supplementary figure of regions of Indonesia
region_centroids <- region_shp %>%
  st_centroid() %>%
  cbind(st_coordinates(.)) %>%
  mutate(region = factor(region_shp$region,levels=c("Sumatra","Java & Bali","Kalimantan","Nusa Tenggara",
                                                    "Sulawesi","Maluku","Papua"))) %>%
  st_drop_geometry()

region_shp$region <- factor(region_shp$region,levels=c("Sumatra","Java & Bali","Kalimantan","Nusa Tenggara",
                                                       "Sulawesi","Maluku","Papua"))

fig_s1 <- region_shp %>% 
  ggplot() +
  geom_sf(data = other_shp, fill = "gray80", linewidth = 0.2) +  
  geom_sf(aes(fill = region), color = 'black', linewidth = 0.2,
          show.legend = c(fill = TRUE, color = FALSE)) +
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

#### S2: Between province phase lags statistics, including 2.5th-97.5th percentile ranges of pairwise differences.
fig_s2 <- province_phase_stats %>%
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
    col = NULL,
    title = "Province phase lags with confidence intervals"
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

#### S3: Wavelet-based lag vs DLNM-based lag. A) DLNM-based lag with ’highest-risk’ approach; B) DLNM-based lag with ’first-significant’ approach.
### supp
phase_heatmap_data_pct90_lag_compare_maximum <- 
  phase_heatmap_data_over85 %>% left_join(cumul_RR_highest_maximum %>% 
                                            filter(percentile==90) %>% 
                                            dplyr::select(idadmin1,clim_var,lag_highest))

phase_heatmap_data_pct90_lag_compare_first <- 
  phase_heatmap_data_over85 %>% left_join(cumul_RR_highest_first %>% 
                                            filter(percentile==90) %>% 
                                            dplyr::select(idadmin1,clim_var,lag_highest))

plot_compare_lag_pre_maximum <- phase_heatmap_data_pct90_lag_compare_maximum %>% filter(clim_var == "Precipitation") %>% 
  ggplot(aes(x=phase_lag_months,y=lag_highest)) +
  geom_point() +
  geom_abline(col="red",linetype=2) +
  ylim(0,4) +
  xlim(0,4) +
  theme_Publication(base_size = 10) +
  labs(x="Wavelet-based lag",y="DLNM-based lag",title="Precipitation",tag="A")

plot_compare_lag_tem_maximum <- phase_heatmap_data_pct90_lag_compare_maximum %>% filter(clim_var == "Temperature") %>% 
  ggplot(aes(x=phase_lag_months,y=lag_highest)) +
  geom_point() +
  geom_abline(col="red",linetype=2) +
  ylim(0,4) +
  xlim(0,4) +
  theme_Publication(base_size = 10) +
  labs(x="Wavelet-based lag",y="DLNM-based lag",title="Temperature")

plot_compare_lag_pre_first <- phase_heatmap_data_pct90_lag_compare_first %>% filter(clim_var == "Precipitation") %>% 
  ggplot(aes(x=phase_lag_months,y=lag_highest)) +
  geom_point() +
  geom_abline(col="red",linetype=2) +
  ylim(0,4) +
  xlim(0,4) +
  theme_Publication(base_size = 10) +
  labs(x="Wavelet-based lag",y="DLNM-based lag",title="Precipitation",tag="B")

plot_compare_lag_tem_first <- phase_heatmap_data_pct90_lag_compare_first %>% filter(clim_var == "Temperature") %>% 
  ggplot(aes(x=phase_lag_months,y=lag_highest)) +
  geom_point() +
  geom_abline(col="red",linetype=2) +
  ylim(0,4) +
  xlim(0,4) +
  theme_Publication(base_size = 10) +
  labs(x="Wavelet-based lag",y="DLNM-based lag",title="Temperature")

plot_compare_lag_maximum <- plot_compare_lag_pre_maximum | plot_compare_lag_tem_maximum
plot_compare_lag_first <- plot_compare_lag_pre_first | plot_compare_lag_tem_first
plot_compare_lag <- plot_compare_lag_maximum / plot_compare_lag_first

#### S4: Year-by-Year Peak Month Variability (all provinces)
cat("\n\n================================================================================\n")
cat("FIGURE S6: Year-by-Year Peak Month Variability\n")
cat("================================================================================\n\n")

cat("Creating heatmap showing peak months by province and epidemic year...\n")

# Extract year-by-year peak months for all provinces
peak_months_by_year_list <- list()

for (i in seq_len(34)) {
  
  # Get province info
  prov_info <- admin1_indonesia[i,]
  
  # Extract peak months by year
  peak_months_result <- find_peak_months(province_wavelet$phase_angles[[i]])
  
  # Get peak data
  peak_data_yearly <- peak_months_result$peak_data %>%
    filter(!Is_Anomalous) %>%  # Exclude anomalous periods
    mutate(
      province = prov_info$admin1,
      idadmin1 = prov_info$idadmin1,
      region = prov_info$region
    ) %>%
    dplyr::select(province, idadmin1, region, Period, Month_Name, Month_JanDec)
  
  peak_months_by_year_list[[i]] <- peak_data_yearly
}

# Combine all provinces
peak_months_by_year_df <- bind_rows(peak_months_by_year_list)

# Get longitude for ordering provinces (west to east)
province_longitude <- admin1_shp %>%
  st_centroid() %>%
  cbind(st_coordinates(.)) %>%
  st_drop_geometry() %>%
  dplyr::select(shapeName, longitude = X) %>%
  left_join(admin1_EN, by = "shapeName") %>%
  arrange(longitude)

# Join with peak month data and order by longitude
peak_months_for_plot <- peak_months_by_year_df %>%
  left_join(admin1_EN, by = "idadmin1") %>%
  left_join(province_longitude %>% dplyr::select(idadmin1, longitude), by = "idadmin1") %>%
  arrange(longitude) %>%
  mutate(
    shapeName = factor(shapeName, levels = rev(province_longitude$shapeName)),
    Month_Name = factor(Month_Name,
                        levels = c("July","August","September","October", "November", "December",
                                   "January", "February", "March", "April", "May","June"))
  )

# Create heatmap
fig_s4 <- ggplot(peak_months_for_plot,
                              aes(x = Period, y = shapeName, fill = Month_Name)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_d(
    option = "turbo",
    na.value = "grey90",
    name = "Peak Month"
  ) +
  labs(
    title = "Year-to-Year Variability in Peak Outbreak Timing",
    subtitle = "Provinces ordered by longitude (west to east)",
    x = "Epidemic Year (July-June)",
    y = NULL
  ) +
  theme_Publication() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey40")
  ) +
  guides(fill = guide_legend(ncol = 1))

#### S5: Year-by-Year Peak Month Variability (western provinces/Sumatra, Java & Bali)
fig_s5 <- peak_months_for_plot %>%
  filter(region %in% c("SUMATRA","JAVA & BALI")) %>%
  ggplot(aes(x = Period, y = shapeName, fill = Month_Name)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_d(
    option = "turbo",
    na.value = "grey90",
    name = "Peak Month"
  ) +
  labs(
    title = "Year-to-Year Variability in Peak Outbreak Timing",
    subtitle = "Provinces ordered by longitude (west to east)",
    x = "Epidemic Year (July-June)",
    y = NULL
  ) +
  theme_Publication() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey40")
  ) +
  guides(fill = guide_legend(ncol = 1))

################################################################################
# Summary of Generated Objects
################################################################################

cat("\n\n================================================================================\n")
cat("SUPPLEMENTARY MATERIALS GENERATION COMPLETE\n")
cat("================================================================================\n\n")

cat("TABLES CREATED:\n")
cat("  - table_s1:  Province reference index\n")
cat("  - table_s2:  Annual incidence by epidemic year\n")
cat("  - table_s3:  Silhouette statistics (all provinces)\n")
cat("  - table_s4:  Silhouette statistics (western + Kalimantan)\n")
cat("  - table_s5:  Silhouette statistics (western only)\n")
cat("  - table_s6:  Province dissimilarity within regions\n")
cat("  - table_s7:  Regional climate index correlations\n")
cat("  - table_s8:  COVID sensitivity - phase lags\n")
cat("  - table_s9:  COVID sensitivity - phase coherence\n")
cat("  - table_s10: COVID sensitivity - DLNM RR\n")
cat("  - table_s11: Threshold sensitivity analysis\n\n")

cat("FIGURES CREATED:\n")
cat("  - fig_s1: Geographical regions of Indonesia\n")
cat("  - fig_s2: Province phase lags with confidence intervals\n")
cat("  - fig_s3: Wavelet vs DLNM lag comparison (plot_compare_lag)\n")
cat("  - fig_s4: Year-by-year peak month variability (all provinces)\n")
cat("  - fig_s5: Year-by-year peak month variability (western provinces)\n\n")

cat("All objects are available in the global environment.\n")
cat("To render the supplementary materials PDF, run:\n")
cat("  quarto render supplementary-materials.qmd\n\n")

