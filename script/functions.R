################################################################################
# Custom Functions for Spatiotemporal Dengue-Climate Analysis
#
# This file contains all custom functions used in the analysis:
#
# VISUALIZATION:
#   - theme_Publication(): Clean ggplot2 theme for manuscripts
#   - scale_fill_Publication(): Colorblind-friendly fill palette
#   - scale_colour_Publication(): Colorblind-friendly color palette
#
# DATA PREPARATION:
#   - impute_incidence_rates(): Hierarchical median imputation
#   - prepare_for_wavelet(): Format data for wavelet analysis
#
# WAVELET ANALYSIS:
#   - conduct_wavelet_analysis(): CWT with Morlet wavelet
#   - mean_circular_diff(): Circular difference calculation
#   - create_phase_matrix(): Pairwise phase differences
#   - calculate_phase_lags(): Average lag per province
#   - calculate_phase_stats(): Phase lag with 95% CI
#   - find_peak_months(): Peak month with anomaly detection
#   - fiscal_year_vector(): July-June fiscal year labels
#
# CLIMATE ANALYSIS:
#   - prepare_wavelet_data(): Format for wavelet coherence
#   - calculate_correlation(): Bootstrap correlation with CI
#
# DLNM:
#   - dlnm_run(): Fit DLNM for 5 climate variables
#   - dlnm_significance_lag(): Extract significant effects
#   - dlnm_lag_visualise(): Generate lag-response curves
#   - dlnm_contour_visualise(): Generate contour plots
#
# NOTE: All packages are loaded in data-processing.R
################################################################################

################################################################################
# 1. VISUALIZATION FUNCTIONS
################################################################################

#' Publication-quality ggplot2 theme
#' @param base_size Base font size (default: 12)
#' @param base_family Font family (default: "sans")
#' @return A ggplot2 theme object
#' @source https://rpubs.com/koundy/71792
theme_Publication <- function(base_size=12, base_family="sans") {
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold", hjust = 0.5),
            # size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(1.5,1.0,1.0,1.0),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))

}

scale_fill_Publication <- function(...){
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

################################################################################
# 2. DATA PREPARATION FUNCTIONS
################################################################################

#' Impute missing incidence rates with median of same month in other years
impute_incidence_rates <- function(data) {
  # calculate medians by district and month
  median_rates <- data %>%
    filter(!is.na(incidence_rate)) %>%
    group_by(idadmin1, idadmin2, month_factor) %>%
    summarise(
      median_incidence = median(incidence_rate, na.rm = TRUE),
      .groups = "drop"
    )

  # for cases where we don't have any data for a specific district-month combination
  # calculate medians by province and month as a fallback
  province_median_rates <- data %>%
    filter(!is.na(incidence_rate)) %>%
    group_by(idadmin1, month_factor) %>%
    summarise(
      province_median_incidence = median(incidence_rate, na.rm = TRUE),
      .groups = "drop"
    )

  # for extreme cases, calculate country-level medians by month
  country_median_rates <- data %>%
    filter(!is.na(incidence_rate)) %>%
    group_by(month_factor) %>%
    summarise(
      country_median_incidence = median(incidence_rate, na.rm = TRUE),
      .groups = "drop"
    )

  # perform the imputation with a multi-level fallback strategy
  imputed_data <- data %>%
    # join with district-level medians
    left_join(median_rates, by = c("idadmin1", "idadmin2", "month_factor")) %>%
    # join with province-level medians
    left_join(province_median_rates, by = c("idadmin1", "month_factor")) %>%
    # join with country-level medians
    left_join(country_median_rates, by = "month_factor") %>%
    # impute using the most specific available median: i.e, if district available, use district
    mutate(
      imputed_incidence_rate = case_when(
        !is.na(incidence_rate) ~ incidence_rate,
        !is.na(median_incidence) ~ median_incidence,
        !is.na(province_median_incidence) ~ province_median_incidence,
        !is.na(country_median_incidence) ~ country_median_incidence,
        TRUE ~ 0  # last resort if no medians available
      ),
      imputation_source = case_when(
        !is.na(incidence_rate) ~ "original",
        !is.na(median_incidence) ~ "district_median",
        !is.na(province_median_incidence) ~ "province_median",
        !is.na(country_median_incidence) ~ "country_median",
        TRUE ~ "zero_default"
      )
    ) %>%
    # calculate imputed cases
    mutate(
      imputed_cases = case_when(
        !is.na(cases) ~ cases,
        TRUE ~ round(imputed_incidence_rate * pop / 100000)
      )
    ) %>%
    # clean up unused columns
    dplyr::select(-median_incidence, -province_median_incidence, -country_median_incidence)

  return(imputed_data)
}

#' Prepare time series for wavelet analysis
#' @param data Data frame with time series data
#' @param group_var Column name for grouping (e.g., "admin1")
#' @param var Column name for the variable to analyze
#' @param min_cases Minimum total cases to include group (default: 100)
prepare_for_wavelet <- function(data, group_var, var, min_cases = 100) {
  # create wide format for time series analysis
  ts_data <- data %>%
    dplyr::select(!!sym(group_var), date, !!sym(var)) %>%
    pivot_wider(
      names_from = !!sym(group_var),
      values_from = !!sym(var)
    )
  
  # fill any missing values (0 cases)
  ts_data[is.na(ts_data)] <- 0
  
  valid_groups <- data %>%
    pull(!!sym(group_var)) %>% unique()
  
  return(list(
    ts_data = ts_data,
    valid_groups = valid_groups
  ))
}

################################################################################
# 3. WAVELET ANALYSIS FUNCTIONS
################################################################################

#' Conduct wavelet analysis on a time series
conduct_wavelet_analysis <- function(ts_data, group_names) {
  # extract date column
  dates <- ts_data$date
  n_dates <- length(dates)

  # create empty lists to store results
  wt_results <- list()
  phase_angles <- list()

  # loop through each region/province/district
  for (i in 1:length(group_names)) {
    group <- group_names[i]

    # extract and prepare time series
    ts <- ts_data[[group]]

    # perform wavelet transform
    wt <- wt(cbind(1:n_dates, as.numeric(ts)),
             dt = 1, # 1 month time step
             mother = "morlet")

    # store all wavelet analysis results
    wt_results[[group]] <- wt

    # find period closest to 12 months (annual cycle)
    annual_idx <- which.min(abs(wt$period - 12))

    # extract phase angles for annual component
    phase_angles[[group]] <- wt$phase[annual_idx, ]
  }

  return(list(
    wt_results = wt_results,
    phase_angles = phase_angles,
    dates = dates
  ))
}

#' Calculate mean circular difference between phase angles
mean_circular_diff <- function(phase_a, phase_b) {
  # calculate circular differences
  diff <- phase_b - phase_a
  # adjust to be in range -pi to pi
  diff <- (diff + pi) %% (2*pi) - pi
  # return mean
  return(mean(diff, na.rm = TRUE))
}

#' Create phase difference matrix for a set of areas
create_phase_matrix <- function(phase_angles, group_names) {
  n_groups <- length(group_names)

  # initialize matrices
  phase_diffs <- matrix(0, nrow = n_groups, ncol = n_groups)
  colnames(phase_diffs) <- rownames(phase_diffs) <- group_names

  # calculate all pairwise phase differences
  for (i in 1:n_groups) {
    for (j in 1:n_groups) {
      if (i != j) {
        phase_diff <- mean_circular_diff(
          phase_angles[[group_names[i]]],
          phase_angles[[group_names[j]]]
        )
        phase_diffs[i, j] <- phase_diff
      }
    }
  }

  # convert phase differences to months (2*pi = 12 months)
  phase_months <- phase_diffs * 12 / (2*pi)

  return(phase_months)
}

#' Calculate average phase lag for each area
calculate_phase_lags <- function(phase_matrix) {
  # for each area, calculate average phase difference from all other areas
  n_areas <- nrow(phase_matrix)
  phase_lags <- numeric(n_areas)
  names(phase_lags) <- rownames(phase_matrix)

  for (i in 1:n_areas) {
    # average phase difference from all other areas
    phase_lags[i] <- mean(phase_matrix[i, -i], na.rm = TRUE)
  }

  return(phase_lags)
}

#' Calculate Phase Lag Statistics for Each Province
#'
#' Calculates the median phase lag and the range of phase differences
#' for each province relative to all other provinces.
#'
#' @param phase_matrix Matrix of pairwise phase differences (in months)
#' @return Data frame with median phase lag and 2.5th-97.5th percentile range
#'
#' @note The lower_ci and upper_ci columns represent the 2.5th and 97.5th
#'   percentiles of phase differences with other provinces, NOT statistical
#'   confidence intervals. They indicate the range of relative timing
#'   differences across the country, not uncertainty in the estimate.
calculate_phase_stats <- function(phase_matrix) {
  n_provinces <- nrow(phase_matrix)
  province_names <- rownames(phase_matrix)

  # create data frame to store results
  results <- data.frame(
    province = character(n_provinces),
    median_phase_lag = numeric(n_provinces),
    lower_ci = numeric(n_provinces),  # 2.5th percentile of phase differences
    upper_ci = numeric(n_provinces),  # 97.5th percentile of phase differences
    stringsAsFactors = FALSE
  )

  # calculate statistics for each province
  for (i in 1:n_provinces) {
    # get phase differences with all other provinces
    phase_diffs <- phase_matrix[i, -i]  # Exclude self-comparison

    # calculate median and range (2.5th-97.5th percentile)
    results$province[i] <- province_names[i]
    results$median_phase_lag[i] <- median(phase_diffs, na.rm = TRUE)
    results$lower_ci[i] <- quantile(phase_diffs, 0.025, na.rm = TRUE)
    results$upper_ci[i] <- quantile(phase_diffs, 0.975, na.rm = TRUE)
  }

  # sort provinces by median phase lag
  results <- results[order(results$median_phase_lag), ]

  # convert province to factor with levels in the sorted order for plotting
  results$province <- factor(results$province, levels = results$province)

  return(results)
}

################################################################################
# 4. DISTANCE AND CORRELATION FUNCTIONS
################################################################################

#' Calculate distance matrix between locations' centroids in kilometers
calculate_distance_km <- function(lon1, lat1, lon2, lat2) {
  # convert degrees to radians
  lon1_rad <- lon1 * pi / 180
  lat1_rad <- lat1 * pi / 180
  lon2_rad <- lon2 * pi / 180
  lat2_rad <- lat2 * pi / 180

  # earth radius in kilometers
  R <- 6371

  # haversine formula
  dlon <- lon2_rad - lon1_rad
  dlat <- lat2_rad - lat1_rad
  a <- sin(dlat/2)^2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  distance <- R * c

  return(distance)
}

#' Fit spline to correlation vs distance data
fit_spline_covariance <- function(data, corr_col) {
  #### sort by distance
  data <- data[order(data$distance), ]

  #### create bins for distance
  n_bins <- 30
  bin_size <- max(data$distance, na.rm = TRUE) / n_bins
  data$bin <- cut(data$distance, breaks = seq(0, max(data$distance, na.rm = TRUE) + bin_size, by = bin_size))

  #### calculate mean correlation by bin
  bin_summary <- data %>%
    group_by(bin) %>%
    summarize(
      mean_dist = mean(distance, na.rm = TRUE),
      mean_corr = mean(!!sym(corr_col), na.rm = TRUE),
      sd_corr = sd(!!sym(corr_col), na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(mean_dist), !is.na(mean_corr))

  #### fit a smooth spline
  spline_fit <- smooth.spline(bin_summary$mean_dist, bin_summary$mean_corr,
                              spar = 0.5)  # Adjust smoothing parameter as needed

  #### generate prediction grid
  pred_grid <- seq(0, max(data$distance, na.rm = TRUE), length.out = 200)
  pred_values <- predict(spline_fit, pred_grid)

  #### calculate country-wide correlation
  country_corr <- mean(data[[corr_col]], na.rm = TRUE)

  #### find where correlation drops below country-wide level
  if(any(pred_values$y < country_corr)) {
    idx <- which(pred_values$y < country_corr)[1]
    crossover_distance <- pred_values$x[idx]
  } else {
    crossover_distance <- NA
  }

  #### bootstrap confidence intervals
  boot_spline <- function(data, i) {
    boot_data <- data[i, ]
    boot_data <- boot_data[order(boot_data$distance), ]

    #### bin the data
    boot_data$bin <- cut(boot_data$distance,
                         breaks = seq(0, max(boot_data$distance, na.rm = TRUE) + bin_size, by = bin_size))

    #### calculate mean by bin
    boot_summary <- boot_data %>%
      group_by(bin) %>%
      summarize(
        mean_dist = mean(distance, na.rm = TRUE),
        mean_corr = mean(!!sym(corr_col), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(!is.na(mean_dist), !is.na(mean_corr))

    #### fit spline
    boot_spline <- try(smooth.spline(boot_summary$mean_dist, boot_summary$mean_corr, spar = 0.5),
                       silent = TRUE)

    if(inherits(boot_spline, "try-error")) {
      return(rep(NA, length(pred_grid)))
    } else {
      return(predict(boot_spline, pred_grid)$y)
    }
  }

  #### run bootstrap
  set.seed(123)
  boot_results <- boot(data, boot_spline, R = 500)

  #### calculate confidence intervals
  ci_matrix <- matrix(NA, nrow = length(pred_grid), ncol = 2)
  for(i in 1:length(pred_grid)) {
    ci <- boot.ci(boot_results, type = "perc", index = i)
    if(!inherits(ci, "try-error") && !is.null(ci$percent)) {
      ci_matrix[i, ] <- ci$percent[4:5]
    }
  }

  #### return results
  return(list(
    distances = pred_grid,
    mean_correlation = pred_values$y,
    lower_ci = ci_matrix[, 1],
    upper_ci = ci_matrix[, 2],
    country_correlation = country_corr,
    crossover_distance = crossover_distance
  ))
}

#' Create correlation vs distance plot
create_correlation_plot <- function(spline_data, title, y_lab) {
  # Create dataframe for plotting
  plot_data <- data.frame(
    distance = spline_data$distances,
    correlation = spline_data$mean_correlation,
    lower_ci = spline_data$lower_ci,
    upper_ci = spline_data$upper_ci
  )

  # Create plot
  p <- ggplot(plot_data, aes(x = distance, y = correlation)) +
    geom_point(data=dist_corr_df,aes(x = distance, y = coherence), col="darkgray", alpha = 0.25) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "gray80", alpha = 0.5) +
    geom_line(color = "blue", size = 1) +
    geom_hline(yintercept = spline_data$country_correlation, color = "red", linetype = "solid", size = 1) +
    geom_vline(xintercept = spline_data$crossover_distance, linetype = "solid") +
    geom_hline(yintercept = 0, color = "red", linetype = 2) +
    annotate("text", x = spline_data$crossover_distance, y = 0.9,
             label = round(spline_data$crossover_distance), hjust = -0.3) +
    labs(
      # title = title,
      x = "Distance (km)",
      y = y_lab
    ) +
    ylim(-0.5, 1) +
    scale_y_continuous(breaks = c(-0.5,-0.25,0,0.25,0.5,0.75,1),limits=c(-0.5,1)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(fill = NA, color = "black"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14)
    )

  return(p)
}

#' Prepare time series data for wavelet coherence analysis
prepare_wavelet_data <- function(merged_data, province_name, climate_var,
                                 incidence_var = "incidence_rate") {

  # Filter data for specific province
  province_data <- merged_data %>%
    filter(admin1 == province_name) %>%
    arrange(date) %>%
    dplyr::select(date, all_of(c(incidence_var, climate_var))) %>%
    na.omit()

  # Create time vector (decimal years for biwavelet)
  province_data$time <- decimal_date(province_data$date)

  # Create matrix format required by biwavelet
  # Column 1: time, Column 2: variable 1, Column 3: variable 2
  wavelet_matrix <- matrix(
    c(province_data$time,
      province_data[[incidence_var]],
      province_data[[climate_var]]),
    ncol = 3
  )

  return(list(
    data = wavelet_matrix,
    province = province_name,
    climate_var = climate_var,
    incidence_var = incidence_var
  ))
}

#' Calculate correlation with confidence intervals
calculate_correlation <- function(x, y, method = "pearson", boot_sample = 1000) {
  if (all(is.na(x)) || all(is.na(y)) || length(x) < 5 || length(y) < 5) {
    return(c(cor = NA, p_value = NA, lower_ci = NA, upper_ci = NA))
  }

  # Remove NAs
  complete_cases <- complete.cases(x, y)
  x <- x[complete_cases]
  y <- y[complete_cases]

  # Calculate correlation
  cor_test <- cor.test(x, y, method = method)
  cor_test_ci <- ci_cor(x, y, method = method, type = "bootstrap", R = boot_sample, seed = 1234)

  return(c(
    cor = cor_test$estimate,
    p_value = cor_test$p.value,
    lower_ci = cor_test_ci$interval[1],
    upper_ci = cor_test_ci$interval[2]
  ))
}

################################################################################
# 5. DLNM FUNCTIONS
################################################################################

#' Extract significant DLNM effects for each climate variable
dlnm_significance_lag <- function(dlnm_output, data_input){

  idadmin1 = unique(data_input$idadmin1)

  #### Median
  #### ONI
  cen_oni = 0

  #### DMI/IOD
  cen_dmi = 0

  #### precipitation
  cen_pre = round(median(data_input$precipitation),-1) # to the nearest tens

  #### temperature
  cen_tem = round(median(data_input$temperature))

  #### relative humidity
  cen_hum = round(median(data_input$rel_humidity))

  #### Filter
  #### ONI
  significance_pos_oni <- dlnm_output$results_oni %>%
    filter(lag_months %in% 0:6) %>%
    filter((round(RR_lower_95ci,1) > 1 & round(RR_upper_95ci,1) > 1)) %>%
    mutate(effect=abs(RR-1)) %>%
    filter(effect==max(effect)) %>%
    mutate(var_value=ifelse(var_value>cen_oni,"High","Low"))

  if(nrow(significance_pos_oni)==0){
    significance_pos_oni=tibble(idadmin1=idadmin1,clim_var="ONI",var_value=NA,lag_months=NA,
                                RR=NA,RR_lower_95ci=NA,RR_upper_95ci=NA,effect=NA)
  }

  significance_neg_oni <- dlnm_output$results_oni %>%
    filter(lag_months %in% 0:6) %>%
    filter((round(RR_lower_95ci,1) < 1 & round(RR_upper_95ci,1) < 1)) %>%
    mutate(effect=abs(RR-1)) %>%
    filter(effect==max(effect)) %>%
    mutate(var_value=ifelse(var_value>cen_oni,"High","Low"))

  if(nrow(significance_neg_oni)==0){
    significance_neg_oni=tibble(idadmin1=idadmin1,clim_var="ONI",var_value=NA,lag_months=NA,
                                RR=NA,RR_lower_95ci=NA,RR_upper_95ci=NA,effect=NA)
  }

  #### IOD/DMI
  significance_pos_dmi <- dlnm_output$results_dmi %>%
    filter(lag_months %in% 0:6) %>%
    filter((round(RR_lower_95ci,1) > 1 & round(RR_upper_95ci,1) > 1)) %>%
    mutate(effect=abs(RR-1)) %>%
    filter(effect==max(effect)) %>%
    mutate(var_value=ifelse(var_value>cen_dmi,"High","Low"))

  if(nrow(significance_pos_dmi)==0){
    significance_pos_dmi=tibble(idadmin1=idadmin1,clim_var="DMI",var_value=NA,lag_months=NA,
                                RR=NA,RR_lower_95ci=NA,RR_upper_95ci=NA,effect=NA)
  }

  significance_neg_dmi <- dlnm_output$results_dmi %>%
    filter(lag_months %in% 0:6) %>%
    filter((round(RR_lower_95ci,1) < 1 & round(RR_upper_95ci,1) < 1)) %>%
    mutate(effect=abs(RR-1)) %>%
    filter(effect==max(effect)) %>%
    mutate(var_value=ifelse(var_value>cen_dmi,"High","Low"))

  if(nrow(significance_neg_dmi)==0){
    significance_neg_dmi=tibble(idadmin1=idadmin1,clim_var="DMI",var_value=NA,lag_months=NA,
                                RR=NA,RR_lower_95ci=NA,RR_upper_95ci=NA,effect=NA)
  }

  #### precipitation
  significance_pos_pre <- dlnm_output$results_pre %>%
    filter(lag_months %in% 0:3) %>%
    filter((round(RR_lower_95ci,1) > 1 & round(RR_upper_95ci,1) > 1)) %>%
    mutate(effect=abs(RR-1)) %>%
    filter(effect==max(effect)) %>%
    mutate(var_value=ifelse(var_value>cen_pre,"High","Low"))

  if(nrow(significance_pos_pre)==0){
    significance_pos_pre=tibble(idadmin1=idadmin1,clim_var="Precipitation",var_value=NA,lag_months=NA,
                                RR=NA,RR_lower_95ci=NA,RR_upper_95ci=NA,effect=NA)
  }

  significance_neg_pre <- dlnm_output$results_pre %>%
    filter(lag_months %in% 0:3) %>%
    filter((round(RR_lower_95ci,1) < 1 & round(RR_upper_95ci,1) < 1)) %>%
    mutate(effect=abs(RR-1)) %>%
    filter(effect==max(effect)) %>%
    mutate(var_value=ifelse(var_value>cen_pre,"High","Low"))

  if(nrow(significance_neg_pre)==0){
    significance_neg_pre=tibble(idadmin1=idadmin1,clim_var="Precipitation",var_value=NA,lag_months=NA,
                                RR=NA,RR_lower_95ci=NA,RR_upper_95ci=NA,effect=NA)
  }

  #### temperature
  significance_pos_tem <- dlnm_output$results_tem %>%
    filter(lag_months %in% 0:3) %>%
    filter((round(RR_lower_95ci,1) > 1 & round(RR_upper_95ci,1) > 1)) %>%
    mutate(effect=abs(RR-1)) %>%
    filter(effect==max(effect)) %>%
    mutate(var_value=ifelse(var_value>cen_tem,"High","Low"))

  if(nrow(significance_pos_tem)==0){
    significance_pos_tem=tibble(idadmin1=idadmin1,clim_var="Temperature",var_value=NA,lag_months=NA,
                                RR=NA,RR_lower_95ci=NA,RR_upper_95ci=NA,effect=NA)
  }

  significance_neg_tem <- dlnm_output$results_tem %>%
    filter(lag_months %in% 0:3) %>%
    filter((round(RR_lower_95ci,1) < 1 & round(RR_upper_95ci,1) < 1)) %>%
    mutate(effect=abs(RR-1)) %>%
    filter(effect==max(effect)) %>%
    mutate(var_value=ifelse(var_value>cen_tem,"High","Low"))

  if(nrow(significance_neg_tem)==0){
    significance_neg_tem=tibble(idadmin1=idadmin1,clim_var="Temperature",var_value=NA,lag_months=NA,
                                RR=NA,RR_lower_95ci=NA,RR_upper_95ci=NA,effect=NA)
  }

  #### Relative humidity
  significance_pos_hum <- dlnm_output$results_hum %>%
    filter(lag_months %in% 0:3) %>%
    filter((round(RR_lower_95ci,1) > 1 & round(RR_upper_95ci,1) > 1)) %>%
    mutate(effect=abs(RR-1)) %>%
    filter(effect==max(effect)) %>%
    mutate(var_value=ifelse(var_value>cen_hum,"High","Low"))

  if(nrow(significance_pos_hum)==0){
    significance_pos_hum=tibble(idadmin1=idadmin1,clim_var="Relative humidity",var_value=NA,lag_months=NA,
                                RR=NA,RR_lower_95ci=NA,RR_upper_95ci=NA,effect=NA)
  }

  significance_neg_hum <- dlnm_output$results_hum %>%
    filter(lag_months %in% 0:3) %>%
    filter((round(RR_lower_95ci,1) < 1 & round(RR_upper_95ci,1) < 1)) %>%
    mutate(effect=abs(RR-1)) %>%
    filter(effect==max(effect)) %>%
    mutate(var_value=ifelse(var_value>cen_hum,"High","Low"))

  if(nrow(significance_neg_hum)==0){
    significance_neg_hum=tibble(idadmin1=idadmin1,clim_var="Relative humidity",var_value=NA,lag_months=NA,
                                RR=NA,RR_lower_95ci=NA,RR_upper_95ci=NA,effect=NA)
  }

  return(list(significance_pos_oni=significance_pos_oni,
              significance_pos_dmi=significance_pos_dmi,
              significance_pos_pre=significance_pos_pre,
              significance_pos_tem=significance_pos_tem,
              significance_pos_hum=significance_pos_hum,
              significance_neg_oni=significance_neg_oni,
              significance_neg_dmi=significance_neg_dmi,
              significance_neg_pre=significance_neg_pre,
              significance_neg_tem=significance_neg_tem,
              significance_neg_hum=significance_neg_hum))

}

#' Visualize lag-response curves for min and max climate values
dlnm_lag_visualise <- function(dlnm_output){

  #### ONI
  min_oni = min(dlnm_output$results_oni$var_value)
  max_oni = max(dlnm_output$results_oni$var_value)

  #### DMI/IOD
  min_dmi = min(dlnm_output$results_dmi$var_value)
  max_dmi = max(dlnm_output$results_dmi$var_value)

  #### precipitation
  min_pre = min(dlnm_output$results_pre$var_value)
  max_pre = max(dlnm_output$results_pre$var_value)

  #### temperature
  min_tem = min(dlnm_output$results_tem$var_value)
  max_tem = max(dlnm_output$results_tem$var_value)

  #### relative humidity
  min_hum = min(dlnm_output$results_hum$var_value)
  max_hum = max(dlnm_output$results_hum$var_value)

  #### Visualise
  #### ONI
  results_oni_filtered <- dlnm_output$results_oni %>%
    filter(var_value %in% c(min_oni,max_oni)) %>%
    mutate(var_value_char=factor(as.character(var_value),levels=c(as.character(c(min_oni,max_oni)))))

  plot_oni <- results_oni_filtered %>%
    ggplot(aes(x=lag_months,y=RR,col=var_value_char,fill=var_value_char)) +
    geom_ribbon(aes(ymin=RR_lower_95ci,ymax=RR_upper_95ci),alpha=0.15,col=NA) +
    geom_line() +
    geom_hline(yintercept=1,col="red",linetype=2) +
    theme_Publication(base_size=10) +
    labs(x="Lag (months)",y="Relative risk",col=NULL,fill=NULL) +
    theme(legend.position = "inside",legend.position.inside = c(0.5,0.8)) +
    scale_x_continuous(breaks=0:6) +
    scale_colour_discrete(palette=hcl.colors(2,"Blue-Red")) +
    scale_fill_discrete(palette=hcl.colors(2,"Blue-Red")) + theme(
      legend.background = element_rect(fill = "transparent", colour = NA), # Transparent legend box
      legend.key = element_rect(fill = "transparent", colour = NA),        # Transparent background for key symbols
      legend.box.background = element_rect(fill = "transparent", colour = NA) # Transparent legend panel
    )

  #### IOD/DMI
  results_dmi_filtered <- dlnm_output$results_dmi %>%
    filter(var_value %in% c(min_dmi,max_dmi)) %>%
    mutate(var_value_char=factor(as.character(var_value),levels=c(as.character(c(min_dmi,max_dmi)))))

  plot_dmi <- results_dmi_filtered %>%
    ggplot(aes(x=lag_months,y=RR,col=var_value_char,fill=var_value_char)) +
    geom_ribbon(aes(ymin=RR_lower_95ci,ymax=RR_upper_95ci),alpha=0.15,col=NA) +
    geom_line() +
    geom_hline(yintercept=1,col="red",linetype=2) +
    theme_Publication(base_size=10) +
    labs(x="Lag (months)",y="Relative risk",col=NULL,fill=NULL) +
    theme(legend.position = "inside",legend.position.inside = c(0.5,0.8)) +
    scale_x_continuous(breaks=0:6) +
    scale_colour_discrete(palette=hcl.colors(2,"Blue-Red")) +
    scale_fill_discrete(palette=hcl.colors(2,"Blue-Red")) + theme(
      legend.background = element_rect(fill = "transparent", colour = NA), # Transparent legend box
      legend.key = element_rect(fill = "transparent", colour = NA),        # Transparent background for key symbols
      legend.box.background = element_rect(fill = "transparent", colour = NA) # Transparent legend panel
    )

  #### precipitation
  results_pre_filtered <- dlnm_output$results_pre %>%
    filter(var_value %in% c(min_pre,max_pre)) %>%
    mutate(var_value_char=factor(as.character(var_value),levels=c(as.character(c(min_pre,max_pre)))))

  plot_pre <- results_pre_filtered %>%
    ggplot(aes(x=lag_months,y=RR,col=var_value_char,fill=var_value_char)) +
    geom_ribbon(aes(ymin=RR_lower_95ci,ymax=RR_upper_95ci),alpha=0.15,col=NA) +
    geom_line() +
    geom_hline(yintercept=1,col="red",linetype=2) +
    theme_Publication(base_size=10) +
    labs(x="Lag (months)",y="Relative risk",col=NULL,fill=NULL) +
    theme(legend.position = "inside",legend.position.inside = c(0.5,0.8)) +
    scale_x_continuous(breaks=0:4) +
    scale_colour_discrete(palette=hcl.colors(2,"Blue-Red")) +
    scale_fill_discrete(palette=hcl.colors(2,"Blue-Red")) + theme(
      legend.background = element_rect(fill = "transparent", colour = NA), # Transparent legend box
      legend.key = element_rect(fill = "transparent", colour = NA),        # Transparent background for key symbols
      legend.box.background = element_rect(fill = "transparent", colour = NA) # Transparent legend panel
    )

  #### temperature
  results_tem_filtered <- dlnm_output$results_tem %>%
    filter(var_value %in% c(min_tem,max_tem)) %>%
    mutate(var_value_char=factor(as.character(var_value),levels=c(as.character(c(min_tem,max_tem)))))

  plot_tem <- results_tem_filtered %>%
    ggplot(aes(x=lag_months,y=RR,col=var_value_char,fill=var_value_char)) +
    geom_ribbon(aes(ymin=RR_lower_95ci,ymax=RR_upper_95ci),alpha=0.15,col=NA) +
    geom_line() +
    geom_hline(yintercept=1,col="red",linetype=2) +
    theme_Publication(base_size=10) +
    labs(x="Lag (months)",y="Relative risk",col=NULL,fill=NULL) +
    theme(legend.position = "inside",legend.position.inside = c(0.5,0.8)) +
    scale_x_continuous(breaks=0:4) +
    scale_colour_discrete(palette=hcl.colors(2,"Blue-Red")) +
    scale_fill_discrete(palette=hcl.colors(2,"Blue-Red")) + theme(
      legend.background = element_rect(fill = "transparent", colour = NA), # Transparent legend box
      legend.key = element_rect(fill = "transparent", colour = NA),        # Transparent background for key symbols
      legend.box.background = element_rect(fill = "transparent", colour = NA) # Transparent legend panel
    )

  #### relative humidity
  results_hum_filtered <- dlnm_output$results_hum %>%
    filter(var_value %in% c(min_hum,max_hum)) %>%
    mutate(var_value_char=factor(as.character(var_value),levels=c(as.character(c(min_hum,max_hum)))))

  plot_hum <- results_hum_filtered %>%
    ggplot(aes(x=lag_months,y=RR,col=var_value_char,fill=var_value_char)) +
    geom_ribbon(aes(ymin=RR_lower_95ci,ymax=RR_upper_95ci),alpha=0.15,col=NA) +
    geom_line() +
    geom_hline(yintercept=1,col="red",linetype=2) +
    theme_Publication(base_size=10) +
    labs(x="Lag (months)",y="Relative risk",col=NULL,fill=NULL) +
    theme(legend.position = "inside",legend.position.inside = c(0.5,0.8)) +
    scale_x_continuous(breaks=0:4) +
    scale_colour_discrete(palette=hcl.colors(2,"Blue-Red")) +
    scale_fill_discrete(palette=hcl.colors(2,"Blue-Red")) + theme(
      legend.background = element_rect(fill = "transparent", colour = NA), # Transparent legend box
      legend.key = element_rect(fill = "transparent", colour = NA),        # Transparent background for key symbols
      legend.box.background = element_rect(fill = "transparent", colour = NA) # Transparent legend panel
    )

  return(list(plot_oni=plot_oni,
              plot_dmi=plot_dmi,
              plot_pre=plot_pre,
              plot_tem=plot_tem,
              plot_hum=plot_hum))

}

#' Generate contour plots for DLNM results
dlnm_contour_visualise <- function(dlnm_output,data_input){

  #### precipitation
  cen_pre = round(median(data_input$precipitation),-1) # to the nearest tens

  #### temperature
  cen_tem = round(median(data_input$temperature))

  #### relative humidity
  cen_hum = round(median(data_input$rel_humidity))

  rr_labels <- c(0.4,0.67,1,1.5,2.5)
  log_breaks <- log(rr_labels)

  contour_oni <- dlnm_output$results_oni %>%
    ggplot(aes(x=lag_months,y=var_value,z=log(RR))) +
    geom_raster(aes(fill=log(RR))) +
    geom_contour(col="black") +
    theme_Publication(base_size=10) +
    labs(x="Lag (months)",y="ONI",fill="RR") +
    geom_hline(yintercept=0,col="red",linetype=2) +
    scale_fill_gradient2(
      low = "#147833",
      mid = "white",
      high = "#582B7C",
      midpoint = 0,
      limits = c(log(0.4), log(2.5)),
      name = "RR",
      breaks = log_breaks,
      labels = rr_labels,
      oob = scales::squish
    )

  contour_dmi <- dlnm_output$results_dmi %>%
    ggplot(aes(x=lag_months,y=var_value,z=log(RR))) +
    geom_raster(aes(fill=log(RR))) +
    geom_contour(col="black") +
    theme_Publication(base_size=10) +
    labs(x="Lag (months)",y="DMI",fill="RR") +
    geom_hline(yintercept=0,col="red",linetype=2) +
    scale_fill_gradient2(
      low = "#147833",
      mid = "white",
      high = "#582B7C",
      midpoint = 0,
      limits = c(log(0.4), log(2.5)),
      name = "RR",
      breaks = log_breaks,
      labels = rr_labels,
      oob = scales::squish
    )

  contour_pre <- dlnm_output$results_pre %>%
    ggplot(aes(x=lag_months,y=var_value,z=log(RR))) +
    geom_raster(aes(fill=log(RR))) +
    geom_contour(col="black") +
    theme_Publication(base_size=10) +
    labs(x="Lag (months)",y="Precipitation",fill="RR") +
    geom_hline(yintercept=cen_pre,col="red",linetype=2) +
    scale_fill_gradient2(
      low = "#147833",
      mid = "white",
      high = "#582B7C",
      midpoint = 0,
      limits = c(log(0.4), log(2.5)),
      name = "RR",
      breaks = log_breaks,
      labels = rr_labels,
      oob = scales::squish
    )

  contour_tem <- dlnm_output$results_tem %>%
    ggplot(aes(x=lag_months,y=var_value,z=log(RR))) +
    geom_raster(aes(fill=log(RR))) +
    geom_contour(col="black") +
    theme_Publication(base_size=10) +
    labs(x="Lag (months)",y="Temperature",fill="RR") +
    geom_hline(yintercept=cen_tem,col="red",linetype=2) +
    scale_fill_gradient2(
      low = "#147833",
      mid = "white",
      high = "#582B7C",
      midpoint = 0,
      limits = c(log(0.4), log(2.5)),
      name = "RR",
      breaks = log_breaks,
      labels = rr_labels,
      oob = scales::squish
    )

  contour_hum <- dlnm_output$results_hum %>%
    ggplot(aes(x=lag_months,y=var_value,z=log(RR))) +
    geom_raster(aes(fill=log(RR))) +
    geom_contour(col="black") +
    theme_Publication(base_size=10) +
    labs(x="Lag (months)",y="Relative humidity",fill="RR") +
    geom_hline(yintercept=cen_hum,col="red",linetype=2) +
    scale_fill_gradient2(
      low = "#147833",
      mid = "white",
      high = "#582B7C",
      midpoint = 0,
      limits = c(log(0.4), log(2.5)),
      name = "RR",
      breaks = log_breaks,
      labels = rr_labels,
      oob = scales::squish
    )

  return(list(contour_oni=contour_oni,
              contour_dmi=contour_dmi,
              contour_pre=contour_pre,
              contour_tem=contour_tem,
              contour_hum=contour_hum))

}

#' Run DLNM for Multiple Climate Variables
#'
#' Fits separate distributed lag non-linear models for each climate variable.
#' Each model is univariate, adjusting only for monthly seasonality.
#'
#' @param data_input Data frame with columns: cases, pop, month, oni, dmi,
#'   precipitation, temperature, rel_humidity
#' @param lag_value_long Maximum lag for ONI/DMI (typically 6 months)
#' @param lag_value_short Maximum lag for local climate (typically 3 months)
#'
#' @return List containing fitted models, predictions, and results tables
#'
#' @note Models are fitted separately for each climate variable. This approach
#'   does not account for correlations between climate variables but allows
#'   clear interpretation of individual dose-response relationships.
dlnm_run <- function(data_input, lag_value_long=6, lag_value_short=4) {

  data_input <- data_input %>%
    mutate(cases = na_if(cases, 0))

  # Filter data for provinces 91 and 94 (start from 2015)
  if(unique(data_input$idadmin1) %in% c(91,94)){
    data_input <- data_input %>%
      filter(year(date) >= 2015)
  }

  # Create time_index for ALL provinces (needed in model formula)
  data_input <- data_input %>%
    mutate(time_index = 1:nrow(data_input))
  
  #### Create DLNM crossbasis object for GLM
  #### ONI
  cb_oni <- crossbasis(
    x = data_input$oni,
    lag = lag_value_long,
    argvar = list(fun = "ns", df = 3),
    arglag = list(fun = "ns", df = 4)
  )
  
  #### IOD/DMI
  cb_dmi <- crossbasis(
    x = data_input$dmi,
    lag = lag_value_long,
    argvar = list(fun = "ns", df = 3),
    arglag = list(fun = "ns", df = 4)
  )
  
  #### precipitation
  cb_pre <- crossbasis(
    x = data_input$precipitation,
    lag = lag_value_short,
    argvar = list(fun = "ns", df = 4),
    arglag = list(fun = "ns", df = 3)
  )
  
  #### temperature
  cb_tem <- crossbasis(
    x = data_input$temperature,
    lag = lag_value_short,
    argvar = list(fun = "ns", df = 4),
    arglag = list(fun = "ns", df = 3)
  )
  
  ##### relative humidity
  cb_hum <- crossbasis(
    x = data_input$rel_humidity,
    lag = lag_value_short,
    argvar = list(fun = "ns", df = 4),
    arglag = list(fun = "ns", df = 3)
  )
  
  year_knots <- 15 # number of year knots for time_index spline
  yearly_multiplier <- 1 # knots for time_index spline
  #### GLM using DLNM crossbasis object - time hardcoded for 15 years
  if(unique(data_input$idadmin1) %in% c(91,94)){
    model_glm_oni <- glm.nb(cases ~ cb_oni + ns(time_index, df = yearly_multiplier*year_knots*2/3) + factor(month) + offset(log(pop)), 
                            link = "log", data = data_input)
    model_glm_dmi <- glm.nb(cases ~ cb_dmi + ns(time_index, df = yearly_multiplier*year_knots*2/3) + factor(month) + offset(log(pop)), 
                            link = "log", data = data_input)
    model_glm_tem <- glm.nb(cases ~ cb_tem + ns(time_index, df = yearly_multiplier*year_knots*2/3) + factor(month) + offset(log(pop)), 
                            link = "log", data = data_input)
    model_glm_pre <- glm.nb(cases ~ cb_pre + ns(time_index, df = yearly_multiplier*year_knots*2/3) + factor(month) + offset(log(pop)), 
                            link = "log", data = data_input)
    model_glm_hum <- glm.nb(cases ~ cb_hum + ns(time_index, df = yearly_multiplier*year_knots*2/3) + factor(month) + offset(log(pop)), 
                            link = "log", data = data_input)
  } else{
    model_glm_oni <- glm.nb(cases ~ cb_oni + ns(time_index, df = yearly_multiplier*year_knots) + factor(month) + offset(log(pop)), 
                            link = "log", data = data_input)
    model_glm_dmi <- glm.nb(cases ~ cb_dmi + ns(time_index, df = yearly_multiplier*year_knots) + factor(month) + offset(log(pop)), 
                            link = "log", data = data_input)
    model_glm_tem <- glm.nb(cases ~ cb_tem + ns(time_index, df = yearly_multiplier*year_knots) + factor(month) + offset(log(pop)), 
                            link = "log", data = data_input)
    model_glm_pre <- glm.nb(cases ~ cb_pre + ns(time_index, df = yearly_multiplier*year_knots) + factor(month) + offset(log(pop)), 
                            link = "log", data = data_input)
    model_glm_hum <- glm.nb(cases ~ cb_hum + ns(time_index, df = yearly_multiplier*year_knots) + factor(month) + offset(log(pop)), 
                            link = "log", data = data_input)
  }
  
  #### Value range for predictions
  #### ONI
  cen_oni = 0
  min_oni = -1
  max_oni = 2
  
  #### DMI
  cen_dmi = 0
  min_dmi = -0.5
  max_dmi = 1
  
  #### precipitation
  cen_pre = round(median(data_input$precipitation),-1) # to the nearest tens
  min_pre = as.numeric(round(quantile(data_input$precipitation,probs=c(0.10)),-1))
  max_pre = as.numeric(round(quantile(data_input$precipitation,probs=c(0.90)),-1))
  
  #### temperature
  cen_tem = round(median(data_input$temperature))
  min_tem = as.numeric(floor(quantile(data_input$temperature,probs=c(0.10))))
  max_tem = as.numeric(ceiling(quantile(data_input$temperature,probs=c(0.90))))
  
  #### relative humidity
  cen_hum = round(median(data_input$rel_humidity))
  min_hum = as.numeric(floor(quantile(data_input$rel_humidity,probs=c(0.10))))
  max_hum = as.numeric(ceiling(quantile(data_input$rel_humidity,probs=c(0.90))))
  
  #### Produce prediction matrices for contour
  pred_glm_oni <- crosspred(
    basis = cb_oni,
    model = model_glm_oni,
    cen = cen_oni,
    at = seq(min_oni, max_oni, length=300),
    bylag = 0.1,
    cumul = TRUE
  )
  
  pred_glm_dmi <- crosspred(
    basis = cb_dmi,
    model = model_glm_dmi,
    cen = cen_dmi,
    at = seq(min_dmi, max_dmi, length=300),
    bylag = 0.1,
    cumul = TRUE
  )
  
  pred_glm_pre <- crosspred(
    basis = cb_pre,
    model = model_glm_pre,
    cen = cen_pre,
    at = seq(min_pre, max_pre, length=300),
    bylag = 0.1,
    cumul = TRUE
  )
  
  pred_glm_tem <- crosspred(
    basis = cb_tem,
    model = model_glm_tem,
    cen = cen_tem,
    at = seq(min_tem, max_tem, length=300),
    bylag = 0.1,
    cumul = TRUE
  )
  
  pred_glm_hum <- crosspred(
    basis = cb_hum,
    model = model_glm_hum,
    cen = cen_hum,
    at = seq(min_hum, max_hum, length=300),
    bylag = 0.1,
    cumul = TRUE
  )
  
  idadmin1 <- unique(data_input$idadmin1)
  
  #### Convert predictions into tables
  results_oni <- data.frame(
    idadmin1 = idadmin1,
    clim_var = "ONI",
    var_value = rep(seq(min_oni, max_oni, length=300), times = length(seq(0,lag_value_long,by=0.1))),
    lag_months = rep(seq(0,lag_value_long,by=0.1), each = 300),
    RR = as.vector(pred_glm_oni$matRRfit),
    RR_lower_95ci = as.vector(pred_glm_oni$matRRlow),
    RR_upper_95ci = as.vector(pred_glm_oni$matRRhigh)
  )
  
  results_dmi <- data.frame(
    idadmin1 = idadmin1,
    clim_var = "DMI",
    var_value = rep(seq(min_dmi, max_dmi, length=300), times = length(seq(0,lag_value_long,by=0.1))),
    lag_months = rep(seq(0,lag_value_long,by=0.1), each = 300),
    RR = as.vector(pred_glm_dmi$matRRfit),
    RR_lower_95ci = as.vector(pred_glm_dmi$matRRlow),
    RR_upper_95ci = as.vector(pred_glm_dmi$matRRhigh)
  )
  
  results_pre <- data.frame(
    idadmin1 = idadmin1,
    clim_var = "Precipitation",
    var_value = rep(seq(min_pre, max_pre, length=300), times = length(seq(0,lag_value_short,by=0.1))),
    lag_months = rep(seq(0,lag_value_short,by=0.1), each = 300),
    RR = as.vector(pred_glm_pre$matRRfit),
    RR_lower_95ci = as.vector(pred_glm_pre$matRRlow),
    RR_upper_95ci = as.vector(pred_glm_pre$matRRhigh)
  )
  
  results_tem <- data.frame(
    idadmin1 = idadmin1,
    clim_var = "Temperature",
    var_value = rep(seq(min_tem, max_tem, length=300), times = length(seq(0,lag_value_short,by=0.1))),
    lag_months = rep(seq(0,lag_value_short,by=0.1), each = 300),
    RR = as.vector(pred_glm_tem$matRRfit),
    RR_lower_95ci = as.vector(pred_glm_tem$matRRlow),
    RR_upper_95ci = as.vector(pred_glm_tem$matRRhigh)
  )
  
  results_hum <- data.frame(
    idadmin1 = idadmin1,
    clim_var = "Relative humidity",
    var_value = rep(seq(min_hum, max_hum, length=300), times = length(seq(0,lag_value_short,by=0.1))),
    lag_months = rep(seq(0,lag_value_short,by=0.1), each = 300),
    RR = as.vector(pred_glm_hum$matRRfit),
    RR_lower_95ci = as.vector(pred_glm_hum$matRRlow),
    RR_upper_95ci = as.vector(pred_glm_hum$matRRhigh)
  )
  
  cumulative_RR <- tibble()
  # regional variable
  for (t in 0:6){
    
    cumulative_RR <- bind_rows(cumulative_RR,
                               tibble(idadmin1=idadmin1,
                                      clim_var="ONI",
                                      percentile=10,
                                      var_value=min_oni,
                                      ref_value=cen_oni,
                                      lag=t,
                                      cumulRR=pred_glm_oni$cumRRfit[1, t+1],
                                      cumulRR_lo=pred_glm_oni$cumRRlow[1, t+1],
                                      cumulRR_hi=pred_glm_oni$cumRRhigh[1, t+1]))
    
    cumulative_RR <- bind_rows(cumulative_RR,
                               tibble(idadmin1=idadmin1,
                                      clim_var="ONI",
                                      percentile=90,
                                      var_value=max_oni,
                                      ref_value=cen_oni,
                                      lag=t,
                                      cumulRR=pred_glm_oni$cumRRfit[300, t+1],
                                      cumulRR_lo=pred_glm_oni$cumRRlow[300, t+1],
                                      cumulRR_hi=pred_glm_oni$cumRRhigh[300, t+1]))
    
    cumulative_RR <- bind_rows(cumulative_RR,
                               tibble(idadmin1=idadmin1,
                                      clim_var="DMI",
                                      percentile=10,
                                      var_value=min_dmi,
                                      ref_value=cen_dmi,
                                      lag=t,
                                      cumulRR=pred_glm_dmi$cumRRfit[1, t+1],
                                      cumulRR_lo=pred_glm_dmi$cumRRlow[1, t+1],
                                      cumulRR_hi=pred_glm_dmi$cumRRhigh[1, t+1]))
    
    cumulative_RR <- bind_rows(cumulative_RR,
                               tibble(idadmin1=idadmin1,
                                      clim_var="DMI",
                                      percentile=90,
                                      var_value=max_dmi,
                                      ref_value=cen_dmi,
                                      lag=t,
                                      cumulRR=pred_glm_dmi$cumRRfit[300, t+1],
                                      cumulRR_lo=pred_glm_dmi$cumRRlow[300, t+1],
                                      cumulRR_hi=pred_glm_dmi$cumRRhigh[300, t+1]))
    
  }
  
  for (t in 0:4){
    
    cumulative_RR <- bind_rows(cumulative_RR,
                               tibble(idadmin1=idadmin1,
                                      clim_var="Precipitation",
                                      percentile=10,
                                      var_value=min_pre,
                                      ref_value=cen_pre,
                                      lag=t,
                                      cumulRR=pred_glm_pre$cumRRfit[1, t+1],
                                      cumulRR_lo=pred_glm_pre$cumRRlow[1, t+1],
                                      cumulRR_hi=pred_glm_pre$cumRRhigh[1, t+1]))
    
    cumulative_RR <- bind_rows(cumulative_RR,
                               tibble(idadmin1=idadmin1,
                                      clim_var="Precipitation",
                                      percentile=90,
                                      var_value=max_pre,
                                      ref_value=cen_pre,
                                      lag=t,
                                      cumulRR=pred_glm_pre$cumRRfit[300, t+1],
                                      cumulRR_lo=pred_glm_pre$cumRRlow[300, t+1],
                                      cumulRR_hi=pred_glm_pre$cumRRhigh[300, t+1]))
    
    cumulative_RR <- bind_rows(cumulative_RR,
                               tibble(idadmin1=idadmin1,
                                      clim_var="Temperature",
                                      percentile=10,
                                      var_value=min_tem,
                                      ref_value=cen_tem,
                                      lag=t,
                                      cumulRR=pred_glm_tem$cumRRfit[1, t+1],
                                      cumulRR_lo=pred_glm_tem$cumRRlow[1, t+1],
                                      cumulRR_hi=pred_glm_tem$cumRRhigh[1, t+1]))
    
    cumulative_RR <- bind_rows(cumulative_RR,
                               tibble(idadmin1=idadmin1,
                                      clim_var="Temperature",
                                      percentile=90,
                                      var_value=max_tem,
                                      ref_value=cen_tem,
                                      lag=t,
                                      cumulRR=pred_glm_tem$cumRRfit[300, t+1],
                                      cumulRR_lo=pred_glm_tem$cumRRlow[300, t+1],
                                      cumulRR_hi=pred_glm_tem$cumRRhigh[300, t+1]))
    
    cumulative_RR <- bind_rows(cumulative_RR,
                               tibble(idadmin1=idadmin1,
                                      clim_var="Relative humidity",
                                      percentile=10,
                                      var_value=min_hum,
                                      ref_value=cen_hum,
                                      lag=t,
                                      cumulRR=pred_glm_hum$cumRRfit[1, t+1],
                                      cumulRR_lo=pred_glm_hum$cumRRlow[1, t+1],
                                      cumulRR_hi=pred_glm_hum$cumRRhigh[1, t+1]))
    
    cumulative_RR <- bind_rows(cumulative_RR,
                               tibble(idadmin1=idadmin1,
                                      clim_var="Relative humidity",
                                      percentile=90,
                                      var_value=max_hum,
                                      ref_value=cen_hum,
                                      lag=t,
                                      cumulRR=pred_glm_hum$cumRRfit[300, t+1],
                                      cumulRR_lo=pred_glm_hum$cumRRlow[300, t+1],
                                      cumulRR_hi=pred_glm_hum$cumRRhigh[300, t+1]))
    
  }
  
  allCumulative_RR <- tibble()
  
  allCumulative_RR <- bind_rows(allCumulative_RR,
                                tibble(idadmin1=idadmin1,
                                       clim_var="ONI",
                                       percentile=10,
                                       var_value=min_oni,
                                       ref_value=cen_oni,
                                       lag=6,
                                       cumulRR=pred_glm_oni$allRRfit[1],
                                       cumulRR_lo=pred_glm_oni$allRRlow[1],
                                       cumulRR_hi=pred_glm_oni$allRRhigh[1]))
  
  allCumulative_RR <- bind_rows(allCumulative_RR,
                                tibble(idadmin1=idadmin1,
                                       clim_var="ONI",
                                       percentile=90,
                                       var_value=max_oni,
                                       ref_value=cen_oni,
                                       lag=6,
                                       cumulRR=pred_glm_oni$allRRfit[300],
                                       cumulRR_lo=pred_glm_oni$allRRlow[300],
                                       cumulRR_hi=pred_glm_oni$allRRhigh[300]))
  
  allCumulative_RR <- bind_rows(allCumulative_RR,
                                tibble(idadmin1=idadmin1,
                                       clim_var="DMI",
                                       percentile=10,
                                       var_value=min_dmi,
                                       ref_value=cen_dmi,
                                       lag=6,
                                       cumulRR=pred_glm_dmi$allRRfit[1],
                                       cumulRR_lo=pred_glm_dmi$allRRlow[1],
                                       cumulRR_hi=pred_glm_dmi$allRRhigh[1]))
  
  allCumulative_RR <- bind_rows(allCumulative_RR,
                                tibble(idadmin1=idadmin1,
                                       clim_var="DMI",
                                       percentile=90,
                                       var_value=max_dmi,
                                       ref_value=cen_dmi,
                                       lag=6,
                                       cumulRR=pred_glm_dmi$allRRfit[300],
                                       cumulRR_lo=pred_glm_dmi$allRRlow[300],
                                       cumulRR_hi=pred_glm_dmi$allRRhigh[300]))
  
  allCumulative_RR <- bind_rows(allCumulative_RR,
                                tibble(idadmin1=idadmin1,
                                       clim_var="Precipitation",
                                       percentile=10,
                                       var_value=min_pre,
                                       ref_value=cen_pre,
                                       lag=4,
                                       cumulRR=pred_glm_pre$allRRfit[1],
                                       cumulRR_lo=pred_glm_pre$allRRlow[1],
                                       cumulRR_hi=pred_glm_pre$allRRhigh[1]))
  
  allCumulative_RR <- bind_rows(allCumulative_RR,
                                tibble(idadmin1=idadmin1,
                                       clim_var="Precipitation",
                                       percentile=90,
                                       var_value=max_pre,
                                       ref_value=cen_pre,
                                       lag=4,
                                       cumulRR=pred_glm_pre$allRRfit[300],
                                       cumulRR_lo=pred_glm_pre$allRRlow[300],
                                       cumulRR_hi=pred_glm_pre$allRRhigh[300]))
  
  allCumulative_RR <- bind_rows(allCumulative_RR,
                                tibble(idadmin1=idadmin1,
                                       clim_var="Temperature",
                                       percentile=10,
                                       var_value=min_tem,
                                       ref_value=cen_tem,
                                       lag=4,
                                       cumulRR=pred_glm_tem$allRRfit[1],
                                       cumulRR_lo=pred_glm_tem$allRRlow[1],
                                       cumulRR_hi=pred_glm_tem$allRRhigh[1]))
  
  allCumulative_RR <- bind_rows(allCumulative_RR,
                                tibble(idadmin1=idadmin1,
                                       clim_var="Temperature",
                                       percentile=90,
                                       var_value=max_tem,
                                       ref_value=cen_tem,
                                       lag=4,
                                       cumulRR=pred_glm_tem$allRRfit[300],
                                       cumulRR_lo=pred_glm_tem$allRRlow[300],
                                       cumulRR_hi=pred_glm_tem$allRRhigh[300]))
  
  allCumulative_RR <- bind_rows(allCumulative_RR,
                                tibble(idadmin1=idadmin1,
                                       clim_var="Relative humidity",
                                       percentile=10,
                                       var_value=min_hum,
                                       ref_value=cen_hum,
                                       lag=4,
                                       cumulRR=pred_glm_hum$allRRfit[1],
                                       cumulRR_lo=pred_glm_hum$allRRlow[1],
                                       cumulRR_hi=pred_glm_hum$allRRhigh[1]))
  
  allCumulative_RR <- bind_rows(allCumulative_RR,
                                tibble(idadmin1=idadmin1,
                                       clim_var="Relative humidity",
                                       percentile=90,
                                       var_value=max_hum,
                                       ref_value=cen_hum,
                                       lag=4,
                                       cumulRR=pred_glm_hum$allRRfit[300],
                                       cumulRR_lo=pred_glm_hum$allRRlow[300],
                                       cumulRR_hi=pred_glm_hum$allRRhigh[300]))
  
  #### Output
  return(list(model_oni=model_glm_oni,
              model_dmi=model_glm_dmi,
              model_pre=model_glm_pre,
              model_tem=model_glm_tem,
              model_hum=model_glm_hum,
              AIC_oni=AIC(model_glm_oni),
              BIC_oni=BIC(model_glm_oni),
              deviance_oni=deviance(model_glm_oni),
              AIC_dmi=AIC(model_glm_dmi),
              BIC_dmi=BIC(model_glm_dmi),
              deviance_dmi=deviance(model_glm_dmi),
              AIC_pre=AIC(model_glm_pre),
              BIC_pre=BIC(model_glm_pre),
              deviance_pre=deviance(model_glm_pre),
              AIC_tem=AIC(model_glm_tem),
              BIC_tem=BIC(model_glm_tem),
              deviance_tem=deviance(model_glm_tem),
              AIC_hum=AIC(model_glm_hum),
              BIC_hum=BIC(model_glm_hum),
              deviance_hum=deviance(model_glm_hum),
              pred_glm_oni=pred_glm_oni,
              pred_glm_dmi=pred_glm_dmi,
              pred_glm_pre=pred_glm_pre,
              pred_glm_tem=pred_glm_tem,
              pred_glm_hum=pred_glm_hum,
              results_oni=results_oni,
              results_dmi=results_dmi,
              results_pre=results_pre,
              results_tem=results_tem,
              results_hum=results_hum,
              cb_oni=cb_oni,
              cb_dmi=cb_dmi,
              cb_pre=cb_pre,
              cb_tem=cb_tem,
              cb_hum=cb_hum,
              cumulative_RR=cumulative_RR,
              allCumulative_RR=allCumulative_RR))
  
}

# Function to find optimal lag for DLNM results based on statistical significance and cumulRR
# For elevated risk: prioritizes significant elevated risk (>1), then any significant, then highest overall
find_optimal_lag <- function(data) {
  # Add significance flag: TRUE if range doesn't contain 1
  data$significance <- !(data$cumulRR_lo <= 1 & data$cumulRR_hi >= 1)
  data$sig_label <- ifelse(data$significance, "s.s.", "n.s.")
  
  # Separate significant and non-significant
  significant <- data[data$significance == TRUE, ]
  non_significant <- data[data$significance == FALSE, ]
  
  # Further filter significant results: only keep those with cumulRR > 1 (elevated risk)
  significant_elevated <- significant[significant$cumulRR > 1, ]
  
  # Determine which lag to return
  if (nrow(significant_elevated) > 0) {
    # If significant elevated results exist, pick the one with highest cumulRR
    optimal_idx <- which.max(significant_elevated$cumulRR)
    result <- significant_elevated[optimal_idx, ]
    result$selection_criterion <- "Highest elevated significant cumulRR (>1)"
  } else {
    # If no significant results, pick highest cumulRR overall
    optimal_idx <- which.max(data$cumulRR)
    result <- data[optimal_idx, ]
    result$selection_criterion <- "Highest cumulRR (no significant found)"
  }
  
  return(result)
}

# Function to find FIRST lag with significant elevated risk (DLNM results)
# Returns the first instance where cumulRR is significantly elevated (>1)
find_optimal_lag_first <- function(data) {
  # Ensure data is sorted by lag
  data <- data[order(data$lag), ]
  
  # Add significance flag: TRUE if range doesn't contain 1
  data$significance <- !(data$cumulRR_lo <= 1 & data$cumulRR_hi >= 1)
  data$sig_label <- ifelse(data$significance, "s.s.", "n.s.")
  
  # Filter for significant AND elevated risk (cumulRR > 1)
  significant_elevated <- data[data$significance == TRUE & data$cumulRR > 1, ]
  
  # Determine which lag to return
  if (nrow(significant_elevated) > 0) {
    # If significant elevated results exist, pick the FIRST one (first lag)
    result <- significant_elevated[1, ]
    result$selection_criterion <- "First lag with significant elevated cumulRR (>1)"
  } else {
    # If no significant elevated results, pick highest cumulRR overall
    optimal_idx <- which.max(data$cumulRR)
    result <- data[optimal_idx, ]
    result$selection_criterion <- "Highest cumulRR (no significant elevated found)"
  }
  
  return(result)
}

################################################################################
# 6. PEAK MONTH DETECTION FUNCTIONS
################################################################################

#' Find peak months from wavelet analysis
#'
#' Identifies peak outbreak months for each epidemic year (July-June) based on
#' wavelet phase angles. Peaks are identified as the time point with phase angle
#' closest to zero within each period.
#'
#' @param phase_angles Numeric vector of phase angle values from wavelet analysis
#' @param start_year Starting year (default: 2010)
#' @param end_year Ending year (default: 2024)
#' @param anomaly_threshold Threshold in radians for flagging anomalous periods
#'   (default: 0.3 radians, approximately 0.57 months). Periods where the minimum
#'   absolute phase angle exceeds this threshold are flagged as anomalous,
#'   indicating weak or unclear annual periodicity. This value was chosen based
#'   on visual inspection of phase angle distributions; sensitivity analyses
#'   with thresholds of 0.2-0.4 radians showed consistent results.
#' @return A list containing peak indices, months, years, and average peak month
find_peak_months <- function(phase_angles, start_year = 2010, end_year = 2024,
                             anomaly_threshold = 0.3) {

  # Validate input
  expected_length <- (end_year - start_year + 1) * 12
  if (length(phase_angles) != expected_length) {
    stop(paste0("Phase angles vector must have length ", expected_length,
                " (12 months x ", end_year - start_year + 1, " years)"))
  }

  # Period year June-July
  period_year <- fiscal_year_vector(start_year, end_year)

  # Remove first 6 and last 6 months (incomplete periods)
  valid_indices <- 7:(length(phase_angles) - 6)
  valid_phase_angles <- phase_angles[valid_indices]
  valid_period_year <- period_year[valid_indices]

  # Calculate absolute values to find values closest to 0
  abs_angles <- abs(valid_phase_angles)

  # Get unique periods
  unique_periods <- unique(valid_period_year)
  n_periods <- length(unique_periods)

  # Find one peak per period
  peak_indices_in_valid <- integer(n_periods)
  peak_values <- numeric(n_periods)
  peak_periods <- character(n_periods)
  peak_abs_values <- numeric(n_periods)

  for (i in seq_along(unique_periods)) {
    period <- unique_periods[i]

    # Get indices for this period
    period_mask <- valid_period_year == period
    period_indices <- which(period_mask)

    # Get absolute angles for this period
    period_abs_angles <- abs_angles[period_mask]

    # Find the minimum absolute value
    min_idx_within_period <- which.min(period_abs_angles)

    # Store the index (relative to valid_indices)
    peak_indices_in_valid[i] <- period_indices[min_idx_within_period]
    peak_values[i] <- valid_phase_angles[period_indices[min_idx_within_period]]
    peak_abs_values[i] <- period_abs_angles[min_idx_within_period]
    peak_periods[i] <- period
  }

  # Flag anomalous periods where abs(phase_angle) > threshold
  is_anomalous <- peak_abs_values > anomaly_threshold

  # Convert back to original indices
  peak_indices <- valid_indices[peak_indices_in_valid]

  # Extract months and years
  peak_months_jan_dec <- ((peak_indices - 1) %% 12) + 1
  peak_years_original <- start_year + ((peak_indices - 1) %/% 12)

  # Set anomalous peaks to NA
  peak_months_jan_dec[is_anomalous] <- NA
  peak_years_original[is_anomalous] <- NA

  # Convert months from Jan-Dec (1-12) to Jul-Jun (1-12) scale
  peak_months_jul_jun <- ifelse(is.na(peak_months_jan_dec), NA,
                                ifelse(peak_months_jan_dec >= 7,
                                       peak_months_jan_dec - 6,
                                       peak_months_jan_dec + 6))

  # Adjust for delayed peaks
  peak_months_jul_jun_adjusted <- peak_months_jul_jun

  for (i in seq_along(peak_months_jul_jun)) {
    if (is_anomalous[i]) {
      next
    }

    period_start_year <- as.integer(substr(peak_periods[i], 1, 4))
    is_jul_dec <- peak_months_jan_dec[i] >= 7

    if (is_jul_dec) {
      if (peak_years_original[i] > period_start_year) {
        peak_months_jul_jun_adjusted[i] <- peak_months_jul_jun[i] + 12
      }
    }
  }

  # Calculate average peak month using circular mean on Jul-Jun scale
  valid_for_avg <- !is_anomalous

  if (sum(valid_for_avg) > 0) {
    month_radians <- (peak_months_jul_jun[valid_for_avg] - 1) * (2 * pi / 12)
    avg_sin <- mean(sin(month_radians))
    avg_cos <- mean(cos(month_radians))
    avg_month_radians <- atan2(avg_sin, avg_cos)

    if (avg_month_radians < 0) {
      avg_month_radians <- avg_month_radians + 2 * pi
    }
    avg_month_jul_jun <- (avg_month_radians * 12 / (2 * pi)) + 1

    if (avg_month_jul_jun > 12) avg_month_jul_jun <- avg_month_jul_jun - 12
    if (avg_month_jul_jun < 1) avg_month_jul_jun <- avg_month_jul_jun + 12

    avg_month_jan_dec <- ifelse(avg_month_jul_jun <= 6,
                                avg_month_jul_jun + 6,
                                avg_month_jul_jun - 6)

    avg_month_rounded <- round(avg_month_jan_dec)
    if (avg_month_rounded < 1) avg_month_rounded <- 1
    if (avg_month_rounded > 12) avg_month_rounded <- 12
    avg_month_name <- month.name[avg_month_rounded]
  } else {
    avg_month_jan_dec <- NA
    avg_month_jul_jun <- NA
    avg_month_name <- NA
  }

  month_names <- ifelse(is.na(peak_months_jan_dec), NA,
                        month.name[peak_months_jan_dec])

  peak_data <- data.frame(
    Index = ifelse(is_anomalous, NA, peak_indices),
    Month_JanDec = peak_months_jan_dec,
    Month_JulJun = peak_months_jul_jun,
    Month_JulJun_Adjusted = peak_months_jul_jun_adjusted,
    Month_Name = month_names,
    Period = peak_periods,
    Year_original = peak_years_original,
    Period_Start_Year = as.integer(substr(peak_periods, 1, 4)),
    Phase_Angle = peak_values,
    Abs_Phase_Angle = peak_abs_values,
    Is_Anomalous = is_anomalous,
    Is_Delayed = ifelse(is_anomalous, NA, peak_months_jul_jun_adjusted > 12)
  )

  results <- list(
    peak_data = peak_data,
    average_peak_month_jan_dec = avg_month_jan_dec,
    average_peak_month_jul_jun = avg_month_jul_jun,
    average_peak_month_name = avg_month_name,
    n_periods = n_periods,
    n_anomalous_periods = sum(is_anomalous),
    n_delayed_peaks = sum(peak_months_jul_jun_adjusted > 12, na.rm = TRUE),
    anomaly_threshold = anomaly_threshold
  )

  return(results)
}

#' Generate July-June fiscal year labels
fiscal_year_vector <- function(start_year, end_year) {
  n_months <- (end_year - start_year + 1) * 12
  result <- character(n_months)
  idx <- 1

  for (year in start_year:end_year) {
    for (month in 1:6) {
      result[idx] <- paste0(year - 1, "-", year)
      idx <- idx + 1
    }
    for (month in 7:12) {
      result[idx] <- paste0(year, "-", year + 1)
      idx <- idx + 1
    }
  }

  return(result)
}

################################################################################
# 7. DTW CLIMATE-LAG FUNCTION
################################################################################

#' Find DTW similarity between cases and climate
dtw_climate_lag_enhanced <- function(cases, climate, max_lag=6) {

  cases_std <- as.vector(scale(cases))
  climate_std <- as.vector(scale(climate))

  # Test POSITIVE relationship
  alignment_pos <- dtw(climate_std, cases_std,
                       keep.internals=TRUE,
                       step.pattern=symmetric2,
                       window.type="slantedband",
                       window.size=max_lag)

  # Test NEGATIVE relationship
  climate_std_neg <- -climate_std

  alignment_neg <- dtw(climate_std_neg, cases_std,
                       keep.internals=TRUE,
                       step.pattern=symmetric2,
                       window.type="slantedband",
                       window.size=max_lag)

  if(alignment_pos$normalizedDistance < alignment_neg$normalizedDistance) {
    best_alignment <- alignment_pos
    relationship <- "positive"
    similarity_sign <- 1
  } else {
    best_alignment <- alignment_neg
    relationship <- "negative"
    similarity_sign <- -1
  }

  alignment_path <- data.frame(
    climate_index = best_alignment$index1,
    cases_index = best_alignment$index2
  )

  alignment_path$lag <- alignment_path$cases_index - alignment_path$climate_index
  modal_lag <- as.numeric(names(sort(table(alignment_path$lag), decreasing=TRUE)[1]))

  dtw_similarity <- 1 / (1 + best_alignment$normalizedDistance)

  if(modal_lag >= 0) {
    cases_subset <- cases[(modal_lag + 1):length(cases)]
    climate_subset <- climate[1:(length(climate) - modal_lag)]
  } else {
    cases_subset <- cases[1:(length(cases) + modal_lag)]
    climate_subset <- climate[(-modal_lag + 1):length(climate)]
  }

  pearson_cor <- cor(climate_subset, cases_subset, use="complete.obs")

  return(list(
    dtw_distance = best_alignment$distance,
    normalized_distance = best_alignment$normalizedDistance,
    dtw_similarity = dtw_similarity,
    optimal_lag = modal_lag,
    relationship = relationship,
    pearson_at_lag = pearson_cor,
    similarity_sign = similarity_sign,
    alignment = best_alignment,
    alignment_path = alignment_path
  ))
}

################################################################################
# 8. CONSTANTS
################################################################################

month_abb_juljun <- c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
                      "Jan", "Feb", "Mar", "Apr", "May", "Jun")

month_name_juljun <- c("July", "August", "September", "October", "November", "December",
                       "January", "February", "March", "April", "May", "June")

################################################################################
# 9. SENSITIVITY ANALYSIS
################################################################################

#' Test Phase Coherence Threshold Sensitivity
#'
#' Tests how the number of qualifying provinces changes across different
#' phase coherence thresholds
#'
#' @param phase_data Data frame with columns: province, climate_var,
#'   phase_lag_months, phase_coherence
#' @param dlnm_significance Data frame with DLNM significance results,
#'   containing columns: shapeName (province), clim_var, significant
#' @param thresholds Vector of threshold values to test (default: 0.70 to 0.90)
#'
#' @return Data frame with columns: threshold, climate_var, n_qualifying,
#'   n_significant_dlnm, pct_significant
#'
#' @examples
#' sensitivity <- test_coherence_sensitivity(
#'   phase_data = climate_dengue_phase_df,
#'   dlnm_significance = significance_positive,
#'   thresholds = c(0.70, 0.75, 0.80, 0.85, 0.90)
#' )
test_coherence_sensitivity <- function(phase_data,
                                       dlnm_significance = NULL,
                                       thresholds = seq(0.70, 0.90, by = 0.05)) {

  # Focus on precipitation and temperature only
  phase_subset <- phase_data %>%
    filter(climate_var %in% c("Precipitation", "Temperature"))

  # Initialize results list
  results_list <- list()

  # Test each threshold
  for (threshold in thresholds) {

    # For each climate variable
    for (clim_var in c("Precipitation", "Temperature")) {

      # Count qualifying provinces (high coherence + positive lag)
      qualifying <- phase_subset %>%
        filter(climate_var == clim_var,
               phase_coherence >= threshold,
               round(phase_lag_months) >= 0)

      n_qualifying <- nrow(qualifying)

      # If DLNM significance data provided, count provinces with significant effects
      n_significant <- NA
      pct_significant <- NA

      if (!is.null(dlnm_significance) && n_qualifying > 0) {

        # Match climate variable names between phase data and DLNM data
        dlnm_clim_var <- clim_var  # Should match format in dlnm_significance

        # Get significant provinces for this climate variable
        sig_provinces <- dlnm_significance %>%
          filter(clim_var == dlnm_clim_var) %>%
          pull(shapeName) %>%
          unique()

        # Count how many qualifying provinces have significant DLNM
        n_significant <- sum(qualifying$province %in% sig_provinces)
        pct_significant <- round(100 * n_significant / n_qualifying, 1)
      }

      # Store results
      results_list[[length(results_list) + 1]] <- tibble(
        threshold = threshold,
        climate_var = clim_var,
        n_qualifying = n_qualifying,
        n_significant_dlnm = n_significant,
        pct_significant = pct_significant
      )
    }
  }

  # Combine results
  results_df <- bind_rows(results_list)

  return(results_df)
}
