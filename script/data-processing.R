################################################################################
# Data Processing Script for Spatiotemporal Dengue-Climate Analysis
#
# This script processes all data needed for the manuscript figures:
#   1. Administrative boundaries and shapefiles
#   2. Dengue surveillance data (province-level, 2010-2024)
#   3. Climate data (ERA5-Land: precipitation, temperature, humidity)
#   4. Climate indices (ONI, DMI via rsoi package)
#   5. Wavelet analysis for phase lag estimation
#   6. DTW clustering for province grouping
#   7. Climate-dengue phase relationships
#   8. DLNM analysis
#
# Run this script first before running any figure scripts.
#
# Data sources:
#   - Dengue: Indonesian Ministry of Health and Open Dengue database
#   - Climate: ERA5-Land via Copernicus Climate Data Store
#   - ONI/DMI: NOAA via rsoi package
################################################################################

#### Set PROJ_LIB for spatial reference system
Sys.setenv(PROJ_LIB = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/sf/proj")

################################################################################
# SECTION 1: Load Required Libraries
################################################################################
# All packages needed for the entire analysis pipeline are loaded here.
# Figure scripts only need to source this file.

#### Core data manipulation
library(tidyverse)    # includes ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats
library(lubridate)    # date/time handling
library(reshape2)     # data reshaping
library(janitor)      # data cleaning
library(snakecase)    # case conversion
library(zoo)          # time series handling
library(readxl)       # Excel file reading

#### Spatial analysis
library(sf)           # spatial data handling
library(spdep)        # spatial autocorrelation (Moran's I)

#### Wavelet analysis
library(biwavelet)    # continuous wavelet transform

#### Time series clustering
library(dtwclust)     # DTW clustering
library(dtw)          # dynamic time warping

#### Statistical modeling
library(dlnm)         # distributed lag non-linear models
library(mgcv)         # generalized additive models
library(splines)      # spline functions
library(MASS)         # negative binomial regression
library(confintr)     # confidence intervals
library(circular)     # circular statistics
library(bpnreg)       # circular regression
library(ncf)          # spatial correlation functions
library(boot)         # bootstrap methods

#### Climate indices
library(rsoi)         # download ONI/DMI indices

#### Visualization
library(ggplot2)      # (included in tidyverse but explicit for clarity)
library(patchwork)    # combine plots
library(cowplot)      # plot arrangement
library(scales)       # scale functions
library(biscale)      # bivariate mapping
library(shadowtext)   # text with shadow
library(colorBlindness) # colorblind-friendly palettes
library(ggthemes)     # additional themes (for theme_foundation)
library(grid)         # grid graphics
library(viridis)      # viridis color scales
library(ggpubr)       # publication-ready plots (get_legend)

#### Tables and output
library(gt)           # publication-ready tables
library(writexl)      # Excel export

#### Load custom functions
source("script/functions.R")

#### Read and prepare administrative data
admin1_indonesia <- readRDS("data/admin1_indonesia.rds") %>% 
  mutate(admin1_no=seq_len(35))

admin2_indonesia <- readRDS("data/admin2_indonesia.rds") %>% 
  dplyr::select(-region) %>% 
  left_join(admin1_indonesia)

idadmin1_indonesia_unique <- unique(admin1_indonesia$idadmin1)
admin1_indonesia_unique <- unique(admin1_indonesia$admin1)

#### Read population data
pop_indonesia_admin2_2009_2024 <- readRDS("data/pop_age_simplegrowth_2009_2024.rds") %>% 
  bind_rows() %>% 
  group_by(idadmin2,year) %>% 
  summarise(pop=sum(pop)) %>% 
  ungroup()

pop_indonesia_admin1_2009_2024 <- pop_indonesia_admin2_2009_2024 %>% 
  left_join(admin2_indonesia) %>% 
  group_by(region,idadmin1,admin1,year) %>% 
  summarise(pop=sum(pop)) %>% 
  ungroup()

pop_indonesia_admin0_2009_2024 <- pop_indonesia_admin1_2009_2024 %>% 
  group_by(year) %>% 
  summarise(pop=sum(pop)) %>% 
  ungroup() %>% 
  mutate(region="NATIONAL",idadmin1=99,admin1="INDONESIA")

pop_indonesia_admin1 <- bind_rows(pop_indonesia_admin1_2009_2024 %>% 
                                    filter(year %in% 2010:2024),
                                  pop_indonesia_admin0_2009_2024 %>% 
                                    filter(year %in% 2010:2024)) %>% 
  mutate(admin1=ifelse(admin1=="SUMATERA UTARA","SUMATRA UTARA",admin1),
         admin1=ifelse(admin1=="SUMATERA BARAT","SUMATRA BARAT",admin1),
         admin1=ifelse(admin1=="SUMATERA SELATAN","SUMATRA SELATAN",admin1),
         region=ifelse(is.na(region),"SUMATRA",region))

#### Read shapefiles
admin1_shp <- read_sf(dsn="data/shapefiles/admin1", layer="admin1_34", stringsAsFactors = FALSE)
other_shp <- read_sf(dsn="data/shapefiles/other", layer="idn_neighbours", stringsAsFactors = FALSE)
region_shp <- read_sf(dsn="data/shapefiles/region", layer="region", stringsAsFactors = FALSE) %>% 
  mutate(region=to_title_case(region)) %>% 
  mutate(region=ifelse(region=="Sumatera","Sumatra",region)) %>% 
  mutate(region=ifelse(region=="Java Bali","Java & Bali",region))

#### Prepare administrative name mappings
admin1_EN <- admin1_shp %>% 
  st_drop_geometry() %>% 
  dplyr::select(shapeName, idadmin1) %>% 
  arrange(idadmin1)
admin1_EN$shapeName <- factor(admin1_EN$shapeName, levels = rev(admin1_EN$shapeName))

admin1_national_EN <- bind_rows(
  admin1_EN %>% mutate(shapeName = as.character(shapeName)),
  tibble(shapeName = "Indonesia", idadmin1 = 99)
) %>% 
  mutate(shapeName = factor(shapeName, levels = c("Indonesia", levels(admin1_EN$shapeName))))

admin1_indonesia_prov <- admin1_indonesia %>% filter(admin1 != "INDONESIA")

#### Process dengue surveillance data at admin1 level
dengue_data_admin1 <- readRDS("data/dengue_data_admin1_indonesia.rds") %>% 
  group_by(region, idadmin1, admin1, year, month) %>% 
  summarise(
    cases = sum(cases, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  left_join(pop_indonesia_admin1) %>% 
  mutate(
    incidence_rate = (cases / pop) * 100000,
    cases_add = ifelse(cases == 0, 0.01, cases),
    cases_NA = ifelse(cases == 0, NA, cases),
    incidence_rate_add = (cases_add / pop) * 100000,
    incidence_rate_NA = (cases_NA / pop) * 100000,
    month_name = factor(month, levels = 1:12, 
                        labels = c("January", "February", "March", "April", "May", "June", 
                                   "July", "August", "September", "October", "November", "December")),
    date = ymd(paste0(year, "-", month, "-01"))
  ) %>% 
  group_by(admin1) %>% 
  mutate(
    incidence_rate_log = log(incidence_rate_NA),
    incidence_rate_scaled = as.numeric(scale(incidence_rate_log))
  ) %>% 
  ungroup() %>% 
  arrange(idadmin1, date)

#### Region level
dengue_data_region <- dengue_data_admin1 %>% 
  group_by(region, year, month) %>% 
  summarise(
    cases = sum(cases, na.rm = TRUE),
    pop = sum(pop)
  ) %>% 
  ungroup() %>% 
  mutate(
    incidence_rate = (cases / pop) * 100000,
    cases_add = ifelse(cases == 0, 0.01, cases),
    incidence_rate_add = (cases_add / pop) * 100000,
    month_name = factor(month, levels = 1:12, 
                        labels = c("January", "February", "March", "April", "May", "June", 
                                   "July", "August", "September", "October", "November", "December")),
    date = ymd(paste0(year, "-", month, "-01"))
  ) %>% 
  group_by(region) %>% 
  mutate(
    incidence_rate_log = log(incidence_rate_add),
    incidence_rate_scaled = as.numeric(scale(incidence_rate_log))
  ) %>% 
  ungroup() %>% 
  arrange(region, date) %>% 
  mutate(region=to_title_case(region)) %>% 
  mutate(region=ifelse(region=="Sumatera","Sumatra",region)) %>% 
  mutate(region=ifelse(region=="Java Bali","Java & Bali",region))

#### Process dengue surveillance data at national level
dengue_data_admin0 <- dengue_data_admin1 %>% 
  filter(idadmin1==99)

#### Combine region and national data
dengue_data_regnat <- bind_rows(dengue_data_region,dengue_data_admin0 %>% mutate(region="National")) %>% 
  mutate(region=factor(region,levels=c("Sumatra","Java & Bali","Kalimantan","Nusa Tenggara",
                                       "Sulawesi","Maluku","Papua","National")))

#### Calculate monthly averages across years
dengue_data_admin1_monthly_avg <- dengue_data_admin1 %>% 
  group_by(region, idadmin1, admin1, month, month_name) %>% 
  summarise(incidence_rate = median(incidence_rate)) %>% 
  ungroup() %>% 
  arrange(idadmin1, month) %>% 
  group_by(admin1) %>% 
  mutate(
    incidence_rate_log = log(incidence_rate),
    incidence_rate_scaled = as.numeric(scale(incidence_rate))
  ) %>% 
  ungroup()

#### Calculate yearly summaries
dengue_data_admin1_yearly <- dengue_data_admin1 %>% 
  group_by(region, idadmin1, admin1, year) %>% 
  summarise(
    cases = sum(cases, na.rm = TRUE),
    pop = mean(pop)
  ) %>% 
  ungroup() %>% 
  mutate(
    incidence_rate = (cases / pop) * 100000
  )

#### Calculate overall 2010-2024 average incidence rates
dengue_data_admin1_overall_2010_2024 <- dengue_data_admin1_yearly %>% 
  group_by(region, idadmin1, admin1) %>% 
  summarise(incidence_rate = mean(incidence_rate)) %>% 
  ungroup() %>% 
  arrange(idadmin1)

#### Calculate summary statistics including coefficient of variation
dengue_data_admin1_overall_summary <- dengue_data_admin1 %>% 
  mutate(month = month(date)) %>% 
  group_by(region, idadmin1, admin1) %>% 
  summarise(
    median_incidence_rate = median(incidence_rate),
    mean_incidence_rate = mean(incidence_rate),
    sd_incidence_rate = sd(incidence_rate),
    cv_incidence_rate = sd_incidence_rate / mean_incidence_rate
  ) %>% 
  ungroup()

#### Prepare data for boxplot analysis by regions
dengue_data_admin1_boxplot <- dengue_data_admin1 %>% 
  mutate(month = factor(format(date, "%B"), levels = month.name)) %>% 
  left_join(admin1_national_EN) %>%
  mutate(shapeName = factor(as.character(shapeName), levels = rev(levels(shapeName)))) %>% 
  mutate(region2 = case_when(
    region %in% c("NUSA TENGGARA", "MALUKU", "PAPUA") ~ "NUSA TENGGARA, MALUKU, & PAPUA",
    region == "SUMATERA" ~ "SUMATERA",
    region == "JAVA & BALI" ~ "JAVA & BALI",
    region == "KALIMANTAN" ~ "KALIMANTAN",
    region == "SULAWESI" ~ "SULAWESI",
    TRUE ~ region
  )) %>% 
  arrange(idadmin1, year)

#### Prepare spatial data for mapping
admin1_ir_shp <- admin1_shp %>% 
  left_join(dengue_data_admin1_overall_2010_2024) %>% 
  left_join(admin1_indonesia)

#### Calculate centroids for map labels
admin1_ir_centroids <- admin1_ir_shp %>%
  st_centroid() %>%
  cbind(st_coordinates(.)) %>%
  mutate(shapeName = admin1_ir_shp$shapeName) %>%
  st_drop_geometry()

#### Prepare bivariate mapping data
admin1_incidence_rate_shp <- admin1_shp %>% 
  left_join(dengue_data_admin1_overall_summary)

admin1_incidence_rate_centroids <- admin1_incidence_rate_shp %>%
  st_centroid() %>%
  cbind(st_coordinates(.)) %>%
  mutate(shapeName = admin1_incidence_rate_shp$shapeName) %>%
  st_drop_geometry()

#### CLIMATE DATA PROCESSING
#### Read and process climate data
#### Downloaded from Copernicus Climate Data Store: https://cds.climate.copernicus.eu/datasets
#### Using KrigR package then aggregated to monthly data based on admin levels

#### Precipitation ERA5 Land data
prec_admin0_era_2009_2024 <- readRDS("data/prec_monthly_admin0.rds") %>% 
  dplyr::select(admin1, date, precipitation_mm_pop_weighted) %>% 
  left_join(admin1_indonesia) %>% 
  left_join(admin1_national_EN) %>% 
  dplyr::select(region, shapeName, idadmin1, admin1, date, precipitation = precipitation_mm_pop_weighted) %>% 
  as_tibble() %>% 
  arrange(idadmin1, date)

prec_admin1_era_2009_2024 <- readRDS("data/prec_monthly_admin1.rds") %>% 
  dplyr::select(idadmin1, date, precipitation_mm_pop_weighted) %>% 
  left_join(admin1_indonesia) %>% 
  left_join(admin1_national_EN) %>% 
  dplyr::select(region, shapeName, idadmin1, admin1, date, precipitation = precipitation_mm_pop_weighted) %>% 
  as_tibble() %>% 
  arrange(idadmin1, date)

prec_admin1_era_2009_2024 <- bind_rows(prec_admin1_era_2009_2024, prec_admin0_era_2009_2024)

#### Temperature ERA5 Land data
temp_admin0_era_2009_2024 <- readRDS("data/temp_monthly_admin0.rds") %>% 
  dplyr::select(admin1, date, temperature_celsius_pop_weighted) %>% 
  left_join(admin1_indonesia) %>% 
  left_join(admin1_national_EN) %>% 
  dplyr::select(region, shapeName, idadmin1, admin1, date, temperature = temperature_celsius_pop_weighted) %>% 
  as_tibble() %>% 
  arrange(idadmin1, date)

temp_admin1_era_2009_2024 <- readRDS("data/temp_monthly_admin1.rds") %>% 
  dplyr::select(shapeName, idadmin1, date, temperature_celsius_pop_weighted) %>% 
  left_join(admin1_indonesia) %>% 
  left_join(admin1_national_EN) %>% 
  dplyr::select(region, shapeName, idadmin1, admin1, date, temperature = temperature_celsius_pop_weighted) %>% 
  as_tibble() %>% 
  arrange(idadmin1, date)

temp_admin1_era_2009_2024 <- bind_rows(temp_admin1_era_2009_2024, temp_admin0_era_2009_2024)

#### Relative humidity ERA5 Land data
humi_admin0_era_2009_2024 <- readRDS("data/humid_monthly_admin0.rds") %>% 
  dplyr::select(admin1, date, rhumidity_pop_weighted) %>% 
  left_join(admin1_indonesia) %>% 
  left_join(admin1_national_EN) %>% 
  dplyr::select(region, shapeName, idadmin1, admin1, date, rel_humidity = rhumidity_pop_weighted) %>% 
  as_tibble() %>% 
  arrange(idadmin1, date)

humi_admin1_era_2009_2024 <- readRDS("data/humid_monthly_admin1.rds") %>% 
  dplyr::select(shapeName, idadmin1, date, rhumidity_pop_weighted) %>% 
  left_join(admin1_indonesia) %>% 
  left_join(admin1_national_EN) %>% 
  dplyr::select(region, shapeName, idadmin1, admin1, date, rel_humidity = rhumidity_pop_weighted) %>% 
  as_tibble() %>% 
  arrange(idadmin1, date)

humi_admin1_era_2009_2024 <- bind_rows(humi_admin1_era_2009_2024, humi_admin0_era_2009_2024)

#### Oceanic Nino Index
#### The Oceanic Nino Index is average sea surface temperature 
#### in the Nino 3.4 region (120W to 170W) averaged over three months
oni_2009_2024 <- download_oni() %>% 
  clean_names() %>% 
  filter(year(date) >= 2009) %>% 
  arrange(date)

iod_2009_2024 <- download_dmi() %>% 
  clean_names() %>% 
  filter(year(date) >= 2009) %>% 
  arrange(date)

oni_iod_2009_2024 <- oni_2009_2024 %>% 
  left_join(iod_2009_2024 %>% dplyr::select(date,dmi))

#### Create lagged climate data up to 3 months
#### Precipitation lagged data
prec_admin1_era_2009_2024_lag1 <- prec_admin1_era_2009_2024 %>% 
  mutate(date = date %m+% months(1)) %>% 
  dplyr::select(region, idadmin1, admin1, shapeName, date, precipitation_era_lag1 = precipitation)

prec_admin1_era_2009_2024_lag2 <- prec_admin1_era_2009_2024 %>% 
  mutate(date = date %m+% months(2)) %>% 
  dplyr::select(region, idadmin1, admin1, shapeName, date, precipitation_era_lag2 = precipitation)

prec_admin1_era_2009_2024_lag3 <- prec_admin1_era_2009_2024 %>% 
  mutate(date = date %m+% months(3)) %>% 
  dplyr::select(region, idadmin1, admin1, shapeName, date, precipitation_era_lag3 = precipitation)

#### Temperature lagged data
temp_admin1_era_2009_2024_lag1 <- temp_admin1_era_2009_2024 %>% 
  mutate(date = date %m+% months(1)) %>% 
  dplyr::select(region, idadmin1, admin1, shapeName, date, temperature_era_lag1 = temperature)

temp_admin1_era_2009_2024_lag2 <- temp_admin1_era_2009_2024 %>% 
  mutate(date = date %m+% months(2)) %>% 
  dplyr::select(region, idadmin1, admin1, shapeName, date, temperature_era_lag2 = temperature)

temp_admin1_era_2009_2024_lag3 <- temp_admin1_era_2009_2024 %>% 
  mutate(date = date %m+% months(3)) %>% 
  dplyr::select(region, idadmin1, admin1, shapeName, date, temperature_era_lag3 = temperature)

#### Relative humidity lagged data
humi_admin1_era_2009_2024_lag1 <- humi_admin1_era_2009_2024 %>% 
  mutate(date = date %m+% months(1)) %>% 
  dplyr::select(region, idadmin1, admin1, shapeName, date, rel_humidity_era_lag1 = rel_humidity)

humi_admin1_era_2009_2024_lag2 <- humi_admin1_era_2009_2024 %>% 
  mutate(date = date %m+% months(2)) %>% 
  dplyr::select(region, idadmin1, admin1, shapeName, date, rel_humidity_era_lag2 = rel_humidity)

humi_admin1_era_2009_2024_lag3 <- humi_admin1_era_2009_2024 %>% 
  mutate(date = date %m+% months(3)) %>% 
  dplyr::select(region, idadmin1, admin1, shapeName, date, rel_humidity_era_lag3 = rel_humidity)

#### ONI lagged data
oni_2009_2024_lag1 <- oni_2009_2024 %>% 
  mutate(date = date %m+% months(1)) %>% 
  dplyr::select(date, oni_lag1 = oni)

oni_2009_2024_lag2 <- oni_2009_2024 %>% 
  mutate(date = date %m+% months(2)) %>% 
  dplyr::select(date, oni_lag2 = oni)

oni_2009_2024_lag3 <- oni_2009_2024 %>% 
  mutate(date = date %m+% months(3)) %>% 
  dplyr::select(date, oni_lag3 = oni)

#### Combine all climate data with lags
prec_admin1_era_2010_2024 <- prec_admin1_era_2009_2024 %>% 
  dplyr::select(region, idadmin1, admin1, shapeName, date, precipitation_era = precipitation) %>% 
  left_join(prec_admin1_era_2009_2024_lag1) %>% 
  left_join(prec_admin1_era_2009_2024_lag2) %>% 
  left_join(prec_admin1_era_2009_2024_lag3) %>% 
  filter(year(date) %in% 2010:2024)

temp_admin1_era_2010_2024 <- temp_admin1_era_2009_2024 %>% 
  dplyr::select(region, idadmin1, admin1, shapeName, date, temperature_era = temperature) %>% 
  left_join(temp_admin1_era_2009_2024_lag1) %>%
  left_join(temp_admin1_era_2009_2024_lag2) %>%
  left_join(temp_admin1_era_2009_2024_lag3) %>%
  filter(year(date) %in% 2010:2024)

humi_admin1_era_2010_2024 <- humi_admin1_era_2009_2024 %>% 
  dplyr::select(region, idadmin1, admin1, shapeName, date, rel_humidity_era = rel_humidity) %>% 
  left_join(humi_admin1_era_2009_2024_lag1) %>%
  left_join(humi_admin1_era_2009_2024_lag2) %>%
  left_join(humi_admin1_era_2009_2024_lag3) %>%
  filter(year(date) %in% 2010:2024)

oni_2010_2024 <- oni_2009_2024 %>% 
  dplyr::select(date, oni) %>% 
  left_join(oni_2009_2024_lag1) %>%
  left_join(oni_2009_2024_lag2) %>%
  left_join(oni_2009_2024_lag3) %>%
  filter(year(date) %in% 2010:2024)

#### Create comprehensive climate dataset
climate_2010_2024 <- prec_admin1_era_2010_2024 %>% 
  left_join(temp_admin1_era_2010_2024) %>% 
  left_join(humi_admin1_era_2010_2024) %>% 
  left_join(oni_2010_2024)

#### Merge climate and dengue incidence data
climate_incidence_data <- climate_2010_2024 %>%
  left_join(dengue_data_admin1 %>% 
              dplyr::select(idadmin1, admin1, date, incidence_rate, incidence_rate_scaled), 
            by = c("idadmin1", "admin1", "date")) %>% 
  group_by(region, idadmin1, admin1) %>% 
  mutate(
    precipitation_era_scaled = as.numeric(scale(log(precipitation_era))),
    precipitation_era_lag1_scaled = as.numeric(scale(log(precipitation_era_lag1))),
    precipitation_era_lag2_scaled = as.numeric(scale(log(precipitation_era_lag2))),
    precipitation_era_lag3_scaled = as.numeric(scale(log(precipitation_era_lag3))),
    temperature_era_scaled = as.numeric(scale(temperature_era)),
    temperature_era_lag1_scaled = as.numeric(scale(temperature_era_lag1)),
    temperature_era_lag2_scaled = as.numeric(scale(temperature_era_lag2)),
    temperature_era_lag3_scaled = as.numeric(scale(temperature_era_lag3)),
    rel_humidity_era_scaled = as.numeric(scale(rel_humidity_era)),
    rel_humidity_era_lag1_scaled = as.numeric(scale(rel_humidity_era_lag1)),
    rel_humidity_era_lag2_scaled = as.numeric(scale(rel_humidity_era_lag2)),
    rel_humidity_era_lag3_scaled = as.numeric(scale(rel_humidity_era_lag3))
  ) %>% 
  ungroup()

#### Create climate data without lags for visualization
climate_era_2009_2024_nolag <- prec_admin1_era_2009_2024 %>% 
  left_join(temp_admin1_era_2009_2024) %>% 
  left_join(humi_admin1_era_2009_2024) %>% 
  left_join(oni_2009_2024) %>% 
  dplyr::select(region, idadmin1, admin1, date,
                precipitation_era = precipitation, temperature_era = temperature, 
                rel_humidity_era = rel_humidity, oni) %>% 
  group_by(admin1) %>% 
  mutate(
    precipitation_era_scaled = as.numeric(scale(precipitation_era)),
    temperature_era_scaled = as.numeric(scale(temperature_era)),
    rel_humidity_era_scaled = as.numeric(scale(rel_humidity_era))
  )

#### WAVELET ANALYSIS
#### Wavelet analysis preparation and execution
# Extract unique admin1 names
admin1_names <- admin1_indonesia %>% pull(admin1) %>% unique()

# Prepare province-level dataset for wavelet analysis
# Filter out 4 provinces in Papua region
province_data <- dengue_data_admin1 %>% filter(idadmin1 != 99)

# Prepare data for wavelet analysis at the province level
province_ts <- prepare_for_wavelet(province_data, "admin1", "incidence_rate_scaled")

#### Run wavelet analysis at the province level
province_wavelet <- conduct_wavelet_analysis(
  province_ts$ts_data,
  province_ts$valid_groups
)

#### Calculate phase differences at the province level
province_phase_months <- create_phase_matrix(
  province_wavelet$phase_angles,
  province_ts$valid_groups
)

#### Calculate average phase lags
province_phase_lags <- calculate_phase_lags(province_phase_months)

#### Create data frames with phase lags for visualization
province_lags_df <- data.frame(
  province = names(province_phase_lags),
  phase_lag = province_phase_lags,
  region = province_data$region[match(names(province_phase_lags), province_data$admin1)]
)

#### Calculate phase statistics from the phase difference matrix
province_phase_stats <- calculate_phase_stats(province_phase_months) 
province_phase_stats <- province_phase_stats %>% 
  left_join(admin1_indonesia %>% 
              dplyr::select(region, idadmin1, province = admin1) %>% 
              mutate(province = factor(province, levels = province_phase_stats$province))) %>% 
  mutate(region=to_title_case(region),
         region=ifelse(region=="Java Bali","Java & Bali",region),
         region=factor(region,levels=c("Sumatra","Java & Bali","Kalimantan","Nusa Tenggara",
                                       "Sulawesi","Maluku","Papua")))

#### Join phase lag data with province shapefile for mapping
admin1_shp_with_lags <- admin1_shp %>%
  left_join(province_phase_stats)

#### Calculate centroids for phase lag map labels
admin1_phase_lags_centroids <- admin1_shp_with_lags %>%
  st_centroid() %>%
  cbind(st_coordinates(.)) %>%
  mutate(shapeName = admin1_shp_with_lags$shapeName) %>%
  st_drop_geometry() %>% 
  left_join(admin1_indonesia) %>% 
  arrange(idadmin1)

#### Find peak months based on province phase
peak_months_admin1_list <- list()
for (i in seq_len(34)){
  
  peak_months <- find_peak_months(province_wavelet$phase_angles[[i]])
  peak_months_admin1_list[[i]] <- admin1_indonesia[i,] %>% 
    mutate(average_peak_month=peak_months$average_peak_month,
           average_peak_month_name=peak_months$average_peak_month_name)
  
}
peak_months_admin1 <- bind_rows(peak_months_admin1_list) %>% 
  left_join(admin1_EN) %>% 
  mutate(average_peak_month_name=factor(average_peak_month_name,levels=c("October","November","December",
                                                                         "January","February","March","April")))

#### Year-by-year analysis
#### Apply to all provinces
all_provinces_peaks <- data.frame()

for(i in seq_len(34)) {
  
  # Find peaks
  peak_results <- find_peak_months(province_wavelet$phase_angles[[i]])
  
  # Store year-by-year data
  province_yearly <- peak_results$peak_data %>%
    mutate(idadmin1 = idadmin1_indonesia_unique[i],
           admin1 = admin1_indonesia_unique[i],
           lat = admin1_phase_lags_centroids$Y[i],
           lon = admin1_phase_lags_centroids$X[i])
  
  all_provinces_peaks <- rbind(all_provinces_peaks, province_yearly) %>% as_tibble()
}

# 1. Year-by-year west-to-east gradient
# Calculate west-east gradient correlation for each epidemic year
# Using two-sided Spearman test (exploratory analysis)
yearly_gradient <- all_provinces_peaks %>%
  filter(Period != "2009-2010") %>%
  group_by(Period) %>%
  summarise(
    cor_longitude = cor.test(lon, Month_JulJun_Adjusted, use = "complete.obs",
                             method = "spearman", exact = FALSE,
                             alternative = "two.sided")$estimate,
    p_value = cor.test(lon, Month_JulJun_Adjusted, use = "complete.obs",
                       method = "spearman", exact = FALSE,
                       alternative = "two.sided")$p.value
  ) %>%
  filter(Period != "2024-2025")

# Add significance categories
yearly_gradient <- yearly_gradient %>%
  mutate(
    significance = case_when(
      p_value < 0.01 ~ "p < 0.01",
      p_value < 0.05 ~ "p < 0.05",
      p_value < 0.10 ~ "p < 0.10",
      TRUE ~ "n.s."
    ),
    year_numeric = as.numeric(substr(Period, 1, 4))
  )

print(yearly_gradient)

all_provinces_peaks_yearly <- all_provinces_peaks %>% 
  left_join(yearly_gradient, by = "Period")

# Calculate west-east gradient for western provinces only (sensitivity analysis)
# Sumatra, Java-Bali, and Kalimantan
yearly_gradient_west <- all_provinces_peaks %>%
  filter(Period != "2009-2010") %>%
  filter(idadmin1 %in% idadmin1_indonesia_unique[c(1:17,20:24)]) %>%
  group_by(Period) %>%
  summarise(
    cor_longitude = cor.test(lon, Month_JulJun_Adjusted, use = "complete.obs",
                             method = "spearman", exact = FALSE,
                             alternative = "two.sided")$estimate,
    p_value = cor.test(lon, Month_JulJun_Adjusted, use = "complete.obs",
                       method = "spearman", exact = FALSE,
                       alternative = "two.sided")$p.value
  ) %>%
  filter(Period != "2024-2025")

# Add significance categories
yearly_gradient_west <- yearly_gradient_west %>%
  mutate(
    significance = case_when(
      p_value < 0.01 ~ "p < 0.01",
      p_value < 0.05 ~ "p < 0.05",
      p_value < 0.10 ~ "p < 0.10",
      TRUE ~ "n.s."
    ),
    year_numeric = as.numeric(substr(Period, 1, 4))
  )

print(yearly_gradient_west)

# Calculate west-east gradient for western provinces only (sensitivity analysis)
# Sumatra, Java-Bali, without Kalimantan
yearly_gradient_west2 <- all_provinces_peaks %>%
  filter(Period != "2009-2010") %>%
  filter(idadmin1 %in% idadmin1_indonesia_unique[c(1:17)]) %>%
  group_by(Period) %>%
  summarise(
    cor_longitude = cor.test(lon, Month_JulJun_Adjusted, use = "complete.obs",
                             method = "spearman", exact = FALSE,
                             alternative = "two.sided")$estimate,
    p_value = cor.test(lon, Month_JulJun_Adjusted, use = "complete.obs",
                       method = "spearman", exact = FALSE,
                       alternative = "two.sided")$p.value
  ) %>%
  filter(Period != "2024-2025")

# Add significance categories
yearly_gradient_west2 <- yearly_gradient_west2 %>%
  mutate(
    significance = case_when(
      p_value < 0.01 ~ "p < 0.01",
      p_value < 0.05 ~ "p < 0.05",
      p_value < 0.10 ~ "p < 0.10",
      TRUE ~ "n.s."
    ),
    year_numeric = as.numeric(substr(Period, 1, 4))
  )

print(yearly_gradient_west2)

# 2. Within-province consistency
province_consistency <- all_provinces_peaks %>%
  filter(Period != "2009-2010") %>% 
  group_by(idadmin1,admin1) %>%
  summarise(
    median_peak = month_name_juljun[round(median(Month_JulJun,na.rm=TRUE))],
    median_peak_JulJun = round(median(Month_JulJun,na.rm=TRUE)),
    sd_peak = sd(Month_JulJun,na.rm=TRUE),
    range_peak = paste(month_abb_juljun[min(Month_JulJun,na.rm=TRUE)], "-", month_abb_juljun[max(Month_JulJun,na.rm=TRUE)]),
    cv = sd(Month_JulJun,na.rm=TRUE) / mean(Month_JulJun,na.rm=TRUE)
  ) %>% 
  ungroup() %>% 
  mutate(median_peak=factor(median_peak,levels=c("October","November","December",
                                                 "January","February","March","April")))

#### Join peak month data with province shapefile for mapping
admin1_shp_with_peak <- admin1_shp %>%
  left_join(province_consistency)

#### Calculate centroids for phase lag map labels
admin1_peak_centroids <- admin1_shp_with_peak %>%
  st_centroid() %>%
  cbind(st_coordinates(.)) %>%
  mutate(shapeName = admin1_shp_with_peak$shapeName) %>%
  st_drop_geometry() %>% 
  left_join(admin1_indonesia) %>% 
  arrange(idadmin1)

# 3. Identify anomalous years (like 2020)
deviation_from_peaks <- all_provinces_peaks %>%
  filter(Period != "2009-2010") %>% 
  left_join(province_consistency %>% dplyr::select(idadmin1,admin1,median_peak,median_peak_JulJun)) %>% 
  group_by(idadmin1, admin1) %>%
  mutate(
    deviation_from_median = abs(Month_JulJun - median_peak_JulJun)
  ) %>%
  dplyr::select(idadmin1,admin1,Period,deviation_from_median,Month_JulJun,Month_Name,median_peak_JulJun,median_peak) %>% 
  ungroup()

# some wrong calculations due to circularity, but manual edit for the manuscript, i.e., Bengkulu 2021-2022 should be 4 months instead of 7 months difference

deviation_from_peaks[deviation_from_peaks$admin1=="BENGKULU" & deviation_from_peaks$Period=="2021-2022",]$deviation_from_median <- 4
deviation_from_peaks[deviation_from_peaks$admin1=="MALUKU UTARA" & deviation_from_peaks$Period=="2017-2018",]$deviation_from_median <- 4

mean_deviation_from_peaks <- deviation_from_peaks %>% 
  group_by(idadmin1, admin1) %>%
  summarise(mean_deviation_from_peaks=round(mean(deviation_from_median,na.rm=TRUE),2),
            median_deviation_from_peaks=median(deviation_from_median,na.rm=TRUE))

range_deviation_from_peaks <- deviation_from_peaks %>% 
  group_by(idadmin1, admin1) %>%
  summarise(range_deviation_from_peaks=paste0(min(deviation_from_median,na.rm=TRUE),"-",max(deviation_from_median,na.rm=TRUE)))

province_consistency <- province_consistency %>% 
  left_join(mean_deviation_from_peaks) %>% 
  left_join(range_deviation_from_peaks)

anomalous_peaks <- deviation_from_peaks %>% 
  filter(deviation_from_median > 3)  # More than 3 months deviation

print(anomalous_peaks)

#### start of DTW analysis
# ============================================================================
# DATA LOADING AND PREPARATION FOR DTW
# ============================================================================

# Assuming data are already loaded
# dengue_data_admin1 <- read.csv("path/to/dengue_data_admin1.csv") %>% as_tibble()
# climate_admin1 <- read.csv("path/to/climate_admin1.csv") %>% as_tibble()

# Imputation for case data
# Calculate median monthly proportion of cases
# In years when there's NA, extrapolate total cases that year using 
# month with highest number of cases divided by median proportion
# Then calculate the NA values using the median proportion of cases that month

dengue_data_admin1_noNA <- dengue_data_admin1 %>% 
  filter(cases != 0)
dengue_data_admin1_noNA_complete_year <- dengue_data_admin1_noNA %>% 
  group_by(year,idadmin1) %>% 
  mutate(count=n()) %>% 
  filter(count == 12)
dengue_data_admin1_monthly_prop <- dengue_data_admin1_noNA_complete_year %>% 
  group_by(region,idadmin1,admin1,month) %>% 
  summarise(cases=sum(cases)) %>% 
  ungroup() %>% 
  group_by(region,idadmin1,admin1) %>% 
  mutate(prop_cases=cases/sum(cases)) %>% 
  ungroup() %>% 
  dplyr::select(-cases)

dengue_data_admin1_NA <- dengue_data_admin1 %>% 
  filter(cases == 0)
dengue_data_admin1_NA_year <- dengue_data_admin1_NA %>% 
  dplyr::select(idadmin1,year) %>% 
  unique()
dengue_data_admin1_NA_highest_cases_in_that_year <- dengue_data_admin1_NA_year %>% 
  left_join(dengue_data_admin1) %>%
  group_by(idadmin1,year) %>% 
  summarise(max_cases=max(cases))
dengue_data_admin1_NA_month_with_highest_cases_in_that_year <- dengue_data_admin1_NA_year %>% 
  left_join(dengue_data_admin1) %>% 
  left_join(dengue_data_admin1_NA_highest_cases_in_that_year) %>% 
  filter(cases == max_cases) %>%  # use feb for 76, and mar for 81
  filter(!(idadmin1 == 76 & month == 4 & year == 2013) & !(idadmin1 == 81 & month == 2 & year == 2019)) %>% 
  left_join(dengue_data_admin1_monthly_prop) %>% 
  mutate(total_cases_est = round(max_cases/prop_cases)) %>% 
  dplyr::select(region,idadmin1,admin1,year,month,max_cases,prop_cases,total_cases_est)

dengue_data_admin1_NA_imputed <- dengue_data_admin1_NA %>% 
  dplyr::select(region,idadmin1,admin1,year,month,pop) %>% 
  left_join(dengue_data_admin1_monthly_prop) %>% 
  left_join(dengue_data_admin1_NA_month_with_highest_cases_in_that_year %>% 
              dplyr::select(region,idadmin1,admin1,year,total_cases_est)) %>% 
  mutate(cases=round(total_cases_est * prop_cases)) %>% 
  dplyr::select(region,idadmin1,admin1,year,month,cases,pop)

dengue_data_admin1_imputed <- bind_rows(dengue_data_admin1_noNA,dengue_data_admin1_NA_imputed) %>% 
  arrange(idadmin1,year,month) %>% 
  group_by(region, idadmin1, admin1, year, month) %>% 
  summarise(
    cases = sum(cases, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  left_join(pop_indonesia_admin1) %>% 
  mutate(
    incidence_rate = (cases / pop) * 100000,
    month_name = factor(month, levels = 1:12, 
                        labels = c("January", "February", "March", "April", "May", "June", 
                                   "July", "August", "September", "October", "November", "December")),
    date = ymd(paste0(year, "-", month, "-01"))
  ) %>% 
  group_by(admin1) %>% 
  mutate(
    incidence_rate_log = log(incidence_rate),
    incidence_rate_scaled = as.numeric(scale(incidence_rate_log))
  ) %>% 
  ungroup() %>% 
  arrange(idadmin1, date)

# Convert to tibble if needed
# dengue_data_admin1 <- dengue_data_admin1 %>% as_tibble()
climate_admin1 <- prec_admin1_era_2009_2024 %>% 
  left_join(temp_admin1_era_2009_2024) %>% 
  left_join(humi_admin1_era_2009_2024) %>% 
  left_join(oni_2009_2024 %>% dplyr::select(date,oni)) %>% 
  left_join(iod_2009_2024 %>% dplyr::select(date,dmi)) %>% 
  mutate(year=year(date),month=month(date)) %>% 
  dplyr::select(region,shapeName,idadmin1,admin1,year,month,date,everything())

# Create year-month identifier and merge datasets
merged_data <- dengue_data_admin1_imputed %>%
  mutate(yearmonth = str_glue("{year}-{sprintf('%02d', month)}")) %>%
  left_join(
    climate_admin1 %>%
      filter(year >= 2010) %>%
      mutate(yearmonth = str_glue("{year}-{sprintf('%02d', month)}")),
    by = c("region", "idadmin1" , "admin1", "year", "month", "yearmonth")
  ) %>%
  filter(year >= 2010) %>%
  dplyr::select(region, idadmin1, admin1, year, month, yearmonth, cases, 
                precipitation, temperature, rel_humidity, oni, dmi)

# ============================================================================
# 2. DATA NORMALIZATION AND TRANSFORMATION
# ============================================================================

normalized_data <- merged_data %>%
  filter(idadmin1 != 99) %>%  # remove national level data
  # Log-transform dengue cases
  mutate(cases_log = log(cases)) %>%
  # Log-transform precipitation (zero-inflated)
  mutate(precipitation_log = log(precipitation)) %>%
  # Z-score normalization within each province
  group_by(admin1) %>%
  mutate(
    cases_norm = scale(cases_log)[, 1],
    precipitation_norm = scale(precipitation_log)[, 1],
    temperature_norm = scale(temperature)[, 1],
    humidity_norm = scale(rel_humidity)[, 1],
    oni_norm = scale(oni)[, 1],
    dmi_norm = scale(dmi)[, 1]
  ) %>%
  ungroup() %>%
  dplyr::select(region, idadmin1, admin1, year, month, yearmonth, cases_norm, precipitation_norm, 
                temperature_norm, humidity_norm, oni_norm, dmi_norm)


#### Clustering
#### Cases
cases_dtw_ready <- normalized_data %>% 
  filter(idadmin1 != 99) %>% 
  arrange(idadmin1,yearmonth) %>% 
  dplyr::select(idadmin1,yearmonth,cases_norm) %>% 
  pivot_wider(names_from=idadmin1,values_from=cases_norm) %>% 
  dplyr::select(-yearmonth) %>% 
  t() %>% 
  tslist()

cases_clust <- tsclust(
  series = cases_dtw_ready, 
  type = "hierarchical",
  k = 2:15, 
  distance = "dtw_basic", 
  # distance="dtw",
  control=hierarchical_control(method="average"),
  seed = 73
)

names(cases_clust) <- paste0("k_", 2:15)
res_cvi <- sapply(cases_clust, cvi, type = "internal") %>% 
  t() %>%
  as_tibble(rownames = "k") %>%
  arrange(-Sil)

cases_sel_clust <- cases_clust[[res_cvi[[6,1]]]]

plot(cases_sel_clust)
abline(h = cases_clust$k_12$height[length(cases_clust$k_12$height) - (12 - 1)], 
       col = "red", lty = 2, lwd = 2, label = paste0("Cut for ", 12, " clusters"))

table(cases_sel_clust@cluster)

cases_admin1_cluster <- tibble(
  idadmin1 = unique(normalized_data$idadmin1),
  admin1 = unique(normalized_data$admin1),
  group = as.character(cases_sel_clust@cluster)
)

#### Clustering only western provinces
#### Sumatra, Java-Bali, and Kalimantan
cases_dtw_ready_west <- normalized_data %>% 
  filter(idadmin1 %in% idadmin1_indonesia_unique[c(1:17,20:24)]) %>% 
  arrange(idadmin1,yearmonth) %>% 
  dplyr::select(idadmin1,yearmonth,cases_norm) %>% 
  pivot_wider(names_from=idadmin1,values_from=cases_norm) %>% 
  dplyr::select(-yearmonth) %>% 
  t() %>% 
  tslist()

cases_clust_west <- tsclust(
  series = cases_dtw_ready_west, 
  type = "hierarchical",
  k = 2:15, 
  distance = "dtw_basic", 
  # distance="dtw",
  control=hierarchical_control(method="average"),
  seed = 73
)

names(cases_clust_west) <- paste0("k_", 2:15)
res_cvi_west <- sapply(cases_clust_west, cvi, type = "internal") %>% 
  t() %>%
  as_tibble(rownames = "k") %>%
  arrange(-Sil)

cases_sel_clust_west <- cases_clust_west[[res_cvi_west[[3,1]]]]

plot(cases_sel_clust_west)
abline(h = cases_clust_west$k_10$height[length(cases_clust_west$k_10$height) - (10 - 1)], 
       col = "red", lty = 2, lwd = 2, label = paste0("Cut for ", 10, " clusters"))

table(cases_sel_clust_west@cluster)

cases_admin1_cluster_west <- tibble(
  idadmin1 = unique(normalized_data$idadmin1)[c(1:17,20:24)],
  admin1 = unique(normalized_data$admin1)[c(1:17,20:24)],
  group = as.character(cases_sel_clust_west@cluster)
)

#### Clustering only western provinces
#### Sumatra, Java-Bali, without Kalimantan
cases_dtw_ready_west2 <- normalized_data %>% 
  filter(idadmin1 %in% idadmin1_indonesia_unique[1:17]) %>% 
  arrange(idadmin1,yearmonth) %>% 
  dplyr::select(idadmin1,yearmonth,cases_norm) %>% 
  pivot_wider(names_from=idadmin1,values_from=cases_norm) %>% 
  dplyr::select(-yearmonth) %>% 
  t() %>% 
  tslist()

cases_clust_west2 <- tsclust(
  series = cases_dtw_ready_west2, 
  type = "hierarchical",
  k = 2:15, 
  distance = "dtw_basic", 
  # distance="dtw",
  control=hierarchical_control(method="average"),
  seed = 73
)

names(cases_clust_west2) <- paste0("k_", 2:15)
res_cvi_west2 <- sapply(cases_clust_west2, cvi, type = "internal") %>% 
  t() %>%
  as_tibble(rownames = "k") %>%
  arrange(-Sil)

cases_sel_clust_west2 <- cases_clust_west2[[res_cvi_west2[[1,1]]]]

plot(cases_sel_clust_west2)
abline(h = cases_clust_west2$k_4$height[length(cases_clust_west2$k_4$height) - (4 - 1)], 
       col = "red", lty = 2, lwd = 2, label = paste0("Cut for ", 4, " clusters"))

table(cases_sel_clust_west2@cluster)

cases_admin1_cluster_west2 <- tibble(
  idadmin1 = unique(normalized_data$idadmin1)[1:17],
  admin1 = unique(normalized_data$admin1)[1:17],
  group = as.character(cases_sel_clust_west2@cluster)
)

#### Silhouette table
write_xlsx(
  list(all_province=res_cvi %>%
         mutate(k = as.numeric(str_remove(k, "^k_"))) %>% 
         dplyr::select(`No. of clusters`=k,Silhouette=Sil),
       west_with_kalimantan=res_cvi_west %>%
         mutate(k = as.numeric(str_remove(k, "^k_"))) %>% 
         dplyr::select(`No. of clusters`=k,Silhouette=Sil),
       west_without_kalimantan=res_cvi_west2 %>%
         mutate(k = as.numeric(str_remove(k, "^k_"))) %>% 
         dplyr::select(`No. of clusters`=k,Silhouette=Sil)),
  "output/sihouette_cluster.xlsx"
)

# Note: Climate-based DTW clustering was explored but not used in final analysis.
# The dengue case clustering (above) better captures epidemiological patterns.

################################################################################
# SECTION: WITHIN-REGION COHESION AND OUTLIER ANALYSIS
################################################################################

# Extract DTW distance matrix from clustering object
# The distance matrix is stored in the clustering object
dtw_dist_matrix <- cases_clust$k_12@distmat
rownames(dtw_dist_matrix) <- colnames(dtw_dist_matrix) <- unique(normalized_data$idadmin1)

# Create province-region mapping
province_region_map <- admin1_indonesia %>%
  filter(idadmin1 != 99) %>%
  dplyr::select(idadmin1, admin1, region) %>%
  mutate(region = to_title_case(region),
         region = ifelse(region == "Java Bali", "Java & Bali", region),
         region = ifelse(region == "Sumatera", "Sumatra", region))

# Function to calculate within-region cohesion metrics
calculate_region_cohesion <- function(dist_matrix, province_region_map) {

  regions <- unique(province_region_map$region)

  cohesion_results <- lapply(regions, function(reg) {
    # Get provinces in this region
    region_provinces <- province_region_map %>%
      filter(region == reg) %>%
      pull(idadmin1) %>%
      as.character()

    n_provinces <- length(region_provinces)

    if (n_provinces < 2) {
      return(tibble(
        region = reg,
        n_provinces = n_provinces,
        mean_dtw_distance = NA_real_,
        sd_dtw_distance = NA_real_,
        cv_dtw_distance = NA_real_,
        min_dtw_distance = NA_real_,
        max_dtw_distance = NA_real_
      ))
    }

    # Extract submatrix for this region
    region_dist <- dist_matrix[region_provinces, region_provinces]

    # Get lower triangle (pairwise distances, excluding diagonal)
    pairwise_distances <- region_dist[lower.tri(region_dist)]

    tibble(
      region = reg,
      n_provinces = n_provinces,
      mean_dtw_distance = mean(pairwise_distances),
      sd_dtw_distance = sd(pairwise_distances),
      cv_dtw_distance = sd(pairwise_distances) / mean(pairwise_distances),
      min_dtw_distance = min(pairwise_distances),
      max_dtw_distance = max(pairwise_distances)
    )
  })

  bind_rows(cohesion_results)
}

# Calculate region cohesion
region_cohesion <- calculate_region_cohesion(dtw_dist_matrix, province_region_map)
print("Within-region cohesion (lower mean DTW distance = more similar patterns):")
print(region_cohesion %>% arrange(mean_dtw_distance))

# Function to identify outlier provinces within each region
identify_region_outliers <- function(dist_matrix, province_region_map, z_threshold = 1.645) {

  regions <- unique(province_region_map$region)

  outlier_results <- lapply(regions, function(reg) {
    # Get provinces in this region
    region_provinces <- province_region_map %>%
      filter(region == reg) %>%
      pull(idadmin1) %>%
      as.character()

    n_provinces <- length(region_provinces)

    if (n_provinces < 3) {
      return(NULL)
    }

    # Extract submatrix for this region
    region_dist <- dist_matrix[region_provinces, region_provinces]

    # Calculate mean distance to other provinces within region for each province
    mean_dist_to_others <- sapply(seq_len(n_provinces), function(i) {
      mean(region_dist[i, -i])
    })
    names(mean_dist_to_others) <- region_provinces

    # Calculate z-scores
    mean_overall <- mean(mean_dist_to_others)
    sd_overall <- sd(mean_dist_to_others)
    z_scores <- (mean_dist_to_others - mean_overall) / sd_overall

    # Identify outliers
    is_outlier <- abs(z_scores) > z_threshold

    tibble(
      region = reg,
      idadmin1 = as.numeric(region_provinces),
      mean_dist_to_region = mean_dist_to_others,
      z_score = z_scores,
      is_outlier = is_outlier
    ) %>%
      left_join(province_region_map %>% dplyr::select(idadmin1, admin1), by = "idadmin1")
  })

  bind_rows(outlier_results)
}

# Identify outliers within regions
region_outliers <- identify_region_outliers(dtw_dist_matrix, province_region_map)
print("Outlier provinces within regions (z-score > 2):")
print(region_outliers %>% filter(is_outlier) %>% arrange(region, desc(z_score)))

# Summary: provinces with highest dissimilarity to their region
province_dissimilarity <- region_outliers %>%
  arrange(desc(mean_dist_to_region)) %>%
  left_join(region_cohesion %>% dplyr::select(region, region_mean = mean_dtw_distance), by = "region") %>%
  mutate(relative_dissimilarity = mean_dist_to_region / region_mean)

print("Top 10 most dissimilar provinces relative to their region:")
print(province_dissimilarity %>% head(10) %>% dplyr::select(region, admin1, mean_dist_to_region, z_score, relative_dissimilarity))

################################################################################
# SECTION: CLIMATE-DENGUE PHASE RELATIONSHIP ANALYSIS (WAVELET-BASED)
################################################################################

# Prepare climate data for wavelet analysis (scaled within province)
climate_for_wavelet <- climate_era_2009_2024_nolag %>%
  filter(idadmin1 != 99, year(date) >= 2010, year(date) <= 2024) %>%
  arrange(idadmin1, date)

# Run wavelet analysis on climate variables to extract phase angles
climate_vars_wavelet <- c("precipitation_era_scaled", "temperature_era_scaled", "rel_humidity_era_scaled")

# Store climate phase angles
climate_phase_angles <- list()

for (clim_var in climate_vars_wavelet) {

  # Prepare data for wavelet
  climate_ts <- prepare_for_wavelet(climate_for_wavelet, "admin1", clim_var)

  # Run wavelet analysis
  climate_wavelet <- conduct_wavelet_analysis(
    climate_ts$ts_data,
    climate_ts$valid_groups
  )

  climate_phase_angles[[clim_var]] <- climate_wavelet$phase_angles

  print(paste0("Wavelet analysis completed for: ", clim_var))
}

# Function to calculate circular correlation between two phase angle series
# Using circular correlation coefficient
circular_correlation <- function(phase1, phase2) {
  # Remove NAs
  valid <- !is.na(phase1) & !is.na(phase2)
  p1 <- phase1[valid]
  p2 <- phase2[valid]

  if (length(p1) < 10) return(NA_real_)

  # Calculate circular correlation (Jammalamadaka & SenGupta)
  sin1 <- sin(p1 - mean(circular(p1, type = "angles", units = "radians")))
  sin2 <- sin(p2 - mean(circular(p2, type = "angles", units = "radians")))

  r <- sum(sin1 * sin2) / sqrt(sum(sin1^2) * sum(sin2^2))

  return(r)
}

# Function to calculate mean phase difference (climate leads dengue by how many months?)
calculate_phase_lag <- function(phase_climate, phase_dengue) {
  # Phase difference: positive = climate leads, negative = dengue leads
  diff <- phase_climate - phase_dengue

  # Wrap to [-pi, pi]
  diff <- ((diff + pi) %% (2 * pi)) - pi
  
  # Use circular mean to get representative phase lag
  # library(circular)
  mean_diff <- mean(circular(diff, type = "angles", units = "radians"))
  
  # Convert to months (2*pi = 12 months)
  mean_diff_months <- as.numeric(mean_diff) * 12 / (2 * pi)

  return(mean_diff_months)
}

# Calculate phase relationships between climate and dengue for each province
climate_dengue_phase <- lapply(seq_len(34), function(i) {

  province_name <- province_ts$valid_groups[i]

  # Get dengue phase angles for this province
  dengue_phase <- province_wavelet$phase_angles[[province_name]]

  # Calculate relationships with each climate variable
  var_results <- lapply(climate_vars_wavelet, function(clim_var) {

    climate_phase <- climate_phase_angles[[clim_var]][[province_name]]

    if (is.null(climate_phase) || length(climate_phase) != length(dengue_phase)) {
      return(NULL)
    }

    # Circular correlation
    # circ_corr <- circular_correlation(climate_phase, dengue_phase)

    # Mean phase lag (months climate leads dengue)
    phase_lag_months <- calculate_phase_lag(climate_phase, dengue_phase)
    
    # Correlation by lag
    phase_lag_round <- round(phase_lag_months)
    data_length <- length(climate_phase)
    if (phase_lag_round > 0) {
      # Positive lag: x2 is earlier than x1 by 'lag' positions
      # Start x1 from position (1 + lag), remove last 'lag' elements from x2
      dengue_phase_subset <- dengue_phase[(1 + phase_lag_round):data_length]
      climate_phase_subset <- climate_phase[1:(data_length - phase_lag_round)]
    } else if (phase_lag_round < 0) {
      # Negative lag: x1 is earlier than x2 by 'abs(lag)' positions
      # Start x2 from position (1 + abs(lag)), remove last 'abs(lag)' elements from x1
      dengue_phase_subset <- dengue_phase[1:(data_length + phase_lag_round)]
      climate_phase_subset <- climate_phase[(1 - phase_lag_round):data_length]
    } else {
      # No lag
      dengue_phase_subset <- dengue_phase
      climate_phase_subset <- climate_phase
    }
    
    # Calculate Spearman correlation
    correlation <- cor(dengue_phase_subset, climate_phase_subset, method = "spearman")

    # Phase coherence (how consistent is the phase relationship over time?)
    phase_diff <- climate_phase - dengue_phase
    phase_diff_wrapped <- ((phase_diff + pi) %% (2 * pi)) - pi
    phase_coherence <- 1 - (sd(phase_diff_wrapped, na.rm = TRUE) / pi)  # 1 = perfectly consistent, 0 = random

    tibble(
      idadmin1 = idadmin1_indonesia_unique[i],
      admin1 = province_name,
      climate_var = gsub("_scaled", "", clim_var),
      correlation = correlation,
      phase_lag_months = phase_lag_months,
      phase_coherence = phase_coherence
    )
  })

  bind_rows(var_results)
})

climate_dengue_phase_df <- bind_rows(climate_dengue_phase) %>%
  left_join(province_region_map %>% dplyr::select(idadmin1, region), by = "idadmin1")

# Summary statistics
print("Climate-Dengue Phase Relationships (wavelet-based):")
print("Positive phase lag = climate leads dengue; Negative = dengue leads climate")
print(climate_dengue_phase_df %>%
        group_by(climate_var) %>%
        summarise(
          mean_phase_lag = mean(phase_lag_months, na.rm = TRUE),
          sd_phase_lag = sd(phase_lag_months, na.rm = TRUE),
          mean_correlation = mean(correlation, na.rm = TRUE),
          mean_coherence = mean(phase_coherence, na.rm = TRUE)
        ))

# By region
print("Phase relationships by region:")
climate_dengue_phase_by_region <- climate_dengue_phase_df %>%
  group_by(region, climate_var) %>%
  summarise(
    mean_phase_lag = mean(phase_lag_months, na.rm = TRUE),
    mean_correlation = mean(correlation, na.rm = TRUE),
    mean_coherence = mean(phase_coherence, na.rm = TRUE),
    .groups = "drop"
  )
print(climate_dengue_phase_by_region)

# Identify provinces where climate-dengue relationship is strongest (highest coherence)
print("Provinces with strongest climate-dengue phase coherence:")
print(climate_dengue_phase_df %>%
        group_by(admin1, region) %>%
        summarise(mean_coherence = mean(phase_coherence, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(mean_coherence)) %>%
        head(10))

# Identify provinces where climate-dengue relationship is strongest (highest coherence)
# Only for positive phase lag (where climate leads dengue)
print("Provinces with strongest climate-dengue phase coherence:")
print(climate_dengue_phase_df %>%
        filter(round(phase_lag_months)>=0) %>%
        group_by(admin1, region) %>%
        summarise(mean_coherence = mean(phase_coherence, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(mean_coherence)) %>%
        head(10))

################################################################################
# PHASE COHERENCE THRESHOLD SENSITIVITY ANALYSIS
################################################################################

# Test sensitivity to phase coherence threshold choice (0.70-0.90)

cat("\n=== Phase Coherence Threshold Sensitivity ===\n\n")

# Note: DLNM significance data will be available after DLNM section runs
# For now, just count qualifying provinces by threshold

coherence_thresholds_tested <- seq(0.70, 0.90, by = 0.05)

threshold_sensitivity_prelim <- tibble()

for (threshold_val in coherence_thresholds_tested) {

  for (clim_var in c("precipitation_era", "temperature_era")) {

    # Count qualifying provinces at this threshold
    n_qual <- climate_dengue_phase_df %>%
      filter(climate_var == clim_var,
             phase_coherence >= threshold_val,
             round(phase_lag_months) >= 0) %>%
      nrow()

    threshold_sensitivity_prelim <- bind_rows(
      threshold_sensitivity_prelim,
      tibble(
        threshold = threshold_val,
        climate_var = clim_var,
        n_qualifying = n_qual
      )
    )
  }
}

cat("Qualifying provinces by threshold:\n")
print(threshold_sensitivity_prelim %>%
        pivot_wider(names_from = climate_var, values_from = n_qualifying))

cat("\nNote: DLNM significance counts will be added in supplementary materials\n")
cat("after DLNM analysis completes.\n\n")

climate_dengue_mean_coherence_df <- climate_dengue_phase_df %>%
  filter(round(phase_lag_months)>=0) %>% 
  group_by(idadmin1, admin1, region) %>%
  summarise(mean_coherence = mean(phase_coherence, na.rm = TRUE), .groups = "drop")

################################################################################
# SECTION: DLNM ANALYSIS
################################################################################

#### Prepare data for DLNM
dengue_data_dlnm <- dengue_data_admin1 %>%
  dplyr::select(region,idadmin1,admin1,date,pop,cases,incidence=incidence_rate) %>%
  left_join(prec_admin1_era_2009_2024 %>% dplyr::select(region,idadmin1,admin1,date,precipitation)) %>%
  left_join(temp_admin1_era_2009_2024 %>% dplyr::select(region,idadmin1,admin1,date,temperature)) %>%
  left_join(humi_admin1_era_2009_2024 %>% dplyr::select(region,idadmin1,admin1,date,rel_humidity)) %>%
  left_join(oni_2009_2024 %>% dplyr::select(date,oni)) %>%
  left_join(iod_2009_2024 %>% dplyr::select(date,dmi)) %>%
  arrange(idadmin1,date) %>%
  mutate(year=year(date),month=month(date))

#### Get unique province IDs for model runs
idadmin1_unique <- unique(dengue_data_dlnm$idadmin1)

################################################################################
# DLNM: Final Model Runs (Spline + Negative Binomial)
################################################################################

# Parameters for final model (selected based on comparison)
parameters_univariate <- expand.grid(
  idadmin1 = idadmin1_unique,
  lag_value_long = 6,
  lag_value_short = 4,
  stringsAsFactors = FALSE
)

# Run final models for all provinces
dlnm_output_list <- list()
contour_plots_list <- list()
lag_plots_list <- list()
significance_table_list <- list()

for (i in seq_len(nrow(parameters_univariate))) {
  
  parameters_model_run <- parameters_univariate[i, ]
  cb_fun <- parameters_model_run$cb_fun
  likelihood <- parameters_model_run$likelihood
  lag_value_long <- parameters_model_run$lag_value_long
  lag_value_short <- parameters_model_run$lag_value_short
  idadmin1_run <- parameters_model_run$idadmin1
  
  data_filter_dlnm <- dengue_data_dlnm %>% filter(idadmin1 == idadmin1_run)
  data_filter_dlnm$time_index <- 1:nrow(data_filter_dlnm)
  
  dlnm_output_list[[i]] <- dlnm_run(
    data_filter_dlnm, lag_value_long, lag_value_short
  )
  contour_plots_list[[i]] <- dlnm_contour_visualise(
    dlnm_output_list[[i]], data_filter_dlnm
  )
  lag_plots_list[[i]] <- dlnm_lag_visualise(dlnm_output_list[[i]])
  significance_table_list[[i]] <- dlnm_significance_lag(
    dlnm_output_list[[i]], data_filter_dlnm
  )
}

cumul_RR_highest_maximum <- tibble()
for (i in seq_len(length(idadmin1_unique))){
  for (cv in c("ONI","DMI","Precipitation","Temperature","Relative humidity")){
    for (pct in c(10,90)){
      
      cumul_RR_df <- dlnm_output_list[[i]]$cumulative_RR %>% 
        filter(clim_var==cv,percentile==pct) %>% 
        mutate(lag_highest=lag)
      
      cumul_RR_highest_maximum <- bind_rows(cumul_RR_highest_maximum,find_optimal_lag(cumul_RR_df))
      
    }
  }
}

cumul_RR_highest_first <- tibble()
for (i in seq_len(length(idadmin1_unique))){
  for (cv in c("ONI","DMI","Precipitation","Temperature","Relative humidity")){
    for (pct in c(10,90)){
      
      cumul_RR_df <- dlnm_output_list[[i]]$cumulative_RR %>% 
        filter(clim_var==cv,percentile==pct) %>% 
        mutate(lag_highest=lag)
      
      cumul_RR_highest_first <- bind_rows(cumul_RR_highest_first,find_optimal_lag_first(cumul_RR_df))
      
    }
  }
}

cumul_RR_all <- tibble()
for (i in seq_len(length(idadmin1_unique))){
  for (cv in c("ONI","DMI","Precipitation","Temperature","Relative humidity")){
    for (pct in c(10,90)){
      
      cumul_RR_all <- bind_rows(cumul_RR_all,dlnm_output_list[[i]]$cumulative_RR %>% 
                                  filter(clim_var==cv,percentile==pct))
      
    }
  }
}

# Create phase_heatmap_data for subsequent filtering
phase_heatmap_data <- climate_dengue_phase_df %>%
  mutate(
    climate_var = case_when(
      climate_var == "precipitation_era" ~ "Precipitation",
      climate_var == "temperature_era" ~ "Temperature",
      climate_var == "rel_humidity_era" ~ "Relative humidity",
      TRUE ~ climate_var
    )
  ) %>%
  left_join(admin1_EN, by = "idadmin1")

phase_heatmap_data_over85 <- phase_heatmap_data %>%
  filter(phase_coherence >= 0.85 & climate_var != "Relative humidity") %>%
  mutate(lag=round(phase_lag_months)) %>%
  filter(lag >= 0) %>% 
  mutate(clim_var=climate_var)

phase_heatmap_data_pct90_over85_RR <- phase_heatmap_data_over85 %>% left_join(cumul_RR_all %>% filter(percentile==90)) %>%
  mutate(significant = if_else(
    (cumulRR_lo > 1 & cumulRR_hi > 1) | (cumulRR_lo < 1 & cumulRR_hi < 1),
    1,
    0
  ))
phase_heatmap_data_pct10_over85_RR <- phase_heatmap_data_over85 %>% left_join(cumul_RR_all %>% filter(percentile==10)) %>%
  mutate(significant = if_else(
    (cumulRR_lo > 1 & cumulRR_hi > 1) | (cumulRR_lo < 1 & cumulRR_hi < 1),
    1,
    0
  ))

phase_heatmap_data_pct90_over85_RR_summary <- phase_heatmap_data_pct90_over85_RR %>%
  group_by(clim_var) %>%
  summarise(
    n = n(),
    n_significant = sum(significant == 1),
    pct_significant = round(n_significant / n * 100, 1)
  ) %>% 
  mutate(threshold=0.85) %>% dplyr::select(threshold,everything())

# Sensitivity analysis for other threshold of phase coherence
# 70
phase_heatmap_data_over70 <- phase_heatmap_data %>%
  filter(phase_coherence >= 0.70 & climate_var != "Relative humidity") %>%
  mutate(lag=round(phase_lag_months)) %>%
  filter(lag >= 0) %>% 
  mutate(clim_var=climate_var)

phase_heatmap_data_pct90_over70_RR <- phase_heatmap_data_over70 %>% left_join(cumul_RR_all %>% filter(percentile==90)) %>%
  mutate(significant = if_else(
    (cumulRR_lo > 1 & cumulRR_hi > 1) | (cumulRR_lo < 1 & cumulRR_hi < 1),
    1,
    0
  ))
phase_heatmap_data_pct10_over70_RR <- phase_heatmap_data_over70 %>% left_join(cumul_RR_all %>% filter(percentile==10)) %>%
  mutate(significant = if_else(
    (cumulRR_lo > 1 & cumulRR_hi > 1) | (cumulRR_lo < 1 & cumulRR_hi < 1),
    1,
    0
  ))

phase_heatmap_data_pct90_over70_RR_summary <- phase_heatmap_data_pct90_over70_RR %>%
  group_by(clim_var) %>%
  summarise(
    n = n(),
    n_significant = sum(significant == 1),
    pct_significant = round(n_significant / n * 100, 1)
  ) %>% 
  mutate(threshold=0.70) %>% dplyr::select(threshold,everything())

# 75
phase_heatmap_data_over75 <- phase_heatmap_data %>%
  filter(phase_coherence >= 0.75 & climate_var != "Relative humidity") %>%
  mutate(lag=round(phase_lag_months)) %>%
  filter(lag >= 0) %>% 
  mutate(clim_var=climate_var)

phase_heatmap_data_pct90_over75_RR <- phase_heatmap_data_over75 %>% left_join(cumul_RR_all %>% filter(percentile==90)) %>%
  mutate(significant = if_else(
    (cumulRR_lo > 1 & cumulRR_hi > 1) | (cumulRR_lo < 1 & cumulRR_hi < 1),
    1,
    0
  ))
phase_heatmap_data_pct10_over75_RR <- phase_heatmap_data_over75 %>% left_join(cumul_RR_all %>% filter(percentile==10)) %>%
  mutate(significant = if_else(
    (cumulRR_lo > 1 & cumulRR_hi > 1) | (cumulRR_lo < 1 & cumulRR_hi < 1),
    1,
    0
  ))

phase_heatmap_data_pct90_over75_RR_summary <- phase_heatmap_data_pct90_over75_RR %>%
  group_by(clim_var) %>%
  summarise(
    n = n(),
    n_significant = sum(significant == 1),
    pct_significant = round(n_significant / n * 100, 1)
  ) %>% 
  mutate(threshold=0.75) %>% dplyr::select(threshold,everything())

# 80
phase_heatmap_data_over80 <- phase_heatmap_data %>%
  filter(phase_coherence >= 0.80 & climate_var != "Relative humidity") %>%
  mutate(lag=round(phase_lag_months)) %>%
  filter(lag >= 0) %>% 
  mutate(clim_var=climate_var)

phase_heatmap_data_pct90_over80_RR <- phase_heatmap_data_over80 %>% left_join(cumul_RR_all %>% filter(percentile==90)) %>%
  mutate(significant = if_else(
    (cumulRR_lo > 1 & cumulRR_hi > 1) | (cumulRR_lo < 1 & cumulRR_hi < 1),
    1,
    0
  ))
phase_heatmap_data_pct10_over80_RR <- phase_heatmap_data_over80 %>% left_join(cumul_RR_all %>% filter(percentile==10)) %>%
  mutate(significant = if_else(
    (cumulRR_lo > 1 & cumulRR_hi > 1) | (cumulRR_lo < 1 & cumulRR_hi < 1),
    1,
    0
  ))

phase_heatmap_data_pct90_over80_RR_summary <- phase_heatmap_data_pct90_over80_RR %>%
  group_by(clim_var) %>%
  summarise(
    n = n(),
    n_significant = sum(significant == 1),
    pct_significant = round(n_significant / n * 100, 1)
  ) %>% 
  mutate(threshold=0.80) %>% dplyr::select(threshold,everything())

# 90
phase_heatmap_data_over90 <- phase_heatmap_data %>%
  filter(phase_coherence >= 0.90 & climate_var != "Relative humidity") %>%
  mutate(lag=round(phase_lag_months)) %>%
  filter(lag >= 0) %>% 
  mutate(clim_var=climate_var)

phase_heatmap_data_pct90_over90_RR <- phase_heatmap_data_over90 %>% left_join(cumul_RR_all %>% filter(percentile==90)) %>%
  mutate(significant = if_else(
    (cumulRR_lo > 1 & cumulRR_hi > 1) | (cumulRR_lo < 1 & cumulRR_hi < 1),
    1,
    0
  ))
phase_heatmap_data_pct10_over90_RR <- phase_heatmap_data_over90 %>% left_join(cumul_RR_all %>% filter(percentile==10)) %>%
  mutate(significant = if_else(
    (cumulRR_lo > 1 & cumulRR_hi > 1) | (cumulRR_lo < 1 & cumulRR_hi < 1),
    1,
    0
  ))

phase_heatmap_data_pct90_over90_RR_summary <- phase_heatmap_data_pct90_over90_RR %>%
  group_by(clim_var) %>%
  summarise(
    n = n(),
    n_significant = sum(significant == 1),
    pct_significant = round(n_significant / n * 100, 1)
  ) %>% 
  mutate(threshold=0.90) %>% dplyr::select(threshold,everything())

phase_significance_threshold_sensitivity <- 
  bind_rows(phase_heatmap_data_pct90_over70_RR_summary,
            phase_heatmap_data_pct90_over75_RR_summary,
            phase_heatmap_data_pct90_over80_RR_summary,
            phase_heatmap_data_pct90_over85_RR_summary,
            phase_heatmap_data_pct90_over90_RR_summary)

################################################################################
# DLNM: Extract Significance Tables
################################################################################

# Initialize lists for positive and negative effects
significance_positive_oni_list <- list()
significance_positive_dmi_list <- list()
significance_positive_pre_list <- list()
significance_positive_tem_list <- list()
significance_positive_hum_list <- list()
significance_negative_oni_list <- list()
significance_negative_dmi_list <- list()
significance_negative_pre_list <- list()
significance_negative_tem_list <- list()
significance_negative_hum_list <- list()

for (i in seq_len(nrow(parameters_univariate))) {
  significance_positive_oni_list[[i]] <- significance_table_list[[i]]$significance_pos_oni
  significance_positive_dmi_list[[i]] <- significance_table_list[[i]]$significance_pos_dmi
  significance_positive_pre_list[[i]] <- significance_table_list[[i]]$significance_pos_pre
  significance_positive_tem_list[[i]] <- significance_table_list[[i]]$significance_pos_tem
  significance_positive_hum_list[[i]] <- significance_table_list[[i]]$significance_pos_hum

  significance_negative_oni_list[[i]] <- significance_table_list[[i]]$significance_neg_oni
  significance_negative_dmi_list[[i]] <- significance_table_list[[i]]$significance_neg_dmi
  significance_negative_pre_list[[i]] <- significance_table_list[[i]]$significance_neg_pre
  significance_negative_tem_list[[i]] <- significance_table_list[[i]]$significance_neg_tem
  significance_negative_hum_list[[i]] <- significance_table_list[[i]]$significance_neg_hum
}

# Combine positive effects
significance_positive_oni <- bind_rows(significance_positive_oni_list) %>% left_join(admin1_national_EN)
significance_positive_dmi <- bind_rows(significance_positive_dmi_list) %>% left_join(admin1_national_EN)
significance_positive_pre <- bind_rows(significance_positive_pre_list) %>% left_join(admin1_national_EN)
significance_positive_tem <- bind_rows(significance_positive_tem_list) %>% left_join(admin1_national_EN)
significance_positive_hum <- bind_rows(significance_positive_hum_list) %>% left_join(admin1_national_EN)

significance_positive <- bind_rows(
  significance_positive_oni, significance_positive_dmi,
  significance_positive_pre, significance_positive_tem,
  significance_positive_hum
) %>%
  mutate(clim_var = factor(clim_var, levels = c(
    "ONI", "DMI", "Precipitation", "Temperature", "Relative humidity"
  )))

# Combine negative effects
significance_negative_oni <- bind_rows(significance_negative_oni_list) %>% left_join(admin1_national_EN)
significance_negative_dmi <- bind_rows(significance_negative_dmi_list) %>% left_join(admin1_national_EN)
significance_negative_pre <- bind_rows(significance_negative_pre_list) %>% left_join(admin1_national_EN)
significance_negative_tem <- bind_rows(significance_negative_tem_list) %>% left_join(admin1_national_EN)
significance_negative_hum <- bind_rows(significance_negative_hum_list) %>% left_join(admin1_national_EN)

significance_negative <- bind_rows(
  significance_negative_oni, significance_negative_dmi,
  significance_negative_pre, significance_negative_tem,
  significance_negative_hum
) %>%
  mutate(clim_var = factor(clim_var, levels = c(
    "ONI", "DMI", "Precipitation", "Temperature", "Relative humidity"
  )))
