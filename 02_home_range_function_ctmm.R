# Home range objects for females and save model outputs

# Load libraries
library(tidyverse)

# Source accessory function
source("accessory_functions/02_function_ctmm_seasonal.R")

# Load the cleaned dataset created in the previous step
data <- readr::read_csv(
  "data/processed/gps_no_outliers_final.csv"
) |>
  dplyr::mutate(
    timestamp = lubridate::ymd_hms(timestamp, tz = "UTC"),
    timestamp_local = lubridate::with_tz(timestamp, "America/Sao_Paulo")
  ) |>
  dplyr::rename(individual_local_identifier = track) |>
  dplyr::relocate(individual_local_identifier, .before = dplyr::everything()) |>
  dplyr::filter(lubridate::year(timestamp_local) == 2024) |>
  dplyr::mutate(
    month_num = lubridate::month(timestamp_local),
    calendar_year = lubridate::year(timestamp_local),
    season_year = dplyr::if_else(
      month_num <= 3,
      calendar_year - 1,
      calendar_year
    ),
    season = as.factor(dplyr::if_else(
      month_num >= 10 | month_num <= 3,
      "Wet",
      "Dry"
    )),
    period = paste(season_year, season, sep = "_")
  ) |>
  dplyr::select(-timestamp_local)

# Save one telemetry/model/AKDE set per biological period
process_period_data(data)
