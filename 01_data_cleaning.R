# Clean GPS tracks and export the dataset to be used in all following analyses

# Load libraries
library(ctmm)
library(lubridate)
library(move2)
library(tidyverse)

# Two focal females used in the study
focal_animals <- c("704870A", "704871A")

# Create the output directory needed by the next script
dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

# Load GPS data
gps_data <- read_rds("data/raw/gps_maned_wolf.rds")

# Standardize fields needed by ctmm and convert as telemetry objects
data_telemetry <- gps_data |>
  mutate(
    timestamp = ymd_hms(timestamp, quiet = TRUE),
    timestamp_local = with_tz(timestamp, "America/Sao_Paulo"),
    gps_latitude = as.numeric(gps_latitude),
    gps_longitude = as.numeric(gps_longitude)
  ) |>
  rename(
    location.long = gps_longitude,
    location.lat = gps_latitude,
    individual.local.identifier = animal_id
  ) |>
  filter(
    individual.local.identifier %in% focal_animals,
    !is.na(location.long),
    !is.na(location.lat),
    !is.na(timestamp),
    year(timestamp_local) == 2024
  ) |>
  select(-timestamp_local) |>
  ctmm::as.telemetry()

# Animals' name
animals <- names(data_telemetry)

# Load the accessory function that applies sequential distance and speed filters
source("accessory_functions/01_function_ctmm_remove_outliers.R")

telemetry_no_out <- ctmm_remove_outliers(
  x = data_telemetry,
  individual_names = animals,
  dist = 10000,
  speed = 0.5
)

# Remove the very first days when residual points still fell inside the house and keep only the year for this study
telemetry_trimmed <- list()

for (animal in names(telemetry_no_out)) {
  start_time <- min(mt_time(telemetry_no_out[[animal]]), na.rm = TRUE)

  telemetry_trimmed[[animal]] <- telemetry_no_out[[animal]] |>
    filter(
      mt_time(telemetry_no_out[[animal]]) > (start_time + days(3))
    )
}

# Keep only standard columns that will be reused by the next analysis scripts
class_target <- c(
  "character",
  "numeric",
  "logical",
  "factor",
  "Date",
  "integer64",
  "POSIXct",
  "POSIXlt"
)

clean_telemetry_final <- list()

for (animal in names(telemetry_trimmed)) {
  clean_telemetry_final[[animal]] <- telemetry_trimmed[[animal]] |>
    as.data.frame() |>
    select(where(~ any(class(.) %in% class_target)))
}

# Save the single CSV consumed by the home-range and TEHS workflows
clean_telemetry_final |>
  bind_rows() |>
  write_csv("data/processed/gps_no_outliers_final.csv")
