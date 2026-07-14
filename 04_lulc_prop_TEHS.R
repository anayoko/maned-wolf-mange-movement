# Calculate step-level LULC proportions for maned wolf data

# Libraries ----
library(lubridate)
library(terra)
library(sf)
library(tidyverse)

options(scipen = 999) # avoid scientific notation
terra::terraOptions(progress = 0)

# Inputs and parameters ----
gps_path <- "data/processed/gps_no_outliers_final.csv"
raster_path <- "data/raster/mapbiomas_10m_sp_2023-0000000000-0000000000.tif"
output_dir <- "output/TEHS/mapbiomas"
raster_tag <- basename(raster_path) |> str_remove("\\.tif$")
focal_animals <- c("704870A", "704871A")

# Define parameters to be used in the code
buffer_m <- 30
speed_quantile <- 0.95
max_step_minutes <- 720
min_step_distance_m <- 1
min_lulc_use_prop <- 0.04

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(raster_path)) {
  stop(
    "Missing MapBiomas LULC raster: ", raster_path,
    "\nDownload the 2023, 10-m LULC raster for São Paulo from ",
    "https://brasil.mapbiomas.org/en/ and save it at this path before running this script."
  )
}

# MapBiomas legend (Collection 10 - relevant classes) ----
lulc_legend <- tibble::tribble(
  ~class_code , ~class_name                  ,
   3L         , "Forest Formation"           ,
   4L         , "Savanna Formation"          ,
   5L         , "Mangrove"                   ,
   6L         , "Floodable Forest"           ,
   9L         , "Forest Plantation"          ,
  11L         , "Wetland"                    ,
  12L         , "Grassland"                  ,
  13L         , "Other Non Forest Formation" ,
  15L         , "Pasture"                    ,
  18L         , "Agriculture"                ,
  19L         , "Temporary Crops"            ,
  20L         , "Sugarcane"                  ,
  21L         , "Mosaic of Uses"             ,
  23L         , "Beach, Dune and Sand Spot"  ,
  24L         , "Urban Infrastructure"       ,
  25L         , "Other Non Vegetated Area"   ,
  26L         , "Water"                      ,
  29L         , "Rocky Outcrop"              ,
  30L         , "Mining"                     ,
  31L         , "Aquaculture"                ,
  33L         , "River, Lake and Ocean"      ,
  36L         , "Perennial Crops"            ,
  49L         , "Wooded Restinga"
)

# Read GPS data directly ----
# Read and standardize columns used in the workflow
gps <- read_csv(gps_path, show_col_types = FALSE) |>
  rename(individual = track) |>
  mutate(
    timestamp = ymd_hms(timestamp, tz = "UTC", quiet = TRUE),
  ) |>
  select(individual, timestamp, longitude, latitude) |>
  filter(
    individual %in% focal_animals,
    !is.na(individual),
    !is.na(timestamp),
    !is.na(longitude),
    !is.na(latitude)
  )

# Read raster directly ----
# Use the confirmed MapBiomas raster file without tile selection logic
mapbiomas_raster <- rast(raster_path)
lulc_col <- names(mapbiomas_raster)[1]

# Convert GPS to projected coordinates for movement metrics ----
# Keep original lon/lat and add projected x/y in meters
gps_sf <- gps |>
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  ) |>
  st_transform(5880)

gps_coords <- st_coordinates(gps_sf)

# Add the columns 5880 to calulcate steps in meters
gps_5880 <- gps |>
  mutate(
    x_5880 = gps_coords[, 1],
    y_5880 = gps_coords[, 2]
  )

# Build movement steps ----
# Create start/end points and movement metrics for each individual
steps <- gps_5880 |>
  arrange(individual, timestamp) |>
  group_by(individual) |>
  mutate(
    x_sta = x_5880,
    y_sta = y_5880,
    x_end = lead(x_5880),
    y_end = lead(y_5880),
    time_sta = timestamp,
    time_end = lead(timestamp),
    time_min = as.numeric(difftime(time_end, time_sta, units = "mins")),
    dist_m = sqrt((x_end - x_sta)^2 + (y_end - y_sta)^2),
    speed_mps = dist_m / (time_min * 60),
    month_num = month(time_sta),
    calendar_year = year(time_sta),
    season_year = if_else(
      month_num <= 3,
      calendar_year - 1L,
      calendar_year
    ),
    season = if_else(month_num >= 10 | month_num <= 3, "Wet", "Dry"),
    period = paste(season_year, season, sep = "_")
  ) |>
  ungroup() |>
  filter(
    !is.na(x_end),
    !is.na(y_end),
    !is.na(time_min),
    time_min > 0,
    !is.na(dist_m),
    !is.na(speed_mps)
  ) |>
  mutate(
    dist_m = if_else(dist_m == 0, 1, dist_m)
  )

# Filter movement outliers ----
# Remove very fast movements and very long temporal gaps
speed_limits <- steps |>
  group_by(individual) |>
  summarise(
    speed_limit = quantile(speed_mps, speed_quantile, na.rm = TRUE),
    .groups = "drop"
  )

steps_filter <- steps |>
  left_join(speed_limits, by = "individual") |>
  filter(
    speed_mps <= speed_limit,
    time_min <= max_step_minutes,
    dist_m > min_step_distance_m
  ) |>
  select(-speed_limit)

if (nrow(steps_filter) == 0) {
  stop("All movement steps were removed by filters. Review thresholds.")
}

# Periods by year and season
period_levels <- steps_filter |>
  group_by(period) |>
  summarise(start_date = min(time_sta), .groups = "drop") |>
  arrange(start_date) |>
  pull(period)

steps_filter <- steps_filter |>
  mutate(
    step_id = row_number(),
    period = factor(period, levels = period_levels, ordered = TRUE)
  )

# Create line buffers for each step ----
# Each step is converted to a 30 m buffer polygon
step_lines_geom <- lapply(seq_len(nrow(steps_filter)), function(i) {
  st_linestring(
    matrix(
      c(
        steps_filter$x_sta[i],
        steps_filter$y_sta[i],
        steps_filter$x_end[i],
        steps_filter$y_end[i]
      ),
      ncol = 2,
      byrow = TRUE
    )
  )
})

# Create buffer
step_buffers <- st_sf(
  step_id = steps_filter$step_id,
  geometry = st_sfc(step_lines_geom, crs = 5880)
) |>
  st_buffer(dist = buffer_m)

step_buffers_vect <- vect(step_buffers)
if (!same.crs(mapbiomas_raster, step_buffers_vect)) {
  step_buffers_vect <- project(step_buffers_vect, mapbiomas_raster)
}

# Crop the raster once to the buffered step extent before extracting pixels
mapbiomas_crop <- terra::crop(
  mapbiomas_raster,
  terra::ext(step_buffers_vect),
  snap = "out"
)

# Extract raster classes inside each step buffer ----
# Keep only non-missing class codes
extracted <- terra::extract(
  mapbiomas_crop[[1]],
  step_buffers_vect,
  ID = TRUE
) |>
  as_tibble() |>
  rename(step_id = ID) |>
  mutate(
    class_code = as.integer(.data[[lulc_col]])
  ) |>
  select(step_id, class_code) |>
  filter(!is.na(class_code))

if (nrow(extracted) == 0) {
  stop("No LULC pixels were extracted from step buffers.")
}

# Recode classes for the analysis setup ----
# Merge Temporary and Perennial Crops
step_lulc <- extracted |>
  left_join(lulc_legend, by = "class_code") |>
  mutate(
    class_name = coalesce(class_name, paste0("Class ", class_code)),
    class_group = case_when(
      class_code %in% c(19L, 36L) ~ "Crops",
      TRUE ~ class_name
    )
  )

# Calculate per-step proportions in long format ----
step_lulc_counts <- step_lulc |>
  count(step_id, class_group, name = "n_pixels")

step_lulc_long_all <- step_lulc_counts |>
  group_by(step_id) |>
  mutate(prop_lulc = n_pixels / sum(n_pixels)) |>
  ungroup()

# Calculate overall usage and apply threshold ----
# Keep only classes above the minimum usage threshold
lulc_usage_overall <- step_lulc_counts |>
  group_by(class_group) |>
  summarise(n_pixels = sum(n_pixels), .groups = "drop") |>
  mutate(
    prop_use = n_pixels / sum(n_pixels),
    keep_class = prop_use >= min_lulc_use_prop
  ) |>
  arrange(desc(prop_use))

kept_groups <- lulc_usage_overall |>
  filter(keep_class) |>
  pull(class_group)

if (length(kept_groups) == 0) {
  stop(
    "No LULC class reached the threshold of ",
    scales::percent(min_lulc_use_prop, accuracy = 0.1),
    "."
  )
}

step_lulc_long <- step_lulc_long_all |>
  filter(class_group %in% kept_groups)

# Build column names for wide output ----
# Convert class labels into safe snake_case column names
lulc_columns <- lulc_usage_overall |>
  filter(keep_class) |>
  mutate(
    lulc_col = paste0(
      "lulc_",
      gsub(
        "_+",
        "_",
        gsub(
          "^_|_$",
          "",
          gsub(
            "[^a-z0-9]+",
            "_",
            tolower(iconv(class_group, to = "ASCII//TRANSLIT"))
          )
        )
      )
    )
  ) |>
  select(class_group, lulc_col)

# Create final step-level dataset (wide format) ----
# One row per step with one column per kept class
step_lulc_wide <- step_lulc_long |>
  left_join(lulc_columns, by = "class_group") |>
  select(step_id, lulc_col, prop_lulc) |>
  pivot_wider(
    names_from = lulc_col,
    values_from = prop_lulc,
    values_fill = 0
  )

final_database <- steps_filter |>
  left_join(step_lulc_wide, by = "step_id") |>
  arrange(individual, time_sta)

lulc_col_names <- lulc_columns$lulc_col
if (length(lulc_col_names) > 0) {
  final_database <- final_database |>
    mutate(
      across(all_of(lulc_col_names), ~ replace_na(.x, 0))
    )
}

# Save only the files consumed by the next TEHS preparation step
write_csv(
  final_database,
  file.path(output_dir, str_glue("step_lulc_prop_common_{raster_tag}.csv"))
)

write_csv(
  lulc_usage_overall,
  file.path(
    output_dir,
    str_glue("step_lulc_overall_use_threshold_{raster_tag}.csv")
  )
)

write_csv(
  lulc_columns,
  file.path(output_dir, str_glue("step_lulc_columns_kept_{raster_tag}.csv"))
)

cat(
  "LULC step extraction finished.\n",
  "Output dir: ",
  output_dir,
  "\n",
  "Threshold: ",
  scales::percent(min_lulc_use_prop, accuracy = 0.1),
  "\n",
  "Groups kept: ",
  paste(kept_groups, collapse = ", "),
  "\n",
  sep = ""
)
