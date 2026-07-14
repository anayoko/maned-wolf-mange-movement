# Time-Explicit Habitat Selection Model - data preparation

# Load libraries
library(terra)
library(sf)
library(tidyverse)

options(scipen = 999)
terra::terraOptions(progress = 0)

# Paths and parameters ----
raster_path <- "data/raster/mapbiomas_10m_sp_2023-0000000000-0000000000.tif"
raster_tag <- basename(raster_path) |> str_remove("\\.tif$")

tehs_dir <- "output/TEHS"
mapbiomas_dir <- file.path(tehs_dir, "mapbiomas")
prep_dir <- file.path(tehs_dir, "prep")
potential_steps_dir <- file.path(prep_dir, "potential_steps")
prep_complete_dir <- file.path(prep_dir, "prep_jags_complete")

buffer_m <- 30
miss_factor <- 1.5
cardinal_names <- c("west", "east", "north", "south")

walk(
  c(prep_dir, potential_steps_dir, prep_complete_dir),
  ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE)
)

if (!file.exists(raster_path)) {
  stop(
    "Missing MapBiomas LULC raster: ", raster_path,
    "\nDownload the 2023, 10-m LULC raster for São Paulo from ",
    "https://brasil.mapbiomas.org/en/ and save it at this path before running this script."
  )
}

# Input files ----
steps_path <- file.path(
  mapbiomas_dir,
  str_glue("step_lulc_prop_common_{raster_tag}.csv")
)
lulc_columns_path <- file.path(
  mapbiomas_dir,
  str_glue("step_lulc_columns_kept_{raster_tag}.csv")
)
lulc_usage_path <- file.path(
  mapbiomas_dir,
  str_glue("step_lulc_overall_use_threshold_{raster_tag}.csv")
)

if (!file.exists(steps_path)) {
  stop("Missing input file: ", steps_path)
}
if (!file.exists(lulc_columns_path)) {
  stop("Missing input file: ", lulc_columns_path)
}
if (!file.exists(lulc_usage_path)) {
  stop("Missing input file: ", lulc_usage_path)
}

# MapBiomas legend ----
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

# Read inputs ----
mapbiomas_raster <- rast(raster_path)
raster_lulc_col <- names(mapbiomas_raster)[1]

steps_common <- read_csv(steps_path, show_col_types = FALSE)
lulc_columns <- read_csv(lulc_columns_path, show_col_types = FALSE)
lulc_usage <- read_csv(lulc_usage_path, show_col_types = FALSE)

required_step_cols <- c(
  "individual",
  "step_id",
  "time_min",
  "dist_m",
  "x_sta",
  "y_sta",
  "x_end",
  "y_end"
)
missing_required <- setdiff(required_step_cols, names(steps_common))
if (length(missing_required) > 0) {
  stop(
    "Missing required columns in step input: ",
    paste(missing_required, collapse = ", ")
  )
}

lulc_cols_kept <- lulc_columns$lulc_col
if (!all(lulc_cols_kept %in% names(steps_common))) {
  stop("Not all LULC columns from lulc_columns file exist in step input.")
}

# Define model covariates (baseline = most used class) ----
baseline_group <- lulc_usage |>
  filter(keep_class) |>
  arrange(desc(prop_use)) |>
  slice(1) |>
  pull(class_group)

baseline_col <- lulc_columns |>
  filter(class_group == baseline_group) |>
  pull(lulc_col)

model_covariates <- lulc_columns |>
  filter(lulc_col != baseline_col) |>
  pull(lulc_col)

covariate_config <- lulc_columns |>
  mutate(
    raster_tag = raster_tag,
    baseline_group = baseline_group,
    baseline_col = baseline_col,
    use_in_model = lulc_col %in% model_covariates,
    covariate_index = cumsum(use_in_model)
  )

write_csv(
  covariate_config,
  file.path(prep_dir, str_glue("tehs_covariate_config_{raster_tag}.csv"))
)

# Build time variables used by JAGS time model ----
fix_interval_min <- median(steps_common$time_min, na.rm = TRUE)
miss_threshold_min <- fix_interval_min * miss_factor

steps_tehs <- steps_common |>
  mutate(
    time1 = time_min,
    dist1 = dist_m,
    Miss = if_else(time1 <= miss_threshold_min, 0, 1)
  )

animals <- steps_tehs |>
  distinct(individual) |>
  pull(individual)

# Prepare per-animal files ----
for (animal in animals) {
  cat("\\n--------------------------------------------\\n")
  cat("Preparing TEHS files for animal: ", animal, "\\n", sep = "")

  steps_animal <- steps_tehs |>
    filter(individual == animal) |>
    arrange(time_sta) |>
    mutate(
      x_end1 = x_sta - dist1,
      y_end1 = y_sta,
      x_end2 = x_sta + dist1,
      y_end2 = y_sta,
      x_end3 = x_sta,
      y_end3 = y_sta + dist1,
      x_end4 = x_sta,
      y_end4 = y_sta - dist1
    )

  if (nrow(steps_animal) == 0) {
    warning("No steps for animal: ", animal)
    next
  }

  # Crop the raster once to the full extent covered by the observed and potential steps for this animal
  animal_x <- c(
    steps_animal$x_sta,
    steps_animal$x_end,
    steps_animal$x_end1,
    steps_animal$x_end2,
    steps_animal$x_end3,
    steps_animal$x_end4
  )

  animal_y <- c(
    steps_animal$y_sta,
    steps_animal$y_end,
    steps_animal$y_end1,
    steps_animal$y_end2,
    steps_animal$y_end3,
    steps_animal$y_end4
  )

  animal_extent_sf <- sf::st_as_sf(
    tibble(
      geometry = sf::st_sfc(
        sf::st_as_sfc(
          sf::st_bbox(
            c(
              xmin = min(animal_x, na.rm = TRUE) - buffer_m,
              xmax = max(animal_x, na.rm = TRUE) + buffer_m,
              ymin = min(animal_y, na.rm = TRUE) - buffer_m,
              ymax = max(animal_y, na.rm = TRUE) + buffer_m
            ),
            crs = sf::st_crs(5880)
          )
        )
      )
    )
  ) |>
    sf::st_transform(crs = terra::crs(mapbiomas_raster))

  mapbiomas_crop_animal <- terra::crop(
    mapbiomas_raster,
    terra::ext(terra::vect(animal_extent_sf)),
    snap = "out"
  )

  write_csv(
    steps_animal,
    file.path(
      potential_steps_dir,
      str_glue("potential_steps_coordinates_{animal}_{raster_tag}.csv")
    )
  )

  # Build all four candidate directions together so the raster is extracted only once per animal
  candidate_steps <- purrr::map_dfr(seq_along(cardinal_names), \(j) {
    x_end_col <- str_glue("x_end{j}")
    y_end_col <- str_glue("y_end{j}")

    steps_animal |>
      transmute(
        step_id,
        path_idx = j,
        x_sta,
        y_sta,
        x_end = .data[[x_end_col]],
        y_end = .data[[y_end_col]]
      )
  })

  candidate_lines_geom <- lapply(seq_len(nrow(candidate_steps)), function(i) {
    st_linestring(
      matrix(
        c(
          candidate_steps$x_sta[i],
          candidate_steps$y_sta[i],
          candidate_steps$x_end[i],
          candidate_steps$y_end[i]
        ),
        ncol = 2,
        byrow = TRUE
      )
    )
  })

  candidate_buffers <- st_sf(
    row_id = seq_len(nrow(candidate_steps)),
    step_id = candidate_steps$step_id,
    path_idx = candidate_steps$path_idx,
    geometry = st_sfc(candidate_lines_geom, crs = 5880)
  ) |>
    st_buffer(dist = buffer_m)

  candidate_buffers_vect <- vect(candidate_buffers)
  if (!same.crs(mapbiomas_crop_animal, candidate_buffers_vect)) {
    candidate_buffers_vect <- project(
      candidate_buffers_vect,
      mapbiomas_crop_animal
    )
  }

  extracted <- terra::extract(
    mapbiomas_crop_animal[[1]],
    candidate_buffers_vect,
    ID = TRUE
  ) |>
    as_tibble() |>
    rename(row_id = ID) |>
    mutate(class_code = as.integer(.data[[raster_lulc_col]])) |>
    select(row_id, class_code) |>
    filter(!is.na(class_code)) |>
    left_join(
      candidate_buffers |>
        st_drop_geometry() |>
        select(row_id, step_id, path_idx),
      by = "row_id"
    )

  if (nrow(extracted) == 0) {
    stop("No LULC extracted for animal ", animal)
  }

  candidate_covariates <- extracted |>
    left_join(lulc_legend, by = "class_code") |>
    mutate(
      class_name = coalesce(class_name, paste0("Class ", class_code)),
      class_group = case_when(
        class_code %in% c(19L, 36L) ~ "Crops",
        TRUE ~ class_name
      )
    ) |>
    count(step_id, path_idx, class_group, name = "n_pixels") |>
    group_by(step_id, path_idx) |>
    mutate(prop_lulc = n_pixels / sum(n_pixels)) |>
    ungroup() |>
    left_join(lulc_columns, by = "class_group") |>
    filter(!is.na(lulc_col)) |>
    group_by(step_id, path_idx, lulc_col) |>
    summarise(prop_lulc = sum(prop_lulc), .groups = "drop") |>
    mutate(path_col = str_glue("{lulc_col}_{path_idx}")) |>
    select(step_id, path_col, prop_lulc) |>
    pivot_wider(
      names_from = path_col,
      values_from = prop_lulc,
      values_fill = 0
    )

  missing_candidate_cols <- setdiff(
    as.vector(outer(lulc_cols_kept, 1:4, paste, sep = "_")),
    names(candidate_covariates)
  )

  if (length(missing_candidate_cols) > 0) {
    candidate_covariates[missing_candidate_cols] <- 0
  }

  # Build full tSSF file (realized + 4 cardinal alternatives)
  realized_covariates <- steps_animal |>
    select(step_id, time1, dist1, Miss, all_of(lulc_cols_kept)) |>
    rename_with(~ str_glue("{.x}_0"), all_of(lulc_cols_kept))

  tssf_data_full <- realized_covariates |>
    left_join(candidate_covariates, by = "step_id")

  covariate_cols_all_paths <- as.vector(
    outer(lulc_cols_kept, 0:4, paste, sep = "_")
  )

  tssf_data_full <- tssf_data_full |>
    mutate(across(all_of(covariate_cols_all_paths), ~ replace_na(.x, 0)))

  # Remove the baseline class across all alternatives and save only the table used by the next model-fitting script
  baseline_drop_cols <- str_glue("{baseline_col}_{0:4}")

  tssf_data_model <- tssf_data_full |>
    select(-any_of(baseline_drop_cols))

  write_csv(
    tssf_data_model,
    file.path(
      prep_complete_dir,
      str_glue("tSSF_data_{animal}_{raster_tag}.csv")
    )
  )
}

cat(
  "\\nTEHS prep finished.\\n",
  "Raster tag: ",
  raster_tag,
  "\\n",
  "Output dir: ",
  prep_dir,
  "\\n",
  "Baseline covariate: ",
  baseline_col,
  "\\n",
  "Model covariates: ",
  paste(model_covariates, collapse = ", "),
  "\\n",
  sep = ""
)
