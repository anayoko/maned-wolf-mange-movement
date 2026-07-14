# Consolidate TEHS posterior draws into the two files consumed directly by the final figure script

# Load library
library(tidyverse)

options(scipen = 999)

raster_path <- "data/raster/mapbiomas_10m_sp_2023-0000000000-0000000000.tif"
raster_tag <- basename(raster_path) |> str_remove("\\.tif$")

tehs_dir <- "output/TEHS"
prep_dir <- file.path(tehs_dir, "prep")
models_dir <- file.path(tehs_dir, "models")
time_model_results_dir <- file.path(models_dir, "time_model_results")
selection_model_results_dir <- file.path(models_dir, "selection_model_results")
results_dir <- file.path(tehs_dir, "results")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Read the model covariate order so the posterior columns stay aligned with the TEHS design matrix
covariate_config_path <- file.path(
  prep_dir,
  str_glue("tehs_covariate_config_{raster_tag}.csv")
)

if (!file.exists(covariate_config_path)) {
  stop("Missing covariate config: ", covariate_config_path)
}

covariate_config <- read_csv(covariate_config_path, show_col_types = FALSE)
model_covariates <- covariate_config |>
  filter(use_in_model) |>
  arrange(covariate_index) |>
  pull(lulc_col)

selection_files <- list.files(
  path = selection_model_results_dir,
  pattern = str_glue("^sims.list_.*_{raster_tag}[.]csv$"),
  full.names = TRUE
)

if (length(selection_files) == 0) {
  stop(
    "No selection model output files found in: ",
    selection_model_results_dir
  )
}

animals <- basename(selection_files) |>
  str_remove("^sims[.]list_") |>
  str_remove(str_glue("_{raster_tag}[.]csv$"))

selection_sims_list <- list()
time_sims_list <- list()

for (animal in animals) {
  selection_sims_path <- file.path(
    selection_model_results_dir,
    str_glue("sims.list_{animal}_{raster_tag}.csv")
  )
  time_sims_path <- file.path(
    time_model_results_dir,
    str_glue("sims.list_{animal}_{raster_tag}.csv")
  )

  if (file.exists(selection_sims_path)) {
    selection_sims_list[[animal]] <- read_csv(
      selection_sims_path,
      show_col_types = FALSE
    )
  }

  if (file.exists(time_sims_path)) {
    time_sims_list[[animal]] <- read_csv(
      time_sims_path,
      show_col_types = FALSE
    )
  }
}

# Save one posterior table for the habitat-selection model
selection_posterior <- map2_dfr(
  selection_sims_list,
  names(selection_sims_list),
  function(df, animal) {
    beta_cols <- names(df) |>
      keep(
        ~ str_detect(.x, "^betas\\[[0-9]+\\]$") |
          str_detect(.x, "^betas[.][0-9]+$") |
          str_detect(.x, "^betas[0-9]+$")
      )

    if (length(beta_cols) > 0) {
      beta_index <- case_when(
        str_detect(beta_cols, "^betas\\[[0-9]+\\]$") ~ as.numeric(str_extract(
          beta_cols,
          "[0-9]+"
        )),
        str_detect(beta_cols, "^betas[.][0-9]+$") ~ as.numeric(str_extract(
          beta_cols,
          "[0-9]+$"
        )),
        TRUE ~ as.numeric(str_remove(beta_cols, "^betas"))
      )
      beta_cols <- beta_cols[order(beta_index)]
    }

    df |>
      select(all_of(beta_cols)) |>
      set_names(model_covariates[seq_len(length(beta_cols))]) |>
      mutate(individual = animal, .before = 1)
  }
)

write_csv(
  selection_posterior,
  file.path(results_dir, str_glue("TEHS_selection_model_{raster_tag}.csv"))
)

# Save one posterior table for the time model
time_posterior <- map2_dfr(
  time_sims_list,
  names(time_sims_list),
  function(df, animal) {
    if (!"b0" %in% names(df)) {
      return(tibble())
    }

    b_effect_cols <- names(df) |>
      keep(
        ~ str_detect(.x, "^b\\[[0-9]+\\]$") |
          str_detect(.x, "^b[0-9]+$") |
          str_detect(.x, "^b[.][0-9]+$")
      ) |>
      discard(~ .x == "b0")

    if (length(b_effect_cols) > 0) {
      b_effect_index <- case_when(
        str_detect(b_effect_cols, "^b\\[[0-9]+\\]$") ~ as.numeric(str_extract(
          b_effect_cols,
          "[0-9]+"
        )),
        str_detect(b_effect_cols, "^b[.][0-9]+$") ~ as.numeric(str_extract(
          b_effect_cols,
          "[0-9]+"
        )),
        TRUE ~ as.numeric(str_remove(b_effect_cols, "^b"))
      )
      b_effect_cols <- b_effect_cols[order(b_effect_index)]
    }

    b_draws <- tibble(intercept = df$b0)

    if (length(b_effect_cols) > 0) {
      b_effect_df <- df |>
        select(all_of(b_effect_cols))

      names(b_effect_df) <- model_covariates[seq_len(ncol(b_effect_df))]
      b_draws <- bind_cols(b_draws, b_effect_df)
    }

    b_draws |>
      mutate(individual = animal, .before = 1)
  }
)

write_csv(
  time_posterior,
  file.path(results_dir, str_glue("TEHS_time_model_{raster_tag}.csv"))
)
