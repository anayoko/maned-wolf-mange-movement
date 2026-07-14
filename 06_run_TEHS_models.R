# Time-Explicit Habitat Selection Model - run time model and TEHS model

.libPaths(c("../r_libs", "r_libs", .libPaths()))

# Load libraries
library(jagsUI)
library(tidyverse)

options(scipen = 999)

# Paths and parameters ----
raster_path <- "data/raster/mapbiomas_10m_sp_2023-0000000000-0000000000.tif"
raster_tag <- basename(raster_path) |> str_remove("\\.tif$")

tehs_dir <- "output/TEHS"
prep_dir <- file.path(tehs_dir, "prep")
potential_steps_dir <- file.path(prep_dir, "potential_steps")
prep_complete_dir <- file.path(prep_dir, "prep_jags_complete")

models_dir <- file.path(tehs_dir, "models")
time_model_results_dir <- file.path(models_dir, "time_model_results")
selection_model_results_dir <- file.path(models_dir, "selection_model_results")

walk(
  c(models_dir, time_model_results_dir, selection_model_results_dir),
  ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE)
)

# Model files ----
time_model_file <- "accessory_functions/03_resist_avg_jags.R"
selection_model_file <- "accessory_functions/04_jags_tssf.R"

if (!file.exists(time_model_file)) {
  stop("Missing JAGS model file: ", time_model_file)
}
if (!file.exists(selection_model_file)) {
  stop("Missing JAGS model file: ", selection_model_file)
}

# Read covariate config ----
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

if (length(model_covariates) == 0) {
  stop("No model covariates available after baseline removal.")
}

# Discover animals from prepared files ----
animals <- list.files(
  path = prep_complete_dir,
  pattern = str_glue("^tSSF_data_.*_{raster_tag}[.]csv$")
) |>
  discard(~ str_detect(.x, "^tSSF_data_full_")) |>
  str_remove("^tSSF_data_") |>
  str_remove(str_glue("_{raster_tag}[.]csv$"))

if (length(animals) == 0) {
  stop("No prepared tSSF files found in: ", prep_complete_dir)
}

# 1) Time model ----
jags_time_results <- list()

for (animal in animals) {
  cat("\\n--------------------------------------------\\n")
  cat("Running TIME model for animal: ", animal, "\\n", sep = "")

  potential_path <- file.path(
    potential_steps_dir,
    str_glue("potential_steps_coordinates_{animal}_{raster_tag}.csv")
  )
  time_output_path <- file.path(
    time_model_results_dir,
    str_glue("sims.list_{animal}_{raster_tag}.csv")
  )

  if (!file.exists(potential_path)) {
    warning("Missing potential steps file for animal: ", animal)
    next
  }

  if (file.exists(time_output_path)) {
    message(
      "Time model output already exists for ",
      animal,
      ". Skipping refit."
    )
    next
  }

  time_data <- read_csv(potential_path, show_col_types = FALSE)

  required_time_cols <- c("time1", "dist1", "Miss", model_covariates)
  missing_time_cols <- setdiff(required_time_cols, names(time_data))
  if (length(missing_time_cols) > 0) {
    warning(
      "Skipping animal ",
      animal,
      " due to missing columns: ",
      paste(missing_time_cols, collapse = ", ")
    )
    next
  }

  dat_jags <- list(
    delta.time = time_data$time1,
    dist = time_data$dist1,
    nobs = nrow(time_data),
    ncov = length(model_covariates),
    Miss = time_data$Miss,
    xmat = data.matrix(time_data[, model_covariates])
  )

  set.seed(2)

  params <- c("b0", "b", "g0", "g1")

  n.iter <- 5000
  n.thin <- 10
  n.burnin <- n.iter / 2
  n.chains <- 3

  time_fit <- jags(
    model.file = time_model_file,
    parameters.to.save = params,
    data = dat_jags,
    n.chains = n.chains,
    n.burnin = n.burnin,
    n.iter = n.iter,
    n.thin = n.thin,
    DIC = FALSE
  )

  jags_time_results[[animal]] <- time_fit

  # Save only the posterior draws consumed by the inspection script
  sims_time_df <- as_tibble(as.data.frame(time_fit$sims.list))

  write_csv(
    sims_time_df,
    time_output_path
  )
}

# 2) TEHS model ----
jags_selection_results <- list()

for (animal in animals) {
  cat("\\n--------------------------------------------\\n")
  cat("Running TEHS model for animal: ", animal, "\\n", sep = "")

  tssf_path <- file.path(
    prep_complete_dir,
    str_glue("tSSF_data_{animal}_{raster_tag}.csv")
  )
  sims_time_path <- file.path(
    time_model_results_dir,
    str_glue("sims.list_{animal}_{raster_tag}.csv")
  )
  selection_output_path <- file.path(
    selection_model_results_dir,
    str_glue("sims.list_{animal}_{raster_tag}.csv")
  )

  if (file.exists(selection_output_path)) {
    message(
      "Selection model output already exists for ",
      animal,
      ". Skipping refit."
    )
    next
  }

  if (!file.exists(tssf_path) || !file.exists(sims_time_path)) {
    warning(
      "Skipping animal ",
      animal,
      " because required input file is missing."
    )
    next
  }

  tssf_data <- read_csv(tssf_path, show_col_types = FALSE)
  sims_time <- read_csv(sims_time_path, show_col_types = FALSE)

  b0_col <- names(sims_time) |>
    keep(~ .x == "b0")

  b_effect_cols <- names(sims_time) |>
    keep(
      ~ str_detect(.x, "^b\\[[0-9]+\\]$") |
        str_detect(.x, "^b[0-9]+$") |
        str_detect(.x, "^b[.][0-9]+$")
    ) |>
    discard(~ .x == "b0")

  if (length(b0_col) == 0) {
    warning("Skipping animal ", animal, ": missing b0 in time model output.")
    next
  }

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

  g_cols <- names(sims_time) |>
    keep(~ str_detect(.x, "^g[0-9]+$")) |>
    (\(x) x[order(as.numeric(str_remove(x, "^g")))])()

  if (length(b_effect_cols) != length(model_covariates)) {
    warning(
      "Skipping animal ",
      animal,
      ": number of b effects does not match model covariates."
    )
    next
  }

  if (length(g_cols) < 2) {
    warning(
      "Skipping animal ",
      animal,
      ": missing g parameters in time model output."
    )
    next
  }

  required_tssf_cols <- c(
    "time1",
    "dist1",
    "Miss",
    str_glue("{model_covariates}_{0}"),
    str_glue("{model_covariates}_{1}"),
    str_glue("{model_covariates}_{2}"),
    str_glue("{model_covariates}_{3}"),
    str_glue("{model_covariates}_{4}")
  )

  missing_tssf_cols <- setdiff(required_tssf_cols, names(tssf_data))
  if (length(missing_tssf_cols) > 0) {
    warning(
      "Skipping animal ",
      animal,
      " due to missing tSSF columns: ",
      paste(missing_tssf_cols, collapse = ", ")
    )
    next
  }

  betas_mat <- cbind(
    b0 = sims_time[[b0_col]],
    sims_time |>
      select(all_of(b_effect_cols)) |>
      as.matrix()
  )

  gs_mat <- sims_time |>
    select(all_of(g_cols[1:2])) |>
    as.matrix()

  prob <- matrix(NA_real_, nrow(tssf_data), 5)

  for (i in seq_len(nrow(tssf_data))) {
    for (path_idx in 0:4) {
      covar_cols_idx <- str_glue("{model_covariates}_{path_idx}")
      xmat <- c(1, as.numeric(tssf_data[i, covar_cols_idx]))

      mean1 <- tssf_data$dist1[i] * exp(xmat %*% t(betas_mat))

      b1 <- exp(gs_mat[, 1] + gs_mat[, 2] * tssf_data$Miss[i])
      a1 <- b1 * mean1

      prob[i, path_idx + 1] <- median(dgamma(tssf_data$time1[i], a1, b1))
    }
  }

  colnames(prob) <- paste0("prob", 0:4)
  dat <- bind_cols(tssf_data, as_tibble(prob))

  covariate_cols_all_paths <- c(
    str_glue("{model_covariates}_0"),
    str_glue("{model_covariates}_1"),
    str_glue("{model_covariates}_2"),
    str_glue("{model_covariates}_3"),
    str_glue("{model_covariates}_4")
  )

  prob_cols <- paste0("prob", 0:4)

  dat <- dat |>
    filter(if_all(all_of(c(covariate_cols_all_paths, prob_cols)), ~ !is.na(.x)))

  if (nrow(dat) == 0) {
    warning(
      "Skipping animal ",
      animal,
      ": no complete rows for TEHS model after NA filter."
    )
    next
  }

  dat_jags <- list(
    nobs = nrow(dat),
    ncov = length(model_covariates),
    xmat0 = data.matrix(dat[, str_glue("{model_covariates}_0")]),
    xmat1 = data.matrix(dat[, str_glue("{model_covariates}_1")]),
    xmat2 = data.matrix(dat[, str_glue("{model_covariates}_2")]),
    xmat3 = data.matrix(dat[, str_glue("{model_covariates}_3")]),
    xmat4 = data.matrix(dat[, str_glue("{model_covariates}_4")]),
    pmov0 = dat$prob0,
    pmov1 = dat$prob1,
    pmov2 = dat$prob2,
    pmov3 = dat$prob3,
    pmov4 = dat$prob4,
    y = rep(1, nrow(dat))
  )

  set.seed(1)

  params <- c("betas")
  n.iter <- 5000
  n.thin <- 10
  n.burnin <- 1000
  n.chains <- 3

  selection_fit <- jags(
    model.file = selection_model_file,
    parameters.to.save = params,
    data = dat_jags,
    n.chains = n.chains,
    n.burnin = n.burnin,
    n.iter = n.iter,
    n.thin = n.thin,
    DIC = TRUE
  )

  jags_selection_results[[animal]] <- selection_fit

  sims_selection_df <- as_tibble(as.data.frame(selection_fit$sims.list))

  write_csv(
    sims_selection_df,
    selection_output_path
  )
}

cat(
  "\\nTEHS model run finished.\\n",
  "Raster tag: ",
  raster_tag,
  "\\n",
  "Animals with time fit: ",
  length(jags_time_results),
  "\\n",
  "Animals with selection fit: ",
  length(jags_selection_results),
  "\\n",
  sep = ""
)
