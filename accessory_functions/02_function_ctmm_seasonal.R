# Process each period and save outputs

process_period_data <- function(data) {
  periods <- unique(data$period)

  for (period_name in periods) {
    # Keep one period at a time and rename the coordinate columns for ctmm
    data_filter <- data |>
      filter(period == period_name) |>
      mutate(timestamp = lubridate::as_datetime(timestamp, tz = "UTC")) |>
      group_by(individual_local_identifier) |>
      arrange(timestamp) |>
      ungroup() |>
      select(
        individual.local.identifier = individual_local_identifier,
        timestamp,
        location.long = longitude,
        location.lat = latitude
      )

    if (nrow(data_filter) < 10) {
      message(sprintf(
        "\nSkipping period %s due to insufficient data.\n",
        period_name
      ))
      next
    }

    message(sprintf("\nStarting analysis for period: %s\n", period_name))

    # Convert the focal-animal data to telemetry objects for model fitting
    data_telemetry <- ctmm::as.telemetry(
      data_filter,
      timezone = "America/Sao_Paulo"
    )

    if (inherits(data_telemetry, "telemetry")) {
      animal_name <- unique(data_filter$individual.local.identifier)
      data_telemetry <- list(data_telemetry)
      names(data_telemetry) <- animal_name
    }

    data_identity <- names(data_telemetry)

    # Save telemetry objects because they are used later
    dir.create(
      "data/processed/telemetry_object",
      recursive = TRUE,
      showWarnings = FALSE
    )
    readr::write_rds(
      data_telemetry,
      sprintf(
        "data/processed/telemetry_object/data_telemetry_%s.rds",
        period_name
      )
    )

    # Fit one ctmm model per animal in the period
    message(sprintf("\nFitting models for period: %s\n", period_name))
    fit_list <- list()

    for (animal in data_identity) {
      my_tel <- data_telemetry[[animal]]
      my_guess <- ctmm::ctmm.guess(my_tel, interactive = FALSE)
      fit_list[[animal]] <- ctmm::ctmm.select(my_tel, my_guess, verbose = FALSE)
    }

    # Save fitted models because the final figure script can reuse them when it needs annual objects that are not yet on disk
    dir.create(
      "output/models/fit_models",
      recursive = TRUE,
      showWarnings = FALSE
    )
    readr::write_rds(
      fit_list,
      sprintf("output/models/fit_models/fit_list_%s.rds", period_name)
    )

    # Save AKDE objects because the overlap script and figure script both depend on them
    message(sprintf("\nCalculating AKDE for period: %s\n", period_name))
    akde_data <- ctmm::akde(data_telemetry, CTMM = fit_list, weights = TRUE)

    dir.create("output/models/fit_akde", recursive = TRUE, showWarnings = FALSE)
    readr::write_rds(
      akde_data,
      file = sprintf("output/models/fit_akde/list_akde_%s.rds", period_name)
    )
  }
}
