# Filter telemetry data using ctmm distance and speed diagnostics
# The function returns only the cleaned move2 objects needed by the next script

ctmm_remove_outliers <- function(
  x,
  individual_names = names(x),
  dist = 3000,
  speed = 1.5
) {
  telemetry_mov2_no_out <- list()

  for (animal in individual_names) {
    cli::cli_h1("{animal} started")

    # Read the telemetry object for one animal and diagnose extreme relocations
    telemetry_data <- x[[animal]]
    initial_outliers <- ctmm::outlie(telemetry_data)

    # Apply the distance filter first to remove spatial spikes
    telemetry_no_out_dist <- telemetry_data[initial_outliers$distance < dist, ]

    # Recalculate outliers after the distance filter and then remove fast steps
    distance_filtered_outliers <- ctmm::outlie(telemetry_no_out_dist)
    telemetry_no_out_dist_speed <- telemetry_no_out_dist[
      distance_filtered_outliers$speed < speed,
    ]

    # Convert the filtered telemetry back to move2 so the main script can keep
    # Using a tidy table-oriented workflow
    telemetry_mov2_no_out[[animal]] <- move2::mt_as_move2(
      telemetry_no_out_dist_speed
    )

    cli::cli_alert_success("Outlier filtering finished for {animal}")
  }

  telemetry_mov2_no_out
}
