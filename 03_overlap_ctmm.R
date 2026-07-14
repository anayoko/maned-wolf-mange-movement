# Calculate seasonal home-range overlap between the female maned wolves.
# The exported table is used in Figure 1

# Load libraries
library(ctmm)
library(tidyverse)

# Load output data
akde_files <- list.files(
  "output/models/fit_akde",
  pattern = "list_akde_.*\\.rds$",
  full.names = TRUE
)

final_overlap_results <- akde_files |>
  map_dfr(
    \(file_path) {
      period_name <- stringr::str_extract(
        basename(file_path),
        "(?<=list_akde_).*(?=\\.rds)"
      )
      period_name <- stringr::str_replace(period_name, "Rainy", "Wet")

      if (stringr::str_detect(period_name, "Total$")) {
        return(tibble())
      }

      akde_list <- readRDS(file_path)

      if (inherits(akde_list[[1]], "list")) {
        akde_list <- purrr::list_flatten(akde_list)
      }

      # Skip periods with fewer than two focal animals because overlap cannot
      # be calculated for a single home range
      if (length(akde_list) < 2) {
        return(tibble())
      }

      overlap_res <- ctmm::overlap(akde_list)

      map_dfr(
        c("low", "est", "high"),
        \(ci_level) {
          ci_matrix <- overlap_res$CI[,, ci_level]
          upper_index <- which(
            upper.tri(ci_matrix, diag = FALSE),
            arr.ind = TRUE
          )

          tibble(
            animal_1 = rownames(ci_matrix)[upper_index[, 1]],
            animal_2 = colnames(ci_matrix)[upper_index[, 2]],
            value = ci_matrix[upper_index],
            ci_level = ci_level,
            period = period_name
          )
        }
      )
    }
  ) |>
  pivot_wider(names_from = ci_level, values_from = value) |>
  mutate(
    period = factor(
      period,
      levels = c(
        "2023_Wet",
        "2024_Dry",
        "2024_Wet"
      ),
      ordered = TRUE
    )
  ) |>
  arrange(period, animal_1, animal_2) |>
  mutate(period = as.character(period))

dir.create("output/tables/overlap", recursive = TRUE, showWarnings = FALSE)

final_overlap_results |>
  write_csv("output/tables/overlap/overlap_results_period_individuals.csv")
