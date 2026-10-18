# Supplementary MCB analysis: FLAIC and IComb versus ETS base forecasts
# Run from repository root AFTER 01_data_setup.R and 05_evaluate.R.
# Mirrors the Friedman/Nemenyi rank calculations in 07b_analysis_workbench.R.
# Including the Base method increases the number of methods (and changes the CD).

data_setup <- readRDS("rds/flaic_01_data_setup.rds")
list2env(data_setup, envir = .GlobalEnv)
rm(data_setup)
evaluation_results <- readRDS("rds/flaic_05_evaluation_complete.rds")
list2env(evaluation_results, envir = .GlobalEnv)
rm(evaluation_results)
dir.create("supplementary", showWarnings = FALSE)

mcb_horizons <- c(1L, 6L, 12L)
mcb_p <- as.character(p_values)
confidence <- 0.95

series_msfe <- function(errors, hor) {
  n_valid <- length(origins) - (hor - 1L)
  errors <- as.matrix(errors)
  stopifnot(nrow(errors) == length(origins), ncol(errors) == ny)
  x <- errors[seq_len(n_valid), , drop = FALSE]
  stopifnot(all(is.finite(x)))
  colMeans(x^2)
}

compare_to_base <- function(hor, estimator, center, scale) {
  hname <- paste0("h", hor)
  stopifnot(estimator %in% c("ridge", "lasso"))
  icomb_name <- if (estimator == "ridge") "icombs_xy1" else "icombl_xy1"

  base_msfe <- series_msfe(error_array[, hname, , "base"], hor)
  icomb_msfe <- series_msfe(error_array[, hname, , icomb_name], hor)
  flaic_msfe <- vapply(mcb_p, function(p) {
    series_msfe(error_array_flaic[, hname, , estimator, p,
                                  as.character(center), as.character(scale)], hor)
  }, numeric(ny))
  msfe <- cbind(Base = base_msfe, IComb = icomb_msfe, flaic_msfe)
  colnames(msfe) <- c("Base", "IComb", paste0("p=", mcb_p))
  stopifnot(nrow(msfe) == ny, all(is.finite(msfe)))

  # Rows are series (ranking blocks) and columns are methods.
  ranks <- t(apply(msfe, 1L, rank, ties.method = "average"))
  average_rank <- colMeans(ranks)
  n <- nrow(msfe)
  k <- ncol(msfe)
  critical_distance <- stats::qtukey(confidence, nmeans = k, df = Inf) / sqrt(2) *
    sqrt(k * (k + 1) / (6 * n))

  base_rank <- unname(average_rank["Base"])
  delta_rank <- base_rank - average_rank
  status <- ifelse(names(average_rank) == "Base", "Reference: base forecasts",
            ifelse(delta_rank > critical_distance, "Significantly better than base",
              ifelse(delta_rank > 0, "Better ranked, not significantly better",
                     "Not better ranked than base")))
  preprocessing <- paste0(if (center) "Centred" else "Uncentred", ", ",
                          if (scale) "scaled" else "unscaled")

  data.frame(
    horizon = hor, estimator = estimator, PCA_specification = preprocessing,
    center = center, scale = scale, method = names(average_rank),
    avg_rank = as.numeric(average_rank),
    base_avg_rank = base_rank,
    base_minus_method_rank = as.numeric(delta_rank),
    nemenyi_critical_distance = critical_distance,
    friedman_p_value = stats::friedman.test(msfe)$p.value,
    comparison_vs_base = status,
    stringsAsFactors = FALSE
  )
}

specs <- expand.grid(horizon = mcb_horizons, estimator = c("ridge", "lasso"),
                     center = c(TRUE, FALSE), scale = c(FALSE, TRUE),
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
results <- do.call(rbind, lapply(seq_len(nrow(specs)), function(i) {
  compare_to_base(specs$horizon[i], specs$estimator[i],
                  specs$center[i], specs$scale[i])
}))
rownames(results) <- NULL
stopifnot(nrow(results) == 24L * (length(mcb_p) + 2L))

write.csv(results, "supplementary/MCB_with_base_all_comparisons.csv", row.names = FALSE)

exceptions <- subset(results, grepl("^p=", method) &
                     comparison_vs_base != "Significantly better than base")
write.csv(exceptions, "supplementary/MCB_with_base_exceptions.csv", row.names = FALSE)

report <- c(
  "FLAIC THESIS: MCB COMPARISONS INCLUDING ETS BASE FORECASTS",
  "======================================================",
  "",
  "Friedman and Nemenyi comparisons, 112 hierarchy series as ranking blocks.",
  "Horizons: 1, 6 and 12; ridge/lasso separately; four PCA specifications.",
  "Each panel contains Base, matching IComb and eight FLAIC p specifications.",
  "Significantly better than base: rank(Base) - rank(method) > Nemenyi CD.",
  "Not significantly better does NOT imply significantly worse than base.",
  "",
  paste("Total FLAIC configurations:", sum(grepl("^p=", results$method))),
  paste("FLAIC configurations not significantly better than base:", nrow(exceptions)),
  ""
)
if (nrow(exceptions) > 0L) {
  report <- c(report, "EXCEPTIONS (unrounded values)",
              capture.output(print(exceptions[, c(
                "horizon", "estimator", "PCA_specification", "method",
                "avg_rank", "base_avg_rank", "base_minus_method_rank",
                "nemenyi_critical_distance", "comparison_vs_base")],
                row.names = FALSE, digits = 12)))
} else {
  report <- c(report, "No exceptions: all FLAIC specifications significantly better than base.")
}
report <- c(report, "", "Statistical caution: aggregation makes hierarchy series dependent.",
            "These comparisons are reported as descriptive rank-based evidence.")
writeLines(report, "supplementary/MCB_with_base_report.txt")
cat(paste(report, collapse = "\n"), "\n")

# Verify against reported thesis footnote, without silently overwriting results.
expected_keys <- c("6|ridge|Centred, unscaled|p=77",
                   "6|ridge|Uncentred, unscaled|p=77")
actual_keys <- paste(exceptions$horizon, exceptions$estimator,
                     exceptions$PCA_specification, exceptions$method, sep = "|")
if (!setequal(actual_keys, expected_keys)) {
  warning("MCB exceptions differ from those reported in the submitted thesis; investigate before publication.")
}
message("09_mcb_with_base complete; CSV and TXT outputs saved.")
