# Forecasting with linear augmented information combination (FLAIC)

# Standalone generation of raw PRIAL tables (requires saved stage 01 and 05 outputs).
library(tidyverse)
data_setup <- readRDS(file.path("rds", "flaic_01_data_setup.rds"))
list2env(data_setup, envir = .GlobalEnv)
rm(data_setup)
evaluation_results <- readRDS(file.path("rds", "flaic_05_evaluation_complete.rds"))
list2env(evaluation_results, envir = .GlobalEnv)
rm(evaluation_results)

### Raw-error PRIAL tables ------------------------------------------------------------------------------------------------------------------------

dir.create("rds", showWarnings = FALSE)

prial_table_h <- c(paste0("h", 1:12), "1-6", "1-12")
prial_table_groups <- c("Average", "Australia", "States", "Zones", "Regions")
prial_table_p <- dimnames(prial_array_flaic)$npc
prial_table_cs <- tibble::tibble(preprocessing = c("Centred, unscaled", "Centred, scaled", "Uncentred, unscaled", "Uncentred, scaled"),
                                 center = c(TRUE, TRUE, FALSE, FALSE), scale = c(FALSE, TRUE, FALSE, TRUE))

make_prial_rows <- function() {
  conventional <- tibble::tibble(Family = "Conventional reconciliation", row_type = "comparison",
                                 Method = c("BU", "OLS", "WLS$_s$", "WLS$_v$", "MinT$_{shr}$"),
                                 preprocessing = "", p = NA_character_,
                                 comparison_method = c("bu", "ols", "wls_s", "wls_v", "mint_shr"),
                                 penalty = NA_character_, center = NA, scale = NA)
  
  make_reg <- function(penalty, family, icomb, flaic) {
    icomb_row <- tibble::tibble(Family = family, row_type = "comparison", Method = icomb, preprocessing = "", p = NA_character_,
                                comparison_method = ifelse(penalty == "ridge", "icombs_xy1", "icombl_xy1"),
                                penalty = penalty, center = NA, scale = NA)
    flaic_rows <- purrr::map_dfr(seq_len(nrow(prial_table_cs)), function(i)
      tibble::tibble(Family = family, row_type = "flaic", Method = flaic, preprocessing = prial_table_cs$preprocessing[i],
                     p = prial_table_p, comparison_method = NA_character_, penalty = penalty,
                     center = prial_table_cs$center[i], scale = prial_table_cs$scale[i]))
    dplyr::bind_rows(icomb_row, flaic_rows)
  }
  
  dplyr::bind_rows(conventional,
                   make_reg("ridge", "Ridge", "IComb-S", "FLAIC-S"),
                   make_reg("lasso", "Lasso", "IComb-L", "FLAIC-L"))
}

make_prial_table <- function(group) {
  out <- make_prial_rows()
  for (hor in prial_table_h) out[[hor]] <- vapply(seq_len(nrow(out)), function(i) {
    if (out$row_type[i] == "comparison") return(as.numeric(prial_array[hor, group, out$comparison_method[i]]))
    as.numeric(prial_array_flaic[hor, group, out$penalty[i], out$p[i], as.character(out$center[i]), as.character(out$scale[i])])
  }, numeric(1))
  out
}

prial_raw_tables <- stats::setNames(lapply(prial_table_groups, make_prial_table), prial_table_groups)

stopifnot(all(vapply(prial_raw_tables, function(x) all(is.finite(as.matrix(x[, prial_table_h]))), logical(1))))
saveRDS(prial_raw_tables, "rds/table_prial_raw_errors.rds")

message("Raw-error PRIAL tables saved for Quarto.")