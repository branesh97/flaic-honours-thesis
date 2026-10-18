# Forecasting with linear augmented information combination (FLAIC)
# Final run: base forecasts, conventional reconciliation, IComb, PCA construction and FLAIC

library(future)
library(future.apply)
library(dplyr)
library(tidyr)

# load data/settings and functions
data_setup <- readRDS("rds/flaic_01_data_setup.rds")
list2env(data_setup, envir = .GlobalEnv)
rm(data_setup)
source("02_model_functions.R")

future::availableCores()
future::plan(multisession, workers = workers)
run_start_all <- Sys.time()

### Setting up storage ---------------------------------------------------------------------------------------------------------------------------

# 1. comparison forecasts: base + conventional reconciliation + regularised IComb
methods_icomb_reg <- c("icombs_xy1", "icombl_xy1") # icomb ridge and lasso where X and Y are standardised
methods <- c("base", "bu", "ols", "wls_s", "wls_v", "mint_shr", methods_icomb_reg)
fcst_array <- array(NA_real_, dim = c(origin = length(origins), horizon = h, series = ny, method = length(methods)),
                    dimnames = list(origin = as.character(origins), horizon = paste0("h", 1:h), series = colnames(ally), method = methods))
# 2. fitted values from ETS models used as predictors in IComb and FLAIC
fitted_array <- array(NA_real_, dim = c(origin = length(origins), window_size = window_size, series = ny),
                      dimnames = list(origin = as.character(origins), window_size = as.character(1:window_size), series = colnames(ally)))
# 3. IComb lambda information
lambda_info_array <- array(NA_real_, dim = c(origin = length(origins), method = length(methods_icomb_reg), statistic = 6),
                           dimnames = list(origin = as.character(origins), method = methods_icomb_reg,
                                           statistic = c("lambda_max", "lambda_min", "lambda_best", "lambda_best_idx", "path_position", "log10_span")))
# 4. IComb coefficient matrices at selected lambda
coef_array_icomb <- array(vector("list", length(origins) * length(methods_icomb_reg)), dim = c(origin = length(origins), method = length(methods_icomb_reg)),
                          dimnames = list(origin = as.character(origins), method = methods_icomb_reg))
# 5. storage related to principal components
pc_actual_array <- array(vector("list", length(origins) * length(centers) * length(scales)), dim = c(origin = length(origins), center = length(centers), scale = length(scales)),
                         dimnames = list(origin = as.character(origins), center = centers, scale = scales))
pc_fitted_array <- array(vector("list", length(origins) * length(centers) * length(scales)), dim = dim(pc_actual_array), dimnames = dimnames(pc_actual_array))
pc_fcst_array <- array(vector("list", length(origins) * length(centers) * length(scales)), dim = dim(pc_actual_array), dimnames = dimnames(pc_actual_array))
pc_eigval_array <- array(vector("list", length(origins) * length(centers) * length(scales)), dim = dim(pc_actual_array), dimnames = dimnames(pc_actual_array))
pc_eigvec_array <- array(vector("list", length(origins) * length(centers) * length(scales)), dim = dim(pc_actual_array), dimnames = dimnames(pc_actual_array))
pc_propvar_array <- array(vector("list", length(origins) * length(centers) * length(scales)), dim = dim(pc_actual_array), dimnames = dimnames(pc_actual_array))
pc_cumvar_array <- array(vector("list", length(origins) * length(centers) * length(scales)), dim = dim(pc_actual_array), dimnames = dimnames(pc_actual_array))
pc_cormat_lc_array <- array(vector("list", length(origins) * length(centers) * length(scales)), dim = dim(pc_actual_array), dimnames = dimnames(pc_actual_array))
pc_cormat_pc_array <- array(vector("list", length(origins) * length(centers) * length(scales)), dim = dim(pc_actual_array), dimnames = dimnames(pc_actual_array))
pc_scores_array <- array(vector("list", length(origins) * length(centers) * length(scales)), dim = dim(pc_actual_array), dimnames = dimnames(pc_actual_array))
# 6. FLAIC forecasts and lambda information
fcst_array_flaic <- array(NA_real_, dim = c(origin = length(origins), horizon = h, series = ny, penalty = length(penalties), npc = length(p_values), center = length(centers), scale = length(scales)),
                          dimnames = list(origin = as.character(origins), horizon = paste0("h", 1:h), series = colnames(ally), penalty = penalties,
                                          npc = as.character(p_values), center = as.character(centers), scale = as.character(scales)))
lambda_info_array_flaic <- array(NA_real_, dim = c(origin = length(origins), penalty = length(penalties), npc = length(p_values), center = length(centers), scale = length(scales), statistic = 6),
                                 dimnames = list(origin = as.character(origins), penalty = penalties, npc = as.character(p_values), center = centers, scale = scales,
                                                 statistic = c("lambda_max", "lambda_min", "lambda_best", "lambda_best_idx", "path_position", "log10_span")))

### Forecast reconciliation 1: base forecasts + conventional reconciliation -----------------------------------------------------------------------

checkpoint_file <- file.path("checkpoints", "flaic_03_part1_checkpoint.rds")
origin_batches <- split(origins, ceiling(seq_along(origins) / checkpoint_every))
completed_origins <- integer(0)
run_start <- Sys.time()

for (b in seq_along(origin_batches)) {
  batch_origins <- origin_batches[[b]]
  message("Starting batch ", b, "/", length(origin_batches), " | origins ", min(batch_origins), "-", max(batch_origins))
  batch_start <- Sys.time()
  batch_results <- future_lapply(batch_origins, bfcst_ets, future.seed = TRUE)
  origin_times <- vapply(batch_results, function(x) x$elapsed_minutes, numeric(1))
  for (res in batch_results) {
    T_chr <- as.character(res$T)
    fcst_array[T_chr, , , "base"] <- res$bfcst
    fcst_array[T_chr, , , "bu"] <- res$bu_fcst
    fcst_array[T_chr, , , "ols"] <- res$ols_fcst
    fcst_array[T_chr, , , "wls_s"] <- res$wls_s_fcst
    fcst_array[T_chr, , , "wls_v"] <- res$wls_v_fcst
    fcst_array[T_chr, , , "mint_shr"] <- res$mint_shr_fcst
    fitted_array[T_chr, , ] <- res$fitted_mat
  }
  completed_origins <- c(completed_origins, batch_origins)
  saveRDS(list(fcst_array = fcst_array, fitted_array = fitted_array, completed_origins = completed_origins), checkpoint_file)
  message("Finished batch ", b, "/", length(origin_batches), " in ", round(as.numeric(difftime(Sys.time(), batch_start, units = "mins")), 2),
          " minutes | origin times: min=", round(min(origin_times), 2), ", mean=", round(mean(origin_times), 2), ", max=", round(max(origin_times), 2), " mins | checkpoint saved")
}
saveRDS(list(fcst_array = fcst_array, fitted_array = fitted_array, completed_origins = completed_origins), "rds/flaic_03_part1_complete.rds")
message("Part 1 complete: ", length(completed_origins), " origins processed.")

### Forecast reconciliation 2: regularised IComb -------------------------------------------------------------------------------------------------

checkpoint_file <- file.path("checkpoints", "flaic_03_part2_checkpoint.rds")
origin_batches <- split(origins, ceiling(seq_along(origins) / checkpoint_every))
completed_origins <- integer(0)
run_start <- Sys.time()

for (b in seq_along(origin_batches)) {
  batch_origins <- origin_batches[[b]]
  message("Starting batch ", b, "/", length(origin_batches), " | origins ", min(batch_origins), "-", max(batch_origins))
  batch_start <- Sys.time()
  batch_results <- future_lapply(batch_origins, icomb_rec, loss_scale = lambda_loss_scale, future.seed = TRUE)
  origin_times <- vapply(batch_results, function(x) x$elapsed_minutes, numeric(1))
  for (res in batch_results) {
    T_chr <- as.character(res$T)
    fcst_array[T_chr, , , "icombs_xy1"] <- res$icombs_xy1_res$preds
    fcst_array[T_chr, , , "icombl_xy1"] <- res$icombl_xy1_res$preds
    coef_array_icomb[[T_chr, "icombs_xy1"]] <- res$icombs_xy1_res$coefs
    coef_array_icomb[[T_chr, "icombl_xy1"]] <- res$icombl_xy1_res$coefs
    info <- res$icombs_xy1_res$info$lambda_info
    lambda_info_array[T_chr, "icombs_xy1", ] <- c(info$lambda_max, info$lambda_min, info$lambda_best, info$lambda_best_idx, info$path_position, info$log10_span)
    info <- res$icombl_xy1_res$info$lambda_info
    lambda_info_array[T_chr, "icombl_xy1", ] <- c(info$lambda_max, info$lambda_min, info$lambda_best, info$lambda_best_idx, info$path_position, info$log10_span)
  }
  completed_origins <- c(completed_origins, batch_origins)
  saveRDS(list(fcst_array = fcst_array, lambda_info_array = lambda_info_array, coef_array_icomb = coef_array_icomb, lambda_loss_scale = lambda_loss_scale,
               completed_origins = completed_origins), checkpoint_file)
  message("Finished batch ", b, "/", length(origin_batches), " in ", round(as.numeric(difftime(Sys.time(), batch_start, units = "mins")), 2),
          " minutes | origin times: min=", round(min(origin_times), 2), ", mean=", round(mean(origin_times), 2), ", max=", round(max(origin_times), 2), " mins | checkpoint saved")
}
saveRDS(list(fcst_array = fcst_array, lambda_info_array = lambda_info_array, coef_array_icomb = coef_array_icomb, lambda_loss_scale = lambda_loss_scale,
             completed_origins = completed_origins), "rds/flaic_03_part2_complete.rds")
message("Part 2 complete: ", length(completed_origins), " origins processed.")

### Forecast reconciliation 3: generating principal-component series -----------------------------------------------------------------------------

checkpoint_file <- file.path("checkpoints", "flaic_03_part3_checkpoint.rds")
task_grid <- expand.grid(T = origins, center_op = centers, scale_op = scales)
task_ids <- seq_len(nrow(task_grid))
grid_batches <- split(task_ids, ceiling(seq_along(task_ids) / checkpoint_every))
completed_grid <- integer(0)
run_start <- Sys.time()

for (b in seq_along(grid_batches)) {
  batch_ids <- grid_batches[[b]]
  batch_grid <- task_grid[batch_ids, , drop = FALSE]
  message("Starting batch ", b, "/", length(grid_batches), " | tasks ", min(batch_ids), "-", max(batch_ids), " | T ", min(batch_grid$T), "-", max(batch_grid$T),
          " | center = ", paste(unique(batch_grid$center_op), collapse = ","), " | scale = ", paste(unique(batch_grid$scale_op), collapse = ","))
  batch_start <- Sys.time()
  batch_results <- future_lapply(seq_len(nrow(batch_grid)), function(k) {
    gen_pc(T = batch_grid$T[k], center = batch_grid$center_op[k], scale = batch_grid$scale_op[k], use_pc_scores = use_pc_scores)
  }, future.seed = TRUE)
  ind_run_times <- vapply(batch_results, function(x) x$elapsed_minutes, numeric(1))
  for (res in batch_results) {
    T_chr <- as.character(res$T)
    center_chr <- as.character(res$center)
    scale_chr <- as.character(res$scale)
    pc_actual_array[[T_chr, center_chr, scale_chr]] <- res$pca_comp
    pc_scores_array[[T_chr, center_chr, scale_chr]] <- res$pc_scores
    pc_fitted_array[[T_chr, center_chr, scale_chr]] <- res$pc_fitted_matrix
    pc_fcst_array[[T_chr, center_chr, scale_chr]] <- res$pc_fcst_matrix
    pc_eigval_array[[T_chr, center_chr, scale_chr]] <- res$eigval
    pc_eigvec_array[[T_chr, center_chr, scale_chr]] <- res$eigvec
    pc_propvar_array[[T_chr, center_chr, scale_chr]] <- res$propvar
    pc_cumvar_array[[T_chr, center_chr, scale_chr]] <- res$cumvar
    pc_cormat_lc_array[[T_chr, center_chr, scale_chr]] <- res$cor_mat_lc
    pc_cormat_pc_array[[T_chr, center_chr, scale_chr]] <- res$cor_mat_pc
  }
  completed_grid <- c(completed_grid, batch_ids)
  saveRDS(list(use_pc_scores = use_pc_scores, pc_actual_array = pc_actual_array, pc_scores_array = pc_scores_array, pc_fitted_array = pc_fitted_array,
               pc_fcst_array = pc_fcst_array, pc_eigval_array = pc_eigval_array, pc_eigvec_array = pc_eigvec_array, pc_propvar_array = pc_propvar_array,
               pc_cumvar_array = pc_cumvar_array, pc_cormat_lc_array = pc_cormat_lc_array, pc_cormat_pc_array = pc_cormat_pc_array, completed_grid = completed_grid), checkpoint_file)
  message("Finished batch ", b, "/", length(grid_batches), " in ", round(as.numeric(difftime(Sys.time(), batch_start, units = "mins")), 2),
          " minutes | individual runs: min=", round(min(ind_run_times), 2), ", mean=", round(mean(ind_run_times), 2), ", max=", round(max(ind_run_times), 2), " mins | checkpoint saved")
}
saveRDS(list(use_pc_scores = use_pc_scores, pc_actual_array = pc_actual_array, pc_scores_array = pc_scores_array, pc_fitted_array = pc_fitted_array,
             pc_fcst_array = pc_fcst_array, pc_eigval_array = pc_eigval_array, pc_eigvec_array = pc_eigvec_array, pc_propvar_array = pc_propvar_array,
             pc_cumvar_array = pc_cumvar_array, pc_cormat_lc_array = pc_cormat_lc_array, pc_cormat_pc_array = pc_cormat_pc_array, completed_grid = completed_grid),
        "rds/flaic_03_part3_complete.rds")
message("Part 3 complete: ", length(completed_grid), " PC runs processed.")

### Forecast reconciliation 4: FLAIC -------------------------------------------------------------------------------------------------------------

# ensure all final p values are supported by the rolling PCA objects
stopifnot(max(p_values) <= nb, all(vapply(pc_fitted_array, function(x) !is.null(x) && ncol(x) >= max(p_values), logical(1))),
          all(vapply(pc_fcst_array, function(x) !is.null(x) && ncol(x) >= max(p_values), logical(1))))
message("PC storage supports all final FLAIC p values.")

flaic_specs <- tidyr::crossing(center = centers, scale = scales, p = p_values)
completed_origins <- integer(0)
run_start <- Sys.time()

for (T in origins) {
  T_chr <- as.character(T)
  origin_start <- Sys.time()
  message("Starting FLAIC origin ", T, " | ", match(T, origins), "/", length(origins), " | ", nrow(flaic_specs), " specifications")
  results_T <- future_lapply(seq_len(nrow(flaic_specs)), function(k) {
    flaic_rec(T = T, p = flaic_specs$p[k], center = flaic_specs$center[k], scale = flaic_specs$scale[k], loss_scale = lambda_loss_scale)
  }, future.seed = TRUE, future.scheduling = Inf)
  coefficients_T <- vector("list", length(results_T))
  for (i in seq_along(results_T)) {
    res <- results_T[[i]]
    center_chr <- as.character(res$center)
    scale_chr <- as.character(res$scale)
    p_chr <- as.character(res$p)
    fcst_array_flaic[T_chr, , , "ridge", p_chr, center_chr, scale_chr] <- res$flaics_xy1_res$preds
    fcst_array_flaic[T_chr, , , "lasso", p_chr, center_chr, scale_chr] <- res$flaicl_xy1_res$preds
    info <- res$flaics_xy1_res$info$lambda_info
    lambda_info_array_flaic[T_chr, "ridge", p_chr, center_chr, scale_chr, ] <- c(info$lambda_max, info$lambda_min, info$lambda_best, info$lambda_best_idx, info$path_position, info$log10_span)
    info <- res$flaicl_xy1_res$info$lambda_info
    lambda_info_array_flaic[T_chr, "lasso", p_chr, center_chr, scale_chr, ] <- c(info$lambda_max, info$lambda_min, info$lambda_best, info$lambda_best_idx, info$path_position, info$log10_span)
    coefficients_T[[i]] <- list(p = res$p, center = res$center, scale = res$scale, ridge = res$flaics_xy1_res$coefs, lasso = res$flaicl_xy1_res$coefs)
  }
  saveRDS(coefficients_T, file.path("coefficients", paste0("flaic_coefficients_origin_", T_chr, ".rds")), compress = TRUE)
  completed_origins <- c(completed_origins, T)
  saveRDS(list(fcst_array_flaic = fcst_array_flaic, lambda_info_array_flaic = lambda_info_array_flaic, p_values = p_values, completed_origins = completed_origins),
          file.path("checkpoints", "flaic_03_part4_checkpoint.rds"))
  message("Finished FLAIC origin ", T, " in ", round(as.numeric(difftime(Sys.time(), origin_start, units = "mins")), 2), " minutes | coefficients saved")
}
stopifnot(length(completed_origins) == length(origins), sum(is.na(fcst_array_flaic)) == 0, sum(is.na(lambda_info_array_flaic)) == 0)
saveRDS(list(fcst_array_flaic = fcst_array_flaic, lambda_info_array_flaic = lambda_info_array_flaic, p_values = p_values, completed_origins = completed_origins),
        "rds/flaic_03_part4_complete.rds")
message("Part 4 complete: FLAIC estimated for all origins and final specifications.")

### Final save ------------------------------------------------------------------------------------------------------------------------------------

run_end_all <- Sys.time()
run_info <- list(run_start = run_start_all, run_end = run_end_all, elapsed_hours = as.numeric(difftime(run_end_all, run_start_all, units = "hours")), session_info = sessionInfo())
saveRDS(list(methods = methods, methods_icomb_reg = methods_icomb_reg, fcst_array = fcst_array, fitted_array = fitted_array, lambda_info_array = lambda_info_array,
             coef_array_icomb = coef_array_icomb, pc_actual_array = pc_actual_array, pc_scores_array = pc_scores_array, pc_fitted_array = pc_fitted_array,
             pc_fcst_array = pc_fcst_array, pc_eigval_array = pc_eigval_array, pc_eigvec_array = pc_eigvec_array, pc_propvar_array = pc_propvar_array,
             pc_cumvar_array = pc_cumvar_array, pc_cormat_lc_array = pc_cormat_lc_array, pc_cormat_pc_array = pc_cormat_pc_array,
             fcst_array_flaic = fcst_array_flaic, lambda_info_array_flaic = lambda_info_array_flaic, run_info = run_info), "rds/flaic_03_forecasts_complete.rds")
message("03_run_forecasts complete. Total elapsed hours: ", round(run_info$elapsed_hours, 2))

# save file for global environment after the complete 03
# Full-session .RData omitted: saved RDS outputs are sufficient.