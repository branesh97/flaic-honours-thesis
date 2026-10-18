# Forecasting with linear augmented information combination (FLAIC)
# Final FLAP benchmark run

# load data/settings, forecasting results and functions
data_setup <- readRDS("rds/flaic_01_data_setup.rds")
list2env(data_setup, envir = .GlobalEnv)
rm(data_setup)
forecast_results <- readRDS("rds/flaic_03_forecasts_complete.rds")
list2env(forecast_results, envir = .GlobalEnv)
rm(forecast_results)
source("02_model_functions.R")

### FLAP benchmark -------------------------------------------------------------------------------------------------------------------------------

flap_p_values <- p_values
fcst_array_flap <- array(NA_real_, dim = c(origin = length(origins), horizon = h, series = ny, p = length(flap_p_values), center = length(centers), scale = length(scales)),
                         dimnames = list(origin = as.character(origins), horizon = paste0("h", seq_len(h)), series = colnames(ally), p = as.character(flap_p_values),
                                         center = as.character(centers), scale = as.character(scales)))
completed_p <- integer(0)
run_start <- Sys.time()

for (p in flap_p_values) {
  p_chr <- as.character(p)
  p_start <- Sys.time()
  message("Starting FLAP benchmark for p = ", p)
  flap_result <- run_flap_benchmark(p_flap = p, verbose = TRUE)
  fcst_array_flap[ , , , p_chr, , ] <- flap_result$fcst_array_star
  completed_p <- c(completed_p, p)
  saveRDS(list(flap_p_values = flap_p_values, fcst_array_flap = fcst_array_flap, completed_p = completed_p), file.path("checkpoints", "flaic_04_flap_checkpoint.rds"))
  message("Finished FLAP p = ", p, " in ", round(as.numeric(difftime(Sys.time(), p_start, units = "mins")), 2), " minutes | checkpoint saved")
}

stopifnot(length(completed_p) == length(flap_p_values), sum(is.na(fcst_array_flap)) == 0)
saveRDS(list(flap_p_values = flap_p_values, fcst_array_flap = fcst_array_flap, completed_p = completed_p), "rds/flaic_04_flap_complete.rds")
message("04_run_flap complete in ", round(as.numeric(difftime(Sys.time(), run_start, units = "hours")), 2), " hours.")

# save file for global environment after the complete 04
# Full-session .RData omitted: saved RDS outputs are sufficient.