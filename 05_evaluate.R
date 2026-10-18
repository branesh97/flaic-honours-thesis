# Forecasting with linear augmented information combination (FLAIC)
# Final forecast evaluation: raw MSFE/PRIAL, FLAP benchmark and scaled-error robustness

# load data/settings and completed forecasts
data_setup <- readRDS("rds/flaic_01_data_setup.rds")
list2env(data_setup, envir = .GlobalEnv)
rm(data_setup)
forecast_results <- readRDS("rds/flaic_03_forecasts_complete.rds")
list2env(forecast_results, envir = .GlobalEnv)
rm(forecast_results)
flap_results <- readRDS("rds/flaic_04_flap_complete.rds")
list2env(flap_results, envir = .GlobalEnv)
rm(flap_results)

### Forecast evaluation 1: actuals and raw forecast errors ----------------------------------------------------------------------------------------

# actual_array stores actuals to compute reconciliation error; actual/target = y_T+h where T is forecast origin and h forecast horizon
actual_array <- array(NA_real_, dim = c(origin = length(origins), horizon = h, series = ny),
                      dimnames = list(origin = as.character(origins), horizon = paste0("h", 1:h), series = colnames(ally)))
for (T in origins) {
  for (hor in 1:h) {
    target <- T + hor
    if (target <= nT) actual_array[as.character(T), paste0("h", hor), ] <- ally[target, ]
  }
}

# forecast errors
error_array <- array(NA_real_, dim = dim(fcst_array), dimnames = dimnames(fcst_array))
error_array_flaic <- array(NA_real_, dim = dim(fcst_array_flaic), dimnames = dimnames(fcst_array_flaic))
error_array_flap <- array(NA_real_, dim = dim(fcst_array_flap), dimnames = dimnames(fcst_array_flap))
for (method in methods) error_array[, , , method] <- actual_array - fcst_array[, , , method]
for (penalty in penalties) {
  for (npc in as.character(p_values)) {
    for (center in centers) {
      for (scale in scales) error_array_flaic[, , , penalty, npc, as.character(center), as.character(scale)] <- actual_array - fcst_array_flaic[, , , penalty, npc, as.character(center), as.character(scale)]
    }
  }
}
for (p in as.character(flap_p_values)) {
  for (center in centers) {
    for (scale in scales) error_array_flap[, , , p, as.character(center), as.character(scale)] <- actual_array - fcst_array_flap[, , , p, as.character(center), as.character(scale)]
  }
}

### Forecast evaluation 2: raw MSFE and PRIAL ----------------------------------------------------------------------------------------------------

horizon_names <- c(paste0("h", 1:h), "1-6", "1-12")
# MSFE storage
msfe_array <- array(NA_real_, dim = c(horizon = h+2, groups = length(group_idx), method = length(methods)),
                    dimnames = list(horizon = horizon_names, groups = names(group_idx), method = methods))
msfe_array_flaic <- array(NA_real_, dim = c(horizon = h+2, groups = length(group_idx), penalty = length(penalties), npc = length(p_values), center = length(centers), scale = length(scales)),
                          dimnames = list(horizon = horizon_names, groups = names(group_idx), penalty = penalties, npc = as.character(p_values), center = centers, scale = scales))
msfe_array_flap <- array(NA_real_, dim = c(horizon = h+2, groups = length(group_idx), p = length(flap_p_values), center = length(centers), scale = length(scales)),
                         dimnames = list(horizon = horizon_names, groups = names(group_idx), p = as.character(flap_p_values), center = centers, scale = scales))

# comparison methods
for (method in methods) {
  for (group in names(group_idx)) {
    for (hor in paste0("h", 1:h)) {
      selected_errors <- error_array[, hor, group_idx[[group]], method]
      msfe_array[hor, group, method] <- mean(selected_errors^2, na.rm = TRUE)
    }
    selected_errors <- error_array[, paste0("h", 1:6), group_idx[[group]], method]
    msfe_array["1-6", group, method] <- mean(selected_errors^2, na.rm = TRUE)
    selected_errors <- error_array[, paste0("h", 1:h), group_idx[[group]], method]
    msfe_array["1-12", group, method] <- mean(selected_errors^2, na.rm = TRUE)
  }
}
# FLAIC
for (penalty in penalties) {
  for (npc in as.character(p_values)) {
    for (center in centers) {
      for (scale in scales) {
        center_chr <- as.character(center)
        scale_chr <- as.character(scale)
        for (hor in paste0("h", 1:h)) {
          for (group in names(group_idx)) {
            selected_errors <- error_array_flaic[, hor, group_idx[[group]], penalty, npc, center_chr, scale_chr]
            msfe_array_flaic[hor, group, penalty, npc, center_chr, scale_chr] <- mean(selected_errors^2, na.rm = TRUE)
          }
        }
        for (group in names(group_idx)) {
          selected_errors <- error_array_flaic[, paste0("h", 1:6), group_idx[[group]], penalty, npc, center_chr, scale_chr]
          msfe_array_flaic["1-6", group, penalty, npc, center_chr, scale_chr] <- mean(selected_errors^2, na.rm = TRUE)
          selected_errors <- error_array_flaic[, paste0("h", 1:h), group_idx[[group]], penalty, npc, center_chr, scale_chr]
          msfe_array_flaic["1-12", group, penalty, npc, center_chr, scale_chr] <- mean(selected_errors^2, na.rm = TRUE)
        }
      }
    }
  }
}
# FLAP
for (p in as.character(flap_p_values)) {
  for (center in centers) {
    for (scale in scales) {
      center_chr <- as.character(center)
      scale_chr <- as.character(scale)
      for (hor in paste0("h", 1:h)) {
        for (group in names(group_idx)) {
          selected_errors <- error_array_flap[, hor, group_idx[[group]], p, center_chr, scale_chr]
          msfe_array_flap[hor, group, p, center_chr, scale_chr] <- mean(selected_errors^2, na.rm = TRUE)
        }
      }
      for (group in names(group_idx)) {
        selected_errors <- error_array_flap[, paste0("h", 1:6), group_idx[[group]], p, center_chr, scale_chr]
        msfe_array_flap["1-6", group, p, center_chr, scale_chr] <- mean(selected_errors^2, na.rm = TRUE)
        selected_errors <- error_array_flap[, paste0("h", 1:h), group_idx[[group]], p, center_chr, scale_chr]
        msfe_array_flap["1-12", group, p, center_chr, scale_chr] <- mean(selected_errors^2, na.rm = TRUE)
      }
    }
  }
}

# PRIAL storage
prial_array <- array(NA_real_, dim = dim(msfe_array), dimnames = dimnames(msfe_array))
prial_array_flaic <- array(NA_real_, dim = dim(msfe_array_flaic), dimnames = dimnames(msfe_array_flaic))
prial_array_flap <- array(NA_real_, dim = dim(msfe_array_flap), dimnames = dimnames(msfe_array_flap))
# comparison methods; base PRIAL is zero by definition
for (method in methods) {
  for (hor in horizon_names) {
    for (group in names(group_idx)) {
      if (method == "base") prial_array[hor, group, method] <- 0 else prial_array[hor, group, method] <- 100 * (msfe_array[hor, group, method]/msfe_array[hor, group, "base"] - 1)
    }
  }
}
# FLAIC
for (penalty in penalties) {
  for (npc in as.character(p_values)) {
    for (center in centers) {
      for (scale in scales) {
        center_chr <- as.character(center)
        scale_chr <- as.character(scale)
        for (hor in horizon_names) {
          for (group in names(group_idx)) prial_array_flaic[hor, group, penalty, npc, center_chr, scale_chr] <- 100 * (msfe_array_flaic[hor, group, penalty, npc, center_chr, scale_chr]/msfe_array[hor, group, "base"] - 1)
        }
      }
    }
  }
}
# FLAP
for (p in as.character(flap_p_values)) {
  for (center in centers) {
    for (scale in scales) {
      center_chr <- as.character(center)
      scale_chr <- as.character(scale)
      for (hor in horizon_names) {
        for (group in names(group_idx)) prial_array_flap[hor, group, p, center_chr, scale_chr] <- 100 * (msfe_array_flap[hor, group, p, center_chr, scale_chr]/msfe_array[hor, group, "base"] - 1)
      }
    }
  }
}

saveRDS(list(actual_array = actual_array, error_array = error_array, error_array_flaic = error_array_flaic, error_array_flap = error_array_flap,
             msfe_array = msfe_array, msfe_array_flaic = msfe_array_flaic, msfe_array_flap = msfe_array_flap,
             prial_array = prial_array, prial_array_flaic = prial_array_flaic, prial_array_flap = prial_array_flap), "rds/flaic_05_raw_evaluation.rds")
message("Raw MSFE and PRIAL evaluation complete.")

### Forecast evaluation 3: scaled-error robustness -----------------------------------------------------------------------------------------------

#### Option 1: seasonal-naive scaling -------------------------------------------------------------------------------------------------------------
# q^2 = forecast_error^2 / d, where d is the in-sample MSE from a seasonal-naive forecast
# this asks how large the forecast error was relative to how well a simple seasonal-naive benchmark could forecast the series

seasonal_lag <- seasonal_period
# d_{T,j}: seasonal-naive scaling denominator; dimensions are origin x series and it is independent of reconciliation method and forecast horizon
d_array <- array(NA_real_, dim = c(origin = length(origins), series = ny), dimnames = list(origin = as.character(origins), series = colnames(ally)))
for (T in origins) {
  T_chr <- as.character(T)
  train_data <- as.matrix(ally[(T - window_size + 1):T, , drop = FALSE])
  snaive_errors <- train_data[(seasonal_lag + 1):window_size, , drop = FALSE] - train_data[1:(window_size - seasonal_lag), , drop = FALSE]
  d_array[T_chr, ] <- colMeans(snaive_errors^2, na.rm = TRUE)
}
# q^2 arrays have exactly the same dimensions as corresponding error arrays
q2_array <- array(NA_real_, dim = dim(error_array), dimnames = dimnames(error_array))
q2_array_flaic <- array(NA_real_, dim = dim(error_array_flaic), dimnames = dimnames(error_array_flaic))
for (T in origins) {
  T_chr <- as.character(T)
  d_T <- d_array[T_chr, ]
  q2_array[T_chr, , , ] <- sweep(error_array[T_chr, , , ]^2, MARGIN = 2, STATS = d_T, FUN = "/")
  q2_array_flaic[T_chr, , , , , , ] <- sweep(error_array_flaic[T_chr, , , , , , ]^2, MARGIN = 2, STATS = d_T, FUN = "/")
}
# mean squared scaled forecast error arrays
mssfe_array <- array(NA_real_, dim = dim(msfe_array), dimnames = dimnames(msfe_array))
mssfe_array_flaic <- array(NA_real_, dim = dim(msfe_array_flaic), dimnames = dimnames(msfe_array_flaic))
for (method in methods) {
  for (group in names(group_idx)) {
    for (hor in paste0("h", 1:h)) mssfe_array[hor, group, method] <- mean(q2_array[, hor, group_idx[[group]], method], na.rm = TRUE)
    mssfe_array["1-6", group, method] <- mean(q2_array[, paste0("h", 1:6), group_idx[[group]], method], na.rm = TRUE)
    mssfe_array["1-12", group, method] <- mean(q2_array[, paste0("h", 1:h), group_idx[[group]], method], na.rm = TRUE)
  }
}
for (penalty in penalties) {
  for (npc in as.character(p_values)) {
    for (center in centers) {
      for (scale in scales) {
        center_chr <- as.character(center)
        scale_chr <- as.character(scale)
        for (hor in paste0("h", 1:h)) {
          for (group in names(group_idx)) mssfe_array_flaic[hor, group, penalty, npc, center_chr, scale_chr] <- mean(q2_array_flaic[, hor, group_idx[[group]], penalty, npc, center_chr, scale_chr], na.rm = TRUE)
        }
        for (group in names(group_idx)) {
          mssfe_array_flaic["1-6", group, penalty, npc, center_chr, scale_chr] <- mean(q2_array_flaic[, paste0("h", 1:6), group_idx[[group]], penalty, npc, center_chr, scale_chr], na.rm = TRUE)
          mssfe_array_flaic["1-12", group, penalty, npc, center_chr, scale_chr] <- mean(q2_array_flaic[, paste0("h", 1:h), group_idx[[group]], penalty, npc, center_chr, scale_chr], na.rm = TRUE)
        }
      }
    }
  }
}
# scaled PRIAL
prial_scaled_array <- array(NA_real_, dim = dim(mssfe_array), dimnames = dimnames(mssfe_array))
prial_scaled_array_flaic <- array(NA_real_, dim = dim(mssfe_array_flaic), dimnames = dimnames(mssfe_array_flaic))
for (method in methods) {
  for (hor in horizon_names) {
    for (group in names(group_idx)) {
      if (method == "base") prial_scaled_array[hor, group, method] <- 0 else prial_scaled_array[hor, group, method] <- 100 * (mssfe_array[hor, group, method]/mssfe_array[hor, group, "base"] - 1)
    }
  }
}
for (penalty in penalties) {
  for (npc in as.character(p_values)) {
    for (center in centers) {
      for (scale in scales) {
        center_chr <- as.character(center)
        scale_chr <- as.character(scale)
        for (hor in horizon_names) {
          for (group in names(group_idx)) prial_scaled_array_flaic[hor, group, penalty, npc, center_chr, scale_chr] <- 100 * (mssfe_array_flaic[hor, group, penalty, npc, center_chr, scale_chr]/mssfe_array[hor, group, "base"] - 1)
        }
      }
    }
  }
}

#### Option 2: variance-standardised scaling ------------------------------------------------------------------------------------------------------
# denominator = training-window variance of each response series, measured as mean squared deviation from the training-window mean

d_var_array <- array(NA_real_, dim = c(origin = length(origins), series = ny), dimnames = list(origin = as.character(origins), series = colnames(ally)))
for (T in origins) {
  T_chr <- as.character(T)
  train_data <- as.matrix(ally[(T - window_size + 1):T, , drop = FALSE])
  train_data_centered <- scale(train_data, center = TRUE, scale = FALSE)
  d_var_array[T_chr, ] <- colMeans(train_data_centered^2, na.rm = TRUE)
}
q2_var_array <- sweep(error_array^2, MARGIN = c(1, 3), STATS = d_var_array, FUN = "/")
q2_var_array_flaic <- sweep(error_array_flaic^2, MARGIN = c(1, 3), STATS = d_var_array, FUN = "/")
mssfe_var_array <- array(NA_real_, dim = dim(msfe_array), dimnames = dimnames(msfe_array))
mssfe_var_array_flaic <- array(NA_real_, dim = dim(msfe_array_flaic), dimnames = dimnames(msfe_array_flaic))
for (method in methods) {
  for (group in names(group_idx)) {
    for (hor in paste0("h", 1:h)) mssfe_var_array[hor, group, method] <- mean(q2_var_array[, hor, group_idx[[group]], method], na.rm = TRUE)
    mssfe_var_array["1-6", group, method] <- mean(q2_var_array[, paste0("h", 1:6), group_idx[[group]], method], na.rm = TRUE)
    mssfe_var_array["1-12", group, method] <- mean(q2_var_array[, paste0("h", 1:h), group_idx[[group]], method], na.rm = TRUE)
  }
}
for (penalty in penalties) {
  for (npc in as.character(p_values)) {
    for (center in centers) {
      for (scale in scales) {
        center_chr <- as.character(center)
        scale_chr <- as.character(scale)
        for (hor in paste0("h", 1:h)) {
          for (group in names(group_idx)) mssfe_var_array_flaic[hor, group, penalty, npc, center_chr, scale_chr] <- mean(q2_var_array_flaic[, hor, group_idx[[group]], penalty, npc, center_chr, scale_chr], na.rm = TRUE)
        }
        for (group in names(group_idx)) {
          mssfe_var_array_flaic["1-6", group, penalty, npc, center_chr, scale_chr] <- mean(q2_var_array_flaic[, paste0("h", 1:6), group_idx[[group]], penalty, npc, center_chr, scale_chr], na.rm = TRUE)
          mssfe_var_array_flaic["1-12", group, penalty, npc, center_chr, scale_chr] <- mean(q2_var_array_flaic[, paste0("h", 1:h), group_idx[[group]], penalty, npc, center_chr, scale_chr], na.rm = TRUE)
        }
      }
    }
  }
}
prial_var_scaled_array <- array(NA_real_, dim = dim(mssfe_var_array), dimnames = dimnames(mssfe_var_array))
prial_var_scaled_array_flaic <- array(NA_real_, dim = dim(mssfe_var_array_flaic), dimnames = dimnames(mssfe_var_array_flaic))
for (method in methods) {
  for (hor in horizon_names) {
    for (group in names(group_idx)) {
      if (method == "base") prial_var_scaled_array[hor, group, method] <- 0 else prial_var_scaled_array[hor, group, method] <- 100 * (mssfe_var_array[hor, group, method]/mssfe_var_array[hor, group, "base"] - 1)
    }
  }
}
for (penalty in penalties) {
  for (npc in as.character(p_values)) {
    for (center in centers) {
      for (scale in scales) {
        center_chr <- as.character(center)
        scale_chr <- as.character(scale)
        for (hor in horizon_names) {
          for (group in names(group_idx)) prial_var_scaled_array_flaic[hor, group, penalty, npc, center_chr, scale_chr] <- 100 * (mssfe_var_array_flaic[hor, group, penalty, npc, center_chr, scale_chr]/mssfe_var_array[hor, group, "base"] - 1)
        }
      }
    }
  }
}

### Final save ------------------------------------------------------------------------------------------------------------------------------------

saveRDS(list(actual_array = actual_array, error_array = error_array, error_array_flaic = error_array_flaic, error_array_flap = error_array_flap,
             msfe_array = msfe_array, msfe_array_flaic = msfe_array_flaic, msfe_array_flap = msfe_array_flap,
             prial_array = prial_array, prial_array_flaic = prial_array_flaic, prial_array_flap = prial_array_flap,
             d_array = d_array, q2_array = q2_array, q2_array_flaic = q2_array_flaic, mssfe_array = mssfe_array, mssfe_array_flaic = mssfe_array_flaic,
             prial_scaled_array = prial_scaled_array, prial_scaled_array_flaic = prial_scaled_array_flaic,
             d_var_array = d_var_array, q2_var_array = q2_var_array, q2_var_array_flaic = q2_var_array_flaic, mssfe_var_array = mssfe_var_array,
             mssfe_var_array_flaic = mssfe_var_array_flaic, prial_var_scaled_array = prial_var_scaled_array,
             prial_var_scaled_array_flaic = prial_var_scaled_array_flaic), "rds/flaic_05_evaluation_complete.rds")
message("05_evaluate complete: raw evaluation, FLAP benchmark evaluation and scaled-error robustness saved.")

# save file for global environment after the complete 05
# Full-session .RData omitted: saved RDS outputs are sufficient.