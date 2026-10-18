# Forecasting with linear augmented information combination (FLAIC)
# Appendix robustness: PRIAL using scaled forecast errors
# Produces Average-hierarchy tables for seasonal-naive-scaled and variance-standardised errors.

library(tidyverse)
dir.create("rds", showWarnings = FALSE)

data_setup <- readRDS("rds/flaic_01_data_setup.rds")
list2env(data_setup, envir = .GlobalEnv)
rm(data_setup)

evaluation_results <- readRDS("rds/flaic_05_evaluation_complete.rds")
list2env(evaluation_results, envir = .GlobalEnv)
rm(evaluation_results)

stopifnot(exists("ally"), exists("origins"), exists("window_size"),
          exists("group_idx"), exists("error_array"), exists("error_array_flaic"))

dn_flaic <- dimnames(error_array_flaic)
p_names <- if (!is.null(names(dn_flaic)) && "npc" %in% names(dn_flaic)) dn_flaic[["npc"]] else dn_flaic[[5]]
analysis_p <- as.integer(p_names)
analysis_p <- analysis_p[!is.na(analysis_p)]
stopifnot(identical(analysis_p, c(1L, 2L, 3L, 5L, 10L, 20L, 40L, 77L)))

analysis_h <- 1:12
horizon_names <- paste0("h", analysis_h)
value_names <- c(horizon_names, "h16", "h112")
avg_idx <- group_idx[["Average"]]
stopifnot(length(avg_idx) == ncol(ally))

method_map <- c(
  BU = "bu", OLS = "ols", `WLS$_s$` = "wls_s", `WLS$_v$` = "wls_v",
  `MinT$_{shr}$` = "mint_shr", `IComb-S` = "icombs_xy1", `IComb-L` = "icombl_xy1"
)
dn_error <- dimnames(error_array)
method_names <- if (!is.null(names(dn_error)) && "method" %in% names(dn_error)) dn_error[["method"]] else dn_error[[4]]
stopifnot(all(unname(method_map) %in% method_names))

pca_specs <- tibble::tribble(
  ~preprocessing, ~center, ~scale,
  "Centred, unscaled", TRUE, FALSE,
  "Centred, scaled", TRUE, TRUE,
  "Uncentred, unscaled", FALSE, FALSE,
  "Uncentred, scaled", FALSE, TRUE
)

make_scaling_denominator <- function(type = c("seasonal_naive", "variance"), seasonal_lag = 12L) {
  type <- match.arg(type)
  d <- matrix(NA_real_, length(origins), ncol(ally),
              dimnames = list(origin = as.character(origins), series = colnames(ally)))

  for (T in origins) {
    train <- as.matrix(ally[(T - window_size + 1L):T, , drop = FALSE])
    if (type == "seasonal_naive") {
      e <- train[(seasonal_lag + 1L):window_size, , drop = FALSE] -
        train[1L:(window_size - seasonal_lag), , drop = FALSE]
      d[as.character(T), ] <- colMeans(e^2, na.rm = TRUE)
    } else {
      x <- scale(train, center = TRUE, scale = FALSE)
      d[as.character(T), ] <- colMeans(x^2, na.rm = TRUE)
    }
  }

  if (any(!is.finite(d)) || any(d <= 0, na.rm = TRUE))
    stop("Scaling denominator contains non-finite or non-positive values.")
  d
}

scaled_loss_comparison <- function(method, denominator) {
  out <- setNames(numeric(length(value_names)), value_names)

  for (hh in analysis_h) {
    hn <- paste0("h", hh)
    e <- error_array[, hn, avg_idx, method, drop = FALSE]
    q2 <- sweep(e^2, c(1, 3), denominator[, avg_idx, drop = FALSE], "/")
    out[hn] <- mean(q2, na.rm = TRUE)
  }

  e <- error_array[, paste0("h", 1:6), avg_idx, method, drop = FALSE]
  out["h16"] <- mean(sweep(e^2, c(1, 3), denominator[, avg_idx, drop = FALSE], "/"), na.rm = TRUE)

  e <- error_array[, paste0("h", 1:12), avg_idx, method, drop = FALSE]
  out["h112"] <- mean(sweep(e^2, c(1, 3), denominator[, avg_idx, drop = FALSE], "/"), na.rm = TRUE)
  out
}

scaled_loss_flaic <- function(penalty, p, center, scale, denominator) {
  pc <- as.character(p); cc <- as.character(center); ss <- as.character(scale)
  out <- setNames(numeric(length(value_names)), value_names)

  for (hh in analysis_h) {
    hn <- paste0("h", hh)
    e <- error_array_flaic[, hn, avg_idx, penalty, pc, cc, ss, drop = FALSE]
    q2 <- sweep(e^2, c(1, 3), denominator[, avg_idx, drop = FALSE], "/")
    out[hn] <- mean(q2, na.rm = TRUE)
  }

  e <- error_array_flaic[, paste0("h", 1:6), avg_idx, penalty, pc, cc, ss, drop = FALSE]
  out["h16"] <- mean(sweep(e^2, c(1, 3), denominator[, avg_idx, drop = FALSE], "/"), na.rm = TRUE)

  e <- error_array_flaic[, paste0("h", 1:12), avg_idx, penalty, pc, cc, ss, drop = FALSE]
  out["h112"] <- mean(sweep(e^2, c(1, 3), denominator[, avg_idx, drop = FALSE], "/"), na.rm = TRUE)
  out
}

make_scaled_prial_average <- function(denominator) {
  base_loss <- scaled_loss_comparison("base", denominator)

  conventional <- purrr::map_dfr(c("BU", "OLS", "WLS$_s$", "WLS$_v$", "MinT$_{shr}$"), function(label) {
    loss <- scaled_loss_comparison(unname(method_map[label]), denominator)
    tibble::tibble(Family = "Conventional reconciliation", Method = label,
                   preprocessing = "", p = NA_integer_,
                   !!!as.list(100 * (loss / base_loss - 1)))
  })

  icomb_s <- scaled_loss_comparison(unname(method_map["IComb-S"]), denominator)
  icomb_l <- scaled_loss_comparison(unname(method_map["IComb-L"]), denominator)

  ridge <- dplyr::bind_rows(
    tibble::tibble(Family = "Ridge", Method = "IComb-S", preprocessing = "", p = NA_integer_,
                   !!!as.list(100 * (icomb_s / base_loss - 1))),
    purrr::pmap_dfr(pca_specs, function(preprocessing, center, scale) {
      purrr::map_dfr(analysis_p, function(p) {
        loss <- scaled_loss_flaic("ridge", p, center, scale, denominator)
        tibble::tibble(Family = "Ridge", Method = "FLAIC-S", preprocessing = preprocessing, p = p,
                       !!!as.list(100 * (loss / base_loss - 1)))
      })
    })
  )

  lasso <- dplyr::bind_rows(
    tibble::tibble(Family = "Lasso", Method = "IComb-L", preprocessing = "", p = NA_integer_,
                   !!!as.list(100 * (icomb_l / base_loss - 1))),
    purrr::pmap_dfr(pca_specs, function(preprocessing, center, scale) {
      purrr::map_dfr(analysis_p, function(p) {
        loss <- scaled_loss_flaic("lasso", p, center, scale, denominator)
        tibble::tibble(Family = "Lasso", Method = "FLAIC-L", preprocessing = preprocessing, p = p,
                       !!!as.list(100 * (loss / base_loss - 1)))
      })
    })
  )

  out <- dplyr::bind_rows(conventional, ridge, lasso)
  stopifnot(nrow(out) == 71L)
  out
}

# Seasonal-naive-scaled PRIAL
d_snaive <- make_scaling_denominator("seasonal_naive")
prial_scaled_snaive_average <- make_scaled_prial_average(d_snaive)

# Variance-standardised PRIAL
d_variance <- make_scaling_denominator("variance")
prial_scaled_variance_average <- make_scaled_prial_average(d_variance)

saveRDS(prial_scaled_snaive_average, "rds/table_prial_scaled_snaive_average.rds")
saveRDS(prial_scaled_variance_average, "rds/table_prial_scaled_variance_average.rds")
saveRDS(
  list(
    seasonal_naive_denominator = d_snaive,
    variance_denominator = d_variance,
    seasonal_naive_table = prial_scaled_snaive_average,
    variance_standardised_table = prial_scaled_variance_average
  ),
  "rds/prial_scaled_error_objects.rds"
)

message("Saved scaled-error PRIAL appendix tables.")

#### Delta-PRIAL heatmaps (Average hierarchy only) -------------------------------------------------------------------------

dir.create("figures", showWarnings = FALSE)

make_delta_prial_long <- function(prial_tbl) {
  value_cols <- paste0("h", 1:12)
  
  icomb_s <- prial_tbl %>% dplyr::filter(Method == "IComb-S") %>% dplyr::slice(1)
  icomb_l <- prial_tbl %>% dplyr::filter(Method == "IComb-L") %>% dplyr::slice(1)
  
  ridge_df <- prial_tbl %>%
    dplyr::filter(Method == "FLAIC-S") %>%
    dplyr::mutate(
      Estimator = "Panel A.\nFLAIC-S",
      dplyr::across(dplyr::all_of(value_cols), ~ icomb_s[[dplyr::cur_column()]] - .x)
    )
  
  lasso_df <- prial_tbl %>%
    dplyr::filter(Method == "FLAIC-L") %>%
    dplyr::mutate(
      Estimator = "Panel B.\nFLAIC-L",
      dplyr::across(dplyr::all_of(value_cols), ~ icomb_l[[dplyr::cur_column()]] - .x)
    )
  
  dplyr::bind_rows(ridge_df, lasso_df) %>%
    dplyr::select(Estimator, preprocessing, p, dplyr::all_of(value_cols)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(value_cols),
      names_to = "h",
      values_to = "delta_prial"
    ) %>%
    dplyr::mutate(
      preprocessing = factor(preprocessing, levels = pca_specs$preprocessing),
      p = factor(p, levels = analysis_p),
      h = factor(sub("^h", "", h), levels = as.character(1:12))
    )
}

delta_snaive_long <- make_delta_prial_long(prial_scaled_snaive_average)
delta_variance_long <- make_delta_prial_long(prial_scaled_variance_average)

delta_lim <- max(
  abs(c(delta_snaive_long$delta_prial, delta_variance_long$delta_prial)),
  na.rm = TRUE
)
delta_lim <- ceiling(delta_lim * 2) / 2

make_delta_prial_heatmap <- function(plot_df, figure_title, show_title = TRUE,
                                     fill_limit = delta_lim,
                                     label_size = 3.1,
                                     label_bold = FALSE) {
  
  plot_df <- plot_df %>%
    dplyr::mutate(
      label = ifelse(abs(delta_prial) < 0.1,
                     sprintf("%.2f", delta_prial),
                     sprintf("%.1f", delta_prial)),
      text_col = ifelse(abs(delta_prial) >= 0.55 * fill_limit, "white", "black")
    )
  
  ggplot(plot_df, aes(x = p, y = h, fill = delta_prial)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    geom_text(aes(label = label, colour = text_col),
              size = label_size,
              fontface = if (label_bold) "bold" else "plain") +
    facet_grid(Estimator ~ preprocessing, switch = "y") +
    scale_fill_gradient2(
      low = "#B2182B", mid = "white", high = "#2166AC",
      midpoint = 0, limits = c(-fill_limit, fill_limit),
      name = expression(Delta * "PRIAL (pp)")
    ) +
    scale_colour_identity() +
    labs(
      title = if (show_title) figure_title else NULL,
      x = "Number of principal components (p)",
      y = "Forecast horizon (h)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold", size = 10),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      plot.title = element_text(face = "plain", size = 15, hjust = 0)
    )
}

# preview versions (optional)
delta_snaive_preview <- make_delta_prial_heatmap(
  delta_snaive_long,
  "PRIAL of FLAIC relative to matching IComb: Average",
  show_title = TRUE
)

delta_variance_preview <- make_delta_prial_heatmap(
  delta_variance_long,
  "PRIAL of FLAIC relative to matching IComb: Average",
  show_title = TRUE
)

# report versions (no internal title)
delta_snaive_fig <- make_delta_prial_heatmap(
  delta_snaive_long,
  "PRIAL of FLAIC relative to matching IComb: Average",
  show_title = FALSE
)

delta_variance_fig <- make_delta_prial_heatmap(
  delta_variance_long,
  "PRIAL of FLAIC relative to matching IComb: Average",
  show_title = FALSE
)

ggsave(
  "figures/figure_prial_delta_scaled_snaive_average.png",
  delta_snaive_fig, width = 13.5, height = 9.2, dpi = 300
)

ggsave(
  "figures/figure_prial_delta_scaled_variance_average.png",
  delta_variance_fig, width = 13.5, height = 9.2, dpi = 300
)

saveRDS(
  list(
    seasonal_naive = delta_snaive_long,
    variance_standardised = delta_variance_long,
    fill_limit = delta_lim
  ),
  "rds/prial_scaled_error_heatmaps.rds"
)

message("Saved scaled-error PRIAL appendix tables and delta-PRIAL heatmaps")