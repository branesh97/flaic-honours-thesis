# Forecasting with linear augmented information combination (FLAIC)
# Analysis workbench: exploratory results, diagnostics, MCB, lasso selection and FLAP tables
# This file is deliberately separate from the final 06/07/08 scripts. PC-only IComb is not included.

library(tidyverse)
library(patchwork)

dir.create("rds", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# load data/settings and completed estimation/evaluation results
data_setup <- readRDS("rds/flaic_01_data_setup.rds")
list2env(data_setup, envir = .GlobalEnv)
rm(data_setup)
forecast_results <- readRDS("rds/flaic_03_forecasts_complete.rds")
list2env(forecast_results, envir = .GlobalEnv)
rm(forecast_results)
flap_results <- readRDS("rds/flaic_04_flap_complete.rds")
list2env(flap_results, envir = .GlobalEnv)
rm(flap_results)
evaluation_results <- readRDS("rds/flaic_05_evaluation_complete.rds")
list2env(evaluation_results, envir = .GlobalEnv)
rm(evaluation_results)

### Common analysis settings ---------------------------------------------------------------------------------------------------------------------

analysis_groups <- c("Average", "Australia", "States", "Zones", "Regions")
analysis_levels <- c("Australia", "States", "Zones", "Regions")
analysis_p <- as.integer(dimnames(prial_array_flaic)$npc)
analysis_h <- seq_len(h)
cs_levels <- c("Centred, unscaled", "Centred, scaled", "Uncentred, unscaled", "Uncentred, scaled")
scale_cols <- c("Unscaled" = "#0072B2", "Scaled" = "#D55E00")
center_ltypes <- c("Centred" = "solid", "Uncentred" = "longdash")
level_cols <- c(Australia = "#000000", States = "#D55E00", Zones = "#0072B2", Regions = "#009E73")

get_icomb_method <- function(penalty) ifelse(penalty == "ridge", "icombs_xy1", "icombl_xy1")
get_estimator_label <- function(penalty) ifelse(penalty == "ridge", "FLAIC-S", "FLAIC-L")
get_cs_label <- function(center, scale) paste0("C=", ifelse(center, "T", "F"), ", S=", ifelse(scale, "T", "F"))
get_preprocessing_label <- function(center, scale) paste0(ifelse(center, "Centred", "Uncentred"), ", ", ifelse(scale, "scaled", "unscaled"))
format_heat_delta <- function(x) ifelse(abs(x) < 0.05, sprintf("%.2f", ifelse(abs(x) < 0.005, 0, x)), sprintf("%.1f", x))

stopifnot(identical(as.character(analysis_p), dimnames(prial_array_flaic)$npc), identical(as.character(analysis_p), dimnames(fcst_array_flaic)$npc),
          identical(as.character(analysis_p), dimnames(prial_array_flap)$p))

### Analysis 1: origin-level FLAIC versus matching IComb heatmaps --------------------------------------------------------------------------------

# At each rolling origin, compare average squared-error loss across the requested hierarchy group.
# Positive gain means FLAIC has lower origin-level loss than the matching IComb estimator.
as_origin_series_matrix <- function(x, series_idx) {
  matrix(as.numeric(x), nrow = length(origins), ncol = length(series_idx),
         dimnames = list(origin = as.character(origins), series = colnames(ally)[series_idx]))
}

make_origin_summary_one <- function(group, hor, p, center, scale, penalty) {
  h_chr <- paste0("h", hor)
  p_chr <- as.character(p)
  c_chr <- as.character(center)
  s_chr <- as.character(scale)
  idx <- group_idx[[group]]
  icomb_method <- get_icomb_method(penalty)
  err_icomb <- as_origin_series_matrix(error_array[, h_chr, idx, icomb_method, drop = FALSE], idx)
  err_flaic <- as_origin_series_matrix(error_array_flaic[, h_chr, idx, penalty, p_chr, c_chr, s_chr, drop = FALSE], idx)
  loss_icomb <- rowMeans(err_icomb^2, na.rm = FALSE)
  loss_flaic <- rowMeans(err_flaic^2, na.rm = FALSE)
  valid <- is.finite(loss_icomb) & is.finite(loss_flaic)
  gain <- loss_icomb[valid] - loss_flaic[valid]
  tibble::tibble(Group = group, horizon = hor, p = p, center = center, scale = scale, penalty = penalty,
                 n_origins = sum(valid), pct_flaic_better = 100 * mean(gain > 0), pct_icomb_better = 100 * mean(gain < 0),
                 pct_tie = 100 * mean(gain == 0), mean_gain = mean(gain), median_gain = median(gain))
}

origin_heat_all <- tidyr::expand_grid(Group = analysis_groups, horizon = analysis_h, p = analysis_p, center = c(FALSE, TRUE), scale = c(FALSE, TRUE), penalty = c("ridge", "lasso")) %>%
  purrr::pmap_dfr(function(Group, horizon, p, center, scale, penalty) make_origin_summary_one(Group, horizon, p, center, scale, penalty)) %>%
  dplyr::mutate(Estimator = factor(ifelse(penalty == "ridge", "Panel A.\nFLAIC-S", "Panel B.\nFLAIC-L"), levels = c("Panel A.\nFLAIC-S", "Panel B.\nFLAIC-L")),
                CS = factor(get_preprocessing_label(center, scale), levels = cs_levels),
                p_factor = factor(p, levels = analysis_p),
                Horizon = factor(horizon, levels = rev(analysis_h)))

stopifnot(all(abs(origin_heat_all$pct_flaic_better + origin_heat_all$pct_icomb_better + origin_heat_all$pct_tie - 100) < 1e-10))
origin_dev <- min(50, max(10, ceiling(max(abs(origin_heat_all$pct_flaic_better - 50), na.rm = TRUE)/5) * 5))
origin_limits <- c(50 - origin_dev, 50 + origin_dev)

make_origin_heatmap <- function(group_name, show_title = TRUE) {
  plot_df <- origin_heat_all %>% dplyr::filter(Group == group_name) %>%
    dplyr::mutate(label = sprintf("%.0f", pct_flaic_better),
                  text_col = ifelse(abs(pct_flaic_better - 50) >= 0.60 * origin_dev, "white", "black"))
  
  ggplot(plot_df, aes(x = p_factor, y = Horizon, fill = pct_flaic_better)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = label, colour = text_col), size = 2.5) +
    facet_grid(rows = vars(Estimator), cols = vars(CS), switch = "y") +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 50, limits = origin_limits, oob = scales::squish,
                         breaks = pretty(origin_limits, n = 5), labels = function(x) paste0(x, "%"), name = "FLAIC wins\n(% rolling forecast origins)") +
    scale_colour_identity() +
    labs(title = if (show_title) paste0("Forecast origin-level FLAIC performance relative to IComb: ", group_name) else NULL,
         x = "Number of principal components (p)", y = "Forecast horizon (h)") +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey90"), strip.placement = "outside",
          strip.text.x = element_text(face = "bold", size = 8.5), strip.text.y.left = element_text(angle = 0, face = "bold", size = 8.5),
          legend.position = "bottom", plot.title = element_text(face = "plain", size = 12, hjust = 0),
          axis.title = element_text(size = 9), axis.text = element_text(size = 8))
}

origin_heatmaps <- stats::setNames(lapply(analysis_groups, function(g) make_origin_heatmap(g, show_title = FALSE)), analysis_groups)
for (g in analysis_groups) ggsave(paste0("figures/figure_origin_win_heatmap_", tolower(g), ".png"), origin_heatmaps[[g]], width = 13.5, height = 9.2, dpi = 300)

### Analysis 2a: incremental PRIAL heatmaps ------------------------------------------------------------------------------------------------------

# Delta PRIAL = PRIAL(IComb) - PRIAL(FLAIC); positive values indicate improvement from component augmentation.
make_delta_prial_one <- function(group, hor, p, center, scale, penalty) {
  h_chr <- paste0("h", hor)
  icomb_method <- get_icomb_method(penalty)
  prial_icomb <- as.numeric(prial_array[h_chr, group, icomb_method])
  prial_flaic <- as.numeric(prial_array_flaic[h_chr, group, penalty, as.character(p), as.character(center), as.character(scale)])
  tibble::tibble(Group = group, horizon = hor, p = p, center = center, scale = scale, penalty = penalty,
                 prial_icomb = prial_icomb, prial_flaic = prial_flaic, delta_prial = prial_icomb - prial_flaic)
}

delta_prial_heat_all <- tidyr::expand_grid(Group = analysis_groups, horizon = analysis_h, p = analysis_p, center = c(FALSE, TRUE), scale = c(FALSE, TRUE), penalty = c("ridge", "lasso")) %>%
  purrr::pmap_dfr(function(Group, horizon, p, center, scale, penalty) make_delta_prial_one(Group, horizon, p, center, scale, penalty)) %>%
  dplyr::mutate(Estimator = factor(ifelse(penalty == "ridge", "Panel A.\nFLAIC-S", "Panel B.\nFLAIC-L"), levels = c("Panel A.\nFLAIC-S", "Panel B.\nFLAIC-L")),
                CS = factor(get_preprocessing_label(center, scale), levels = cs_levels),
                p_factor = factor(p, levels = analysis_p),
                Horizon = factor(horizon, levels = rev(analysis_h)))

stopifnot(all(is.finite(delta_prial_heat_all$delta_prial)))
prial_cap <- ceiling(as.numeric(stats::quantile(abs(delta_prial_heat_all$delta_prial), 0.95, na.rm = TRUE)))
prial_cap <- max(prial_cap, 1)

make_delta_prial_heatmap <- function(group_name, show_title = TRUE) {
  plot_df <- delta_prial_heat_all %>% dplyr::filter(Group == group_name) %>%
    dplyr::mutate(label = format_heat_delta(delta_prial),
                  text_col = ifelse(abs(delta_prial) >= 0.60 * prial_cap, "white", "black"))
  
  ggplot(plot_df, aes(x = p_factor, y = Horizon, fill = delta_prial)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = label, colour = text_col), size = 2.45) +
    facet_grid(rows = vars(Estimator), cols = vars(CS), switch = "y") +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, limits = c(-prial_cap, prial_cap), oob = scales::squish,
                         breaks = seq(-prial_cap, prial_cap, length.out = 5), name = expression(Delta*"PRIAL (pp)")) +
    scale_colour_identity() +
    labs(title = if (show_title) paste0("PRIAL of FLAIC relative to IComb: ", group_name) else NULL,
         x = "Number of principal components (p)", y = "Forecast horizon (h)") +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey90"), strip.placement = "outside",
          strip.text.x = element_text(face = "bold", size = 8.5), strip.text.y.left = element_text(angle = 0, face = "bold", size = 8.5),
          legend.position = "bottom", plot.title = element_text(face = "plain", size = 12, hjust = 0),
          axis.title = element_text(size = 9), axis.text = element_text(size = 8))
}

prial_heatmaps <- stats::setNames(lapply(analysis_groups, function(g) make_delta_prial_heatmap(g, show_title = FALSE)), analysis_groups)
for (g in analysis_groups) ggsave(paste0("figures/figure_prial_heatmap_", tolower(g), ".png"), prial_heatmaps[[g]], width = 13.5, height = 9.2, dpi = 300)

### Analysis 2b: heatmaps across hierarchy levels for a fixed PCA preprocessing specification -----------------------------------------------

make_hierarchy_origin_heatmap <- function(center = TRUE, scale = FALSE, show_title = TRUE, label_size = 2.8, label_bold = TRUE) {
  spec_label <- get_preprocessing_label(center, scale)
  plot_df <- origin_heat_all %>%
    dplyr::filter(Group %in% analysis_levels, .data$center == .env$center, .data$scale == .env$scale) %>%
    dplyr::mutate(Group = factor(Group, levels = analysis_levels),
                  label = sprintf("%.0f", pct_flaic_better),
                  text_col = ifelse(abs(pct_flaic_better - 50) >= 0.60 * origin_dev, "white", "black"))
  
  ggplot(plot_df, aes(x = p_factor, y = Horizon, fill = pct_flaic_better)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = label, colour = text_col), size = label_size, fontface = ifelse(label_bold, "bold", "plain")) +
    facet_grid(rows = vars(Estimator), cols = vars(Group), switch = "y") +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 50, limits = origin_limits, oob = scales::squish,
                         breaks = pretty(origin_limits, n = 5), labels = function(x) paste0(x, "%"), name = "FLAIC wins\n(% of rolling forecast origins)") +
    scale_colour_identity() +
    labs(title = if (show_title) "Forecast origin-level FLAIC performance relative to matching IComb by hierarchy level" else NULL,
         subtitle = if (show_title) paste0(spec_label, " PCA") else NULL,
         x = "Number of principal components (p)", y = "Forecast horizon (h)") +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey90"), strip.placement = "outside",
          strip.text.x = element_text(face = "bold", size = 8.5), strip.text.y.left = element_text(angle = 0, face = "bold", size = 8.5, lineheight = 0.95),
          legend.position = "bottom", plot.title = element_text(face = "plain", size = 12, hjust = 0),
          plot.subtitle = element_text(size = 9.5, hjust = 0), axis.title = element_text(size = 9), axis.text = element_text(size = 8))
}

make_hierarchy_prial_heatmap <- function(center = TRUE, scale = FALSE, show_title = TRUE, label_size = 2.8, label_bold = TRUE) {
  spec_label <- get_preprocessing_label(center, scale)
  plot_df <- delta_prial_heat_all %>%
    dplyr::filter(Group %in% analysis_levels, .data$center == .env$center, .data$scale == .env$scale) %>%
    dplyr::mutate(Group = factor(Group, levels = analysis_levels),
                  label = format_heat_delta(delta_prial),
                  text_col = ifelse(abs(delta_prial) >= 0.60 * prial_cap, "white", "black"))
  
  ggplot(plot_df, aes(x = p_factor, y = Horizon, fill = delta_prial)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = label, colour = text_col), size = label_size, fontface = ifelse(label_bold, "bold", "plain")) +
    facet_grid(rows = vars(Estimator), cols = vars(Group), switch = "y") +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, limits = c(-prial_cap, prial_cap), oob = scales::squish,
                         breaks = seq(-prial_cap, prial_cap, length.out = 5), name = expression(Delta*"PRIAL (pp)")) +
    scale_colour_identity() +
    labs(title = if (show_title) "PRIAL of FLAIC relative to matching IComb by hierarchy level" else NULL,
         subtitle = if (show_title) paste0(spec_label, " PCA") else NULL,
         x = "Number of principal components (p)", y = "Forecast horizon (h)") +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey90"), strip.placement = "outside",
          strip.text.x = element_text(face = "bold", size = 8.5), strip.text.y.left = element_text(angle = 0, face = "bold", size = 8.5, lineheight = 0.95),
          legend.position = "bottom", plot.title = element_text(face = "plain", size = 12, hjust = 0),
          plot.subtitle = element_text(size = 9.5, hjust = 0), axis.title = element_text(size = 9), axis.text = element_text(size = 8))
}

spec_grid <- tibble::tribble(
  ~center, ~scale, ~spec_name,
  TRUE,   FALSE,  "centred_unscaled",
  TRUE,   TRUE,   "centred_scaled",
  FALSE,  FALSE,  "uncentred_unscaled",
  FALSE,  TRUE,   "uncentred_scaled"
)

hierarchy_origin_heatmaps_ <- stats::setNames(vector("list", nrow(spec_grid)), spec_grid$spec_name)
hierarchy_prial_heatmaps_ <- stats::setNames(vector("list", nrow(spec_grid)), spec_grid$spec_name)

for (i in seq_len(nrow(spec_grid))) {
  hierarchy_origin_heatmaps_[[spec_grid$spec_name[i]]] <- make_hierarchy_origin_heatmap(center = spec_grid$center[i], scale = spec_grid$scale[i], show_title = FALSE, label_size = 2.8, label_bold = FALSE)
  hierarchy_prial_heatmaps_[[spec_grid$spec_name[i]]] <- make_hierarchy_prial_heatmap(center = spec_grid$center[i], scale = spec_grid$scale[i], show_title = FALSE, label_size = 2.8, label_bold = FALSE)
}

for (nm in names(hierarchy_origin_heatmaps_)) ggsave(file.path("figures", paste0("figure_hierarchy_origin_heatmap_", nm, ".png")), hierarchy_origin_heatmaps_[[nm]], width = 13.5, height = 9.2, dpi = 300)
for (nm in names(hierarchy_prial_heatmaps_)) ggsave(file.path("figures", paste0("figure_hierarchy_prial_heatmap_", nm, ".png")), hierarchy_prial_heatmaps_[[nm]], width = 13.5, height = 9.2, dpi = 300)

hierarchy_origin_heatmap_preview <- make_hierarchy_origin_heatmap(center = TRUE, scale = FALSE, show_title = TRUE, label_size = 2.8, label_bold = FALSE)
hierarchy_prial_heatmap_preview <- make_hierarchy_prial_heatmap(center = TRUE, scale = FALSE, show_title = TRUE, label_size = 2.8, label_bold = FALSE)

hierarchy_origin_heatmap_preview
hierarchy_prial_heatmap_preview


make_hierarchy_prial_heatmap(center = TRUE, scale = FALSE, show_title = TRUE, label_size = 2.8, label_bold = FALSE)
make_hierarchy_prial_heatmap(center = TRUE, scale = TRUE, show_title = TRUE, label_size = 2.8, label_bold = FALSE)
make_hierarchy_prial_heatmap(center = FALSE, scale = FALSE, show_title = TRUE, label_size = 2.8, label_bold = FALSE)
make_hierarchy_prial_heatmap(center = FALSE, scale = TRUE, show_title = TRUE, label_size = 2.8, label_bold = FALSE)

### Analysis 3: forecast origin-level decomposition of Delta PRIAL -----------------------------------------------------------------------------------------------

# The histogram retained here shows origin-level contributions to the aggregate incremental PRIAL.
# The scatter plot compares how often FLAIC wins with the magnitude-weighted aggregate PRIAL gain.
make_origin_case <- function(hor, p, center, scale, penalty = c("ridge", "lasso"), group = "Average") {
  penalty <- match.arg(penalty)
  h_chr <- paste0("h", hor); p_chr <- as.character(p); c_chr <- as.character(center); s_chr <- as.character(scale)
  idx <- group_idx[[group]]; icomb <- get_icomb_method(penalty)
  as_mat <- function(x) matrix(as.numeric(x), nrow = length(origins), ncol = length(idx))
  e_base <- as_mat(error_array[, h_chr, idx, "base", drop = FALSE])
  e_icomb <- as_mat(error_array[, h_chr, idx, icomb, drop = FALSE])
  e_flaic <- as_mat(error_array_flaic[, h_chr, idx, penalty, p_chr, c_chr, s_chr, drop = FALSE])
  loss_base <- rowMeans(e_base^2); loss_icomb <- rowMeans(e_icomb^2); loss_flaic <- rowMeans(e_flaic^2)
  valid <- is.finite(loss_base) & is.finite(loss_icomb) & is.finite(loss_flaic)
  base_msfe <- as.numeric(msfe_array[h_chr, group, "base"])
  delta_prial <- as.numeric(prial_array[h_chr, group, icomb] - prial_array_flaic[h_chr, group, penalty, p_chr, c_chr, s_chr])
  out <- tibble::tibble(origin = origins[valid], loss_base = loss_base[valid], loss_icomb = loss_icomb[valid], loss_flaic = loss_flaic[valid]) %>%
    dplyr::mutate(d = loss_icomb - loss_flaic, flaic_better = d > 0,
                  prial_contribution = 100 * d/base_msfe,
                  additive_contribution = prial_contribution/dplyr::n(),
                  relative_gain = 100 * d/loss_base)
  stopifnot(abs(mean(out$prial_contribution) - delta_prial) < 1e-8, abs(sum(out$additive_contribution) - delta_prial) < 1e-8)
  attr(out, "delta_prial") <- delta_prial
  out
}

summarise_origin_case <- function(x) {
  tibble::tibble(n_origins = nrow(x), delta_prial = attr(x, "delta_prial"), win_pct = 100 * mean(x$flaic_better),
                 mean_contribution = mean(x$prial_contribution), median_contribution = median(x$prial_contribution),
                 q05_contribution = quantile(x$prial_contribution, 0.05), q95_contribution = quantile(x$prial_contribution, 0.95),
                 mean_relative_gain = mean(x$relative_gain), median_relative_gain = median(x$relative_gain),
                 q05_relative_gain = quantile(x$relative_gain, 0.05), q95_relative_gain = quantile(x$relative_gain, 0.95))
}

plot_origin_contribution <- function(x, hor, p, center, scale, penalty, group = "Average", show_title = TRUE) {
  delta <- attr(x, "delta_prial"); win <- 100 * mean(x$flaic_better); med <- median(x$prial_contribution)
  spec <- get_preprocessing_label(center, scale)
  
  ggplot(x, aes(x = prial_contribution)) +
    geom_histogram(bins = 25, fill = "grey70", colour = "white") +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.6) +
    geom_vline(xintercept = mean(x$prial_contribution), linewidth = 0.7) +
    geom_vline(xintercept = med, linetype = "dashed", linewidth = 0.7) +
    labs(
      title = if (show_title) {
        paste0(
          "Origin-level normalised loss differences: ",
          get_estimator_label(penalty),
          " relative to matching IComb (", group, ")"
        )
      } else NULL,
      subtitle = bquote(.(spec)*" PCA | p = "*.(p)*" | h = "*.(hor)*" | "*Delta*"PRIAL = "*.(sprintf("%.2f", delta))*" pp | FLAIC wins = "*.(sprintf("%.1f", win))*"%"),
      x = "Origin-level normalised loss difference (pp)",
      y = "Number of rolling forecast origins"
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "plain", size = 12),
          plot.subtitle = element_text(size = 9.5))
}

case_grid_average <- delta_prial_heat_all %>% dplyr::filter(Group == "Average") %>%
  dplyr::select(horizon, p, center, scale, penalty, delta_prial) %>%
  dplyr::left_join(origin_heat_all %>% dplyr::filter(Group == "Average") %>% dplyr::select(horizon, p, center, scale, penalty, pct_flaic_better),
                   by = c("horizon", "p", "center", "scale", "penalty")) %>%
  dplyr::mutate(Estimator = factor(ifelse(penalty == "ridge", "FLAIC-S", "FLAIC-L"), levels = c("FLAIC-S", "FLAIC-L")),
                Scaling = factor(ifelse(scale, "Scaled", "Unscaled"), levels = c("Unscaled", "Scaled")),
                Centring = factor(ifelse(center, "Centred", "Uncentred"), levels = c("Centred", "Uncentred")),
                Quadrant = dplyr::case_when(pct_flaic_better > 50 & delta_prial > 0 ~ "Frequent wins, positive PRIAL",
                                            pct_flaic_better <= 50 & delta_prial > 0 ~ "Less-frequent wins, positive PRIAL",
                                            pct_flaic_better > 50 & delta_prial <= 0 ~ "Frequent wins, negative PRIAL",
                                            TRUE ~ "Less-frequent wins, negative PRIAL"))


# Agreement between aggregate Delta PRIAL and origin-level win frequency
origin_agreement <- case_grid_average %>%
  dplyr::mutate(
    outcome = dplyr::case_when(
      pct_flaic_better > 50 & delta_prial > 0 ~ "Agree: FLAIC favoured",
      pct_flaic_better < 50 & delta_prial < 0 ~ "Agree: IComb favoured",
      pct_flaic_better > 50 & delta_prial < 0 ~ "Disagree: FLAIC wins majority",
      pct_flaic_better < 50 & delta_prial > 0 ~ "Disagree: IComb wins majority",
      TRUE ~ "Tie"
    ),
    agreement = dplyr::case_when(
      startsWith(outcome, "Agree:") ~ "Agreement",
      startsWith(outcome, "Disagree:") ~ "Disagreement",
      TRUE ~ "Tie"
    )
  )

# Detailed counts showing which type of agreement/disagreement occurs
origin_agreement %>%
  dplyr::count(Estimator, outcome, name = "cells") %>%
  dplyr::group_by(Estimator) %>%
  dplyr::mutate(percent = round(100 * cells / sum(cells), 1)) %>%
  print(n = Inf)

# Concise summary
agreement_summary <- origin_agreement %>%
  dplyr::group_by(Estimator) %>%
  dplyr::summarise(
    total = dplyr::n(),
    agree = sum(agreement == "Agreement"),
    disagree = sum(agreement == "Disagreement"),
    ties = sum(agreement == "Tie"),
    agreement_pct = round(100 * agree / total, 1),
    disagreement_pct = round(100 * disagree / total, 1),
    .groups = "drop"
  )

print(agreement_summary)

origin_agreement %>%
  dplyr::filter(agreement == "Tie") %>%
  dplyr::mutate(
    reason = dplyr::case_when(
      is.na(pct_flaic_better) | is.na(delta_prial) ~ "Missing value",
      pct_flaic_better == 50 & delta_prial == 0 ~ "Both boundaries",
      pct_flaic_better == 50 ~ "Exactly 50% origin wins",
      delta_prial == 0 ~ "Zero Delta PRIAL"
    )
  ) %>%
  dplyr::count(Estimator, reason)

# plot

figure_origin_prial_scatter <- ggplot(case_grid_average, aes(x = pct_flaic_better, y = delta_prial, colour = Scaling, shape = Centring)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") + geom_vline(xintercept = 50, linetype = "dotted", colour = "grey40") +
  geom_point(alpha = 0.70, size = 1.8) + facet_wrap(~Estimator, nrow = 1) + scale_colour_manual(values = scale_cols) +
  scale_shape_manual(values = c("Centred" = 16, "Uncentred" = 17)) +
  labs(title = "Origin-level win frequency and aggregate incremental PRIAL", x = "% of rolling origins where FLAIC outperforms matching IComb",
       y = expression(Delta*"PRIAL (pp)"), colour = "Scaling", shape = "Centring") +
  theme_bw(base_size = 9.5) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title = element_text(face = "plain", size = 12, hjust = 0))

ggsave(file.path("figures", "figure_origin_win_vs_incremental_prial.png"), figure_origin_prial_scatter, width = 9.5, height = 4.8, dpi = 300)

case_grid_quadrants <- case_grid_average %>% dplyr::count(Estimator, Quadrant) %>% dplyr::group_by(Estimator) %>% dplyr::mutate(pct = 100 * n/sum(n)) %>% dplyr::ungroup()
origin_metric_correlations <- case_grid_average %>% dplyr::group_by(Estimator) %>%
  dplyr::summarise(pearson = cor(pct_flaic_better, delta_prial), spearman = cor(pct_flaic_better, delta_prial, method = "spearman"), .groups = "drop")
print(case_grid_quadrants, n = Inf)
print(origin_metric_correlations, n = Inf)

discordant_example <- case_grid_average %>% dplyr::filter((delta_prial > 0 & pct_flaic_better < 50) | (delta_prial < 0 & pct_flaic_better > 50)) %>%
  dplyr::slice_max(order_by = abs(delta_prial), n = 1, with_ties = FALSE)
if (nrow(discordant_example) == 0L) discordant_example <- case_grid_average %>% dplyr::slice_max(order_by = abs(delta_prial), n = 1, with_ties = FALSE)
illustrative_origin_case <- make_origin_case(discordant_example$horizon, discordant_example$p, discordant_example$center, discordant_example$scale,
                                            discordant_example$penalty, group = "Average")
illustrative_origin_case_summary <- summarise_origin_case(illustrative_origin_case)
figure_origin_contribution_hist <- plot_origin_contribution(illustrative_origin_case, discordant_example$horizon, discordant_example$p,
                                                           discordant_example$center, discordant_example$scale, discordant_example$penalty, "Average")
print(discordant_example)
print(illustrative_origin_case_summary)
ggsave(file.path("figures", "figure_origin_contribution_histogram.png"), figure_origin_contribution_hist, width = 7.8, height = 5.0, dpi = 300)

# chosen cases
# first
x_h1_p5_ridge_cu <- make_origin_case(hor = 1, p = 5, center = TRUE, scale = FALSE, penalty = "ridge", group = "Average")
x_h1_p5_ridge_cu_summary <- summarise_origin_case(x_h1_p5_ridge_cu)
fig_h1_p5_ridge_cu <- plot_origin_contribution(x_h1_p5_ridge_cu, hor = 1, p = 5, center = TRUE, scale = FALSE, penalty = "ridge", group = "Average", show_title = FALSE)

x_h1_p5_ridge_cu_summary
x_h1_p5_ridge_cu_summary$median_relative_gain
fig_h1_p5_ridge_cu

ggsave("figures/figure_origin_contribution_h1_p5_flaic_s_centred_unscaled.png", fig_h1_p5_ridge_cu, width = 9.5, height = 5.5, dpi = 300)

# second
x_h1_p20_lasso_cs <- make_origin_case(hor = 1, p = 20, center = TRUE, scale = TRUE, penalty = "lasso", group = "Average")
x_h1_p20_lasso_cs_summary <- summarise_origin_case(x_h1_p20_lasso_cs)
fig_x_h1_p20_lasso_cs <- plot_origin_contribution(x_h1_p20_lasso_cs, hor = 1, p = 20, center = TRUE, scale = TRUE, penalty = "lasso", group = "Average", show_title = FALSE)

x_h1_p20_lasso_cs_summary
x_h1_p20_lasso_cs_summary$median_relative_gain
x_h1_p20_lasso_cs_summary$q05_relative_gain
x_h1_p20_lasso_cs_summary$q95_relative_gain
fig_x_h1_p20_lasso_cs

ggsave("figures/figure_origin_contribution_h1_p20_lasso_cs.png", fig_x_h1_p20_lasso_cs, width = 9.5, height = 5.5, dpi = 300)

### Selected forecast origin-level PRIAL contribution cases

# Panel A: agreement case — FLAIC wins more often and aggregate Delta PRIAL is positive
x_origin_A <- make_origin_case(hor = 1, p = 5, center = TRUE, scale = FALSE, penalty = "ridge", group = "Average")
summary_origin_A <- summarise_origin_case(x_origin_A)
fig_origin_A <- plot_origin_contribution(x_origin_A, hor = 1, p = 5, center = TRUE, scale = FALSE, penalty = "ridge", group = "Average", show_title = FALSE) +
  labs(title = "Panel A. FLAIC-S")

# Panel B: disagreement case — FLAIC wins more often but aggregate Delta PRIAL is negative
x_origin_B <- make_origin_case(hor = 1, p = 20, center = TRUE, scale = TRUE, penalty = "lasso", group = "Average")
summary_origin_B <- summarise_origin_case(x_origin_B)
fig_origin_B <- plot_origin_contribution(x_origin_B, hor = 1, p = 20, center = TRUE, scale = TRUE, penalty = "lasso", group = "Average", show_title = FALSE) +
  labs(title = "Panel B. FLAIC-L")

figure_origin_contribution_selected <- fig_origin_A + fig_origin_B +
  patchwork::plot_layout(ncol = 2) &
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8.5),
        axis.title = element_text(size = 8.5),
        axis.text = element_text(size = 8))

summary_origin_A
summary_origin_B
figure_origin_contribution_selected

ggsave("figures/figure_origin_contribution_selected_cases.png", figure_origin_contribution_selected, width = 13.5, height = 5.8, dpi = 300)

### Analysis 4: multiple comparisons with the best (MCB) -----------------------------------------------------------------------------------------

mcb_h <- c(1, 6, 12)
mcb_p <- as.character(analysis_p)
mcb_conf <- 0.95

series_msfe_mcb <- function(errors, hor) {
  n_valid <- nrow(errors) - (hor - 1L)
  x <- errors[seq_len(n_valid), , drop = FALSE]
  stopifnot(all(is.finite(x)))
  colMeans(x^2)
}

make_mcb_one <- function(hor, estimator, center, scale, conf_level = mcb_conf, include_base = FALSE) {
  icomb <- series_msfe_mcb(error_array[, paste0("h", hor), , get_icomb_method(estimator)], hor)
  flaic <- sapply(mcb_p, function(p) series_msfe_mcb(error_array_flaic[, paste0("h", hor), , estimator, p, as.character(center), as.character(scale)], hor))
  if (include_base) {
    base <- series_msfe_mcb(error_array[, paste0("h", hor), , "base"], hor)
    msfe <- cbind(Base = base, IComb = icomb, flaic)
    colnames(msfe) <- c("Base", "IComb", paste0("p=", mcb_p))
  } else {
    msfe <- cbind(IComb = icomb, flaic)
    colnames(msfe) <- c("IComb", paste0("p=", mcb_p))
  }
  ranks <- t(apply(msfe, 1, rank, ties.method = "average"))
  avg_rank <- colMeans(ranks)
  n <- nrow(msfe)
  k <- ncol(msfe)
  cd <- stats::qtukey(conf_level, nmeans = k, df = Inf)/sqrt(2) * sqrt(k * (k + 1)/(6 * n))
  best_rank <- min(avg_rank)
  cutoff <- best_rank + cd
  friedman_p <- stats::friedman.test(msfe)$p.value
  out <- tibble::tibble(h = hor, estimator = estimator, center = center, scale = scale, Method = names(avg_rank), avg_rank = as.numeric(avg_rank),
                        best_rank = best_rank, CD = cd, cutoff = cutoff, friedman_p = friedman_p,
                        Status = dplyr::case_when(abs(avg_rank - best_rank) < 1e-10 ~ "Best", avg_rank <= cutoff ~ "Not significantly different from best", TRUE ~ "Significantly worse than best"))
  if (include_base) {
    base_rank <- unname(avg_rank["Base"])
    out <- out %>% dplyr::mutate(base_rank = base_rank, rank_diff_vs_base = base_rank - avg_rank,
                                 Base_comparison = dplyr::case_when(Method == "Base" ~ "Base", avg_rank < base_rank & rank_diff_vs_base > CD ~ "Significantly better than base",
                                                                    avg_rank < base_rank ~ "Better rank, not significantly different", TRUE ~ "Not better than base"))
  }
  out
}

make_mcb_figure <- function(center, scale, title_spec, file_stub = NULL, show_title = FALSE) {
  df <- dplyr::bind_rows(lapply(c("ridge", "lasso"), function(est) dplyr::bind_rows(lapply(mcb_h, function(hh) make_mcb_one(hh, est, center, scale)))))
  method_levels <- c("IComb", paste0("p=", mcb_p))
  df <- df %>% dplyr::mutate(Method = factor(Method, levels = rev(method_levels)),
                             Horizon = factor(paste0("h = ", h), levels = paste0("h = ", mcb_h)),
                             Estimator = factor(ifelse(estimator == "ridge", "Panel A.\nFLAIC-S", "Panel B.\nFLAIC-L"), levels = c("Panel A.\nFLAIC-S", "Panel B.\nFLAIC-L")),
                             Status = factor(Status, levels = c("Best", "Not significantly different from best", "Significantly worse than best")))
  
  meta <- df %>% dplyr::distinct(estimator, h, Estimator, Horizon, best_rank, CD, cutoff, friedman_p)
  cd_txt <- sprintf("%.2f", unique(round(meta$CD, 10))[1])
  n_small <- sum(meta$friedman_p < 0.001)
  number_words <- c("zero", "one", "two", "three", "four", "five", "six")
  panel_txt <- if (n_small == 6) "Friedman p < 0.001 in all six panels" else paste0("Friedman p < 0.001 in ", number_words[n_small + 1], ifelse(n_small == 1, " panel", " panels"))
  
  exceptions <- meta %>% dplyr::filter(friedman_p >= 0.001) %>%
    dplyr::mutate(estimator_label = ifelse(estimator == "ridge", "Ridge", "Lasso"),
                  p_txt = ifelse(friedman_p < 0.01, sprintf("%.4f", friedman_p), sprintf("%.3f", friedman_p)),
                  txt = paste0(estimator_label, ", h = ", h, ": p = ", p_txt))
  
  friedman_txt <- if (nrow(exceptions) == 0) panel_txt else paste0(panel_txt, "; ", paste(exceptions$txt, collapse = "; "))
  stats_txt <- paste0(stringr::str_to_sentence(title_spec), " PCA: 95% Nemenyi critical distance = ", cd_txt, " | ", friedman_txt)
  
  fig <- ggplot(df, aes(x = avg_rank, y = Method)) +
    geom_rect(data = meta, aes(xmin = 1, xmax = cutoff, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey94") +
    geom_vline(data = meta, aes(xintercept = cutoff), inherit.aes = FALSE, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    geom_segment(aes(x = 1, xend = avg_rank, yend = Method), colour = "grey75", linewidth = 0.4) +
    geom_point(aes(shape = Status), size = 2.3, stroke = 0.8) +
    facet_grid(Estimator ~ Horizon, switch = "y") +
    scale_shape_manual(values = c("Best" = 17, "Not significantly different from best" = 16, "Significantly worse than best" = 1)) +
    scale_x_continuous(breaks = seq_len(length(method_levels)), limits = c(1, length(method_levels))) +
    labs(title = if (show_title) paste0("Multiple comparisons with the best for FLAIC under ", title_spec, " PCA") else NULL,
         x = "Average rank (lower = better)", y = NULL, shape = NULL, caption = stats_txt) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          strip.background = element_rect(fill = "grey90"), strip.placement = "outside",
          strip.text.x = element_text(face = "plain", size = 9.5),
          strip.text.y.left = element_text(angle = 0, face = "plain", size = 9.5),
          legend.position = "bottom", legend.text = element_text(size = 8),
          axis.title.x = element_text(size = 9), axis.text = element_text(size = 8.5),
          plot.title = element_text(face = "plain", size = 12.5, hjust = 0),
          plot.caption = element_text(size = 8.5, hjust = 0, face = "plain", margin = margin(t = 9)),
          plot.margin = margin(6, 6, 6, 6))
  
  tests <- meta %>% dplyr::mutate(Friedman_p = friedman_p, Nemenyi_CD = CD) %>%
    dplyr::select(Estimator, Horizon, Friedman_p, Nemenyi_CD)
  
  if (!is.null(file_stub)) ggsave(paste0("figures/", file_stub, ".png"), fig, width = 11.5, height = 7.7, dpi = 300)
  list(plot = fig, tests = tests, data = df, stats_text = stats_txt)
}

mcb_ct_sf <- make_mcb_figure(TRUE, FALSE, "centred, unscaled", "figure_mcb_centred_unscaled")
mcb_ct_st <- make_mcb_figure(TRUE, TRUE, "centred, scaled", "figure_mcb_centred_scaled")
mcb_cf_sf <- make_mcb_figure(FALSE, FALSE, "uncentred, unscaled", "figure_mcb_uncentred_unscaled")
mcb_cf_st <- make_mcb_figure(FALSE, TRUE, "uncentred, scaled", "figure_mcb_uncentred_scaled")

mcb_figures <- list(centred_unscaled = mcb_ct_sf, centred_scaled = mcb_ct_st, uncentred_unscaled = mcb_cf_sf, uncentred_scaled = mcb_cf_st)

mcb_ct_sf$stats_text
mcb_ct_st$stats_text
mcb_cf_sf$stats_text
mcb_cf_st$stats_text

#### Diagnostic MCB including base forecasts -----------------------------------------------------------------------------------------------------
# This diagnostic is supplementary material. A method is significantly better than Base under the Nemenyi comparison when Base rank - method rank > CD.

mcb_base_check_all <- tidyr::expand_grid(horizon = mcb_h, estimator = c("ridge", "lasso"), center = c(TRUE, FALSE), scale = c(FALSE, TRUE)) %>%
  purrr::pmap_dfr(function(horizon, estimator, center, scale) make_mcb_one(horizon, estimator, center, scale, include_base = TRUE)) %>%
  dplyr::mutate(CS = factor(get_preprocessing_label(center, scale), levels = cs_levels))

mcb_base_claim_summary <- mcb_base_check_all %>% dplyr::filter(Method != "Base") %>%
  dplyr::group_by(h, estimator, CS) %>%
  dplyr::summarise(all_nonbase_significantly_better = all(Base_comparison == "Significantly better than base"),
                   icomb_significantly_better = Base_comparison[Method == "IComb"] == "Significantly better than base",
                   all_flaic_significantly_better = all(Base_comparison[grepl("^p=", Method)] == "Significantly better than base"), .groups = "drop")
print(mcb_base_claim_summary, n = Inf)
saveRDS(mcb_base_check_all, "rds/mcb_with_base_diagnostic.rds")

plot_mcb_with_base <- function(estimator = c("ridge", "lasso"), center = TRUE, scale = FALSE, hor = 1) {
  estimator <- match.arg(estimator)
  df <- mcb_base_check_all %>% dplyr::filter(.data$estimator == .env$estimator, .data$center == .env$center, .data$scale == .env$scale, h == hor)
  method_levels <- c("Base", "IComb", paste0("p=", mcb_p))
  df <- df %>% dplyr::mutate(Method = factor(Method, levels = rev(method_levels)), Method_type = dplyr::case_when(as.character(Method) == "Base" ~ "Base", as.character(Method) == "IComb" ~ "IComb", TRUE ~ "FLAIC"))
  meta <- df %>% dplyr::slice(1)
  p_txt <- ifelse(meta$friedman_p < 0.001, "< 0.001", sprintf("= %.3f", meta$friedman_p))
  ggplot(df, aes(x = avg_rank, y = Method)) +
    geom_segment(aes(x = pmax(1, avg_rank - CD/2), xend = pmin(length(method_levels), avg_rank + CD/2), yend = Method), linewidth = 0.45) +
    geom_point(aes(shape = Method_type, colour = Base_comparison), size = 2.5, stroke = 0.8) +
    scale_shape_manual(values = c(Base = 15, IComb = 17, FLAIC = 16)) +
    scale_colour_manual(values = c("Base" = "black", "Significantly better than base" = "#2166AC", "Better rank, not significantly different" = "grey45", "Not better than base" = "#B2182B")) +
    scale_x_continuous(limits = c(1, length(method_levels)), breaks = seq_len(length(method_levels))) +
    labs(title = paste0(ifelse(estimator == "ridge", "Ridge", "Lasso"), " IComb/FLAIC specifications with base forecasts: h = ", hor),
         subtitle = paste0("C=", ifelse(center, "T", "F"), ", S=", ifelse(scale, "T", "F"), " | Friedman test p-value ", p_txt, " | Critical distance = ", sprintf("%.2f", meta$CD)),
         x = "Mean rank (lower is better)", y = NULL, shape = NULL, colour = NULL,
         caption = "Horizontal intervals have total width equal to the Nemenyi critical distance; the companion table gives the direct pairwise rank-distance check against Base.") +
    theme_bw(base_size = 9.5) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), legend.position = "bottom", plot.title = element_text(face = "plain"))
}

mcb_with_base_example <- plot_mcb_with_base("ridge", TRUE, FALSE, 1)
mcb_with_base_example
mcb_with_base_example_2 <- plot_mcb_with_base("ridge", TRUE, FALSE, 6)
mcb_with_base_example_2
mcb_with_base_example_3 <- plot_mcb_with_base("ridge", FALSE, FALSE, 6)
mcb_with_base_example_3
mcb_with_base_example_4 <- plot_mcb_with_base("ridge", FALSE, TRUE, 6)
mcb_with_base_example_4
mcb_with_base_example_5 <- plot_mcb_with_base("ridge", TRUE, TRUE, 6)
mcb_with_base_example_5

# diagnose things


mcb_base_check_all %>%
  dplyr::filter(Method != "Base") %>%
  dplyr::group_by(h, estimator) %>%
  dplyr::summarise(
    comparisons = dplyr::n(),
    significantly_better = sum(
      Base_comparison == "Significantly better than base"
    ),
    not_significantly_better = sum(
      Base_comparison != "Significantly better than base"
    ),
    .groups = "drop"
  ) %>%
  print(n = Inf)

# Examine borderline cases without rounding
mcb_base_check_all %>%
  dplyr::filter(
    h == 6,
    estimator == "ridge",
    Method == "p=77",
    !scale
  ) %>%
  dplyr::mutate(
    excess_over_CD = rank_diff_vs_base - CD
  ) %>%
  dplyr::select(
    center, scale, Method,
    rank_diff_vs_base, CD, excess_over_CD,
    Base_comparison
  ) %>%
  print(n = Inf, digits = 12)



### Analysis 5: line-graph alternatives to the heatmaps ------------------------------------------------------------------------------------------

make_sensitivity_df <- function(groups, centers_use, scales_use, horizons_use = analysis_h) {
  tidyr::expand_grid(horizon = paste0("h", horizons_use), group = groups, penalty = c("ridge", "lasso"), npc = as.character(analysis_p), center = centers_use, scale = scales_use) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(prial_icomb = as.numeric(prial_array[horizon, group, get_icomb_method(penalty)]),
                  prial_flaic = as.numeric(prial_array_flaic[horizon, group, penalty, npc, as.character(center), as.character(scale)]),
                  gain_pp = prial_icomb - prial_flaic) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p = factor(npc, levels = as.character(analysis_p)), Estimator = factor(ifelse(penalty == "ridge", "Ridge", "Lasso"), levels = c("Ridge", "Lasso")),
                  Scaling = factor(ifelse(scale, "Scaled", "Unscaled"), levels = c("Unscaled", "Scaled")),
                  Centring = factor(ifelse(center, "Centred", "Uncentred"), levels = c("Centred", "Uncentred")),
                  CS = factor(get_cs_label(center, scale), levels = cs_levels))
}

sens_all <- make_sensitivity_df("Average", c(FALSE, TRUE), c(FALSE, TRUE), analysis_h)
common_ylim <- function(x) {
  r <- range(c(0, x), na.rm = TRUE)
  pad <- max(diff(r) * 0.06, 0.05)
  c(r[1] - pad, r[2] + pad)
}
ylim_all <- common_ylim(sens_all$gain_pp)

make_sensitivity_horizon_panel <- function(hor, show_y = FALSE, show_x = FALSE) {
  p <- ggplot(dplyr::filter(sens_all, horizon == paste0("h", hor)), aes(x = p, y = gain_pp, colour = Scaling, linetype = Centring, group = interaction(Scaling, Centring))) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45", linewidth = 0.3) + geom_line(linewidth = 0.55) + geom_point(size = 1.25) +
    facet_wrap(~Estimator, nrow = 1) + scale_colour_manual(values = scale_cols, name = "Scaling") + scale_linetype_manual(values = center_ltypes, name = "Centring") +
    coord_cartesian(ylim = ylim_all) + labs(title = paste0("h = ", hor), x = if (show_x) "Number of principal components (p)" else NULL,
                                           y = if (show_y) "Incremental PRIAL gain (pp)\n(positive = improvement)" else NULL) +
    theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey90"), strip.text = element_text(face = "bold", size = 7.4),
                                    plot.title = element_text(face = "plain", size = 8.5, hjust = 0.5), legend.position = "bottom", axis.title = element_text(size = 7.5), axis.text = element_text(size = 7))
  if (!show_y) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p
}

sens_all_plots <- lapply(seq_along(analysis_h), function(i) {
  row_i <- ceiling(i/4)
  col_i <- ((i - 1) %% 4) + 1
  make_sensitivity_horizon_panel(analysis_h[i], show_y = col_i == 1, show_x = row_i == 3)
})

figure_sensitivity_all_h <- patchwork::wrap_plots(sens_all_plots, ncol = 4, guides = "collect") +
  patchwork::plot_annotation(title = "Sensitivity of FLAIC performance to PCA preprocessing specification across forecast horizons",
                             caption = stringr::str_wrap("Incremental PRIAL gain is PRIAL(IComb) minus PRIAL(FLAIC), so positive values indicate improvement. Blue/orange distinguish unscaled/scaled PCA and solid/dashed lines distinguish centred/uncentred PCA. Results are averaged across all 112 hierarchy series.", width = 190),
                             theme = theme(plot.title = element_text(face = "plain", size = 12.5, hjust = 0), plot.caption = element_text(size = 7.5, hjust = 0, lineheight = 1.05, margin = margin(t = 6)))) &
  theme(legend.position = "bottom")
ggsave(file.path("figures", "figure_flaic_pca_sensitivity_all_h.png"), figure_sensitivity_all_h, width = 15.5, height = 8.5, dpi = 300)

sens_hier <- make_sensitivity_df(analysis_levels, c(FALSE, TRUE), c(FALSE, TRUE), c(1, 6, 12)) %>%
  dplyr::mutate(Specification = factor(CS, levels = cs_levels), Horizon = factor(horizon, levels = paste0("h", c(1, 6, 12)), labels = paste0("h = ", c(1, 6, 12))),
                Level = factor(group, levels = analysis_levels))
ylim_hier <- common_ylim(sens_hier$gain_pp)

figure_hierarchy_all_specs <- ggplot(sens_hier, aes(x = p, y = gain_pp, colour = Level, group = Level)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45", linewidth = 0.3) + geom_line(linewidth = 0.50) + geom_point(size = 1.1) +
  facet_grid(rows = vars(Specification), cols = vars(Horizon, Estimator), switch = "y") + scale_colour_manual(values = level_cols, name = "Hierarchy level") +
  coord_cartesian(ylim = ylim_hier) +
  labs(title = "FLAIC performance by hierarchy level under alternative PCA preprocessing specifications", x = "Number of principal components (p)",
       y = "Incremental PRIAL gain (FLAIC vs matching IComb; pp)\n(positive = improvement)",
       caption = stringr::str_wrap("Rows show all four PCA preprocessing specifications, with centred-unscaled PCA shown first. Lines report Australia, States, Zones and Regions; Average is omitted. Ridge FLAIC is compared with IComb-S and Lasso FLAIC with IComb-L.", width = 190)) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey90"), strip.placement = "outside", strip.text.x = element_text(face = "plain", size = 7.5),
        strip.text.y.left = element_text(angle = 0, face = "plain", size = 8), legend.position = "bottom", axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        plot.title = element_text(face = "plain", size = 12.5, hjust = 0), plot.caption = element_text(size = 7.5, hjust = 0, lineheight = 1.05, margin = margin(t = 6)))
ggsave(file.path("figures", "figure_flaic_hierarchy_all_pca_specs.png"), figure_hierarchy_all_specs, width = 15.5, height = 11.2, dpi = 300)

### Analysis 6: lasso predictor selection from saved final-run coefficients ----------------------------------------------------------------------

coef_tol <- 1e-10

get_flaic_coef_spec <- function(x, p, center, scale, estimator = "lasso") {
  idx <- which(vapply(x, function(z) as.integer(z$p) == as.integer(p) && as.logical(z$center) == center && as.logical(z$scale) == scale, logical(1)))
  if (length(idx) != 1L) stop("Could not identify exactly one FLAIC coefficient specification for p=", p, ", center=", center, ", scale=", scale)
  x[[idx]][[estimator]]
}

icomb_drop_origin <- purrr::map_dfr(origins, function(T) {
  B <- coef_array_icomb[[as.character(T), "icombl_xy1"]]
  slopes <- B[-1, , drop = FALSE]
  stopifnot(nrow(slopes) == ny)
  dropped <- rowSums(abs(slopes) > coef_tol) == 0
  tibble::tibble(origin = T, icomb_orig_dropped = sum(dropped))
})

flaic_selection_origin <- vector("list", length(origins) * length(analysis_p) * 4)
pc_retention_detail <- vector("list", length(origins) * length(analysis_p) * 4)
k1 <- 0L
k2 <- 0L

for (T in origins) {
  coef_T <- readRDS(file.path("coefficients", paste0("flaic_coefficients_origin_", T, ".rds")))
  icomb_dropped <- icomb_drop_origin$icomb_orig_dropped[icomb_drop_origin$origin == T]
  for (C in c(TRUE, FALSE)) {
    for (S in c(FALSE, TRUE)) {
      for (p in analysis_p) {
        B <- get_flaic_coef_spec(coef_T, p, C, S, "lasso")
        slopes <- B[-1, , drop = FALSE]
        stopifnot(nrow(slopes) == ny + p)
        hierarchy_dropped <- rowSums(abs(slopes[seq_len(ny), , drop = FALSE]) > coef_tol) == 0
        pc_dropped <- rowSums(abs(slopes[ny + seq_len(p), , drop = FALSE]) > coef_tol) == 0
        k1 <- k1 + 1L
        flaic_selection_origin[[k1]] <- tibble::tibble(origin = T, p = p, center = C, scale = S, CS = get_cs_label(C, S),
                                                       icomb_orig_dropped = icomb_dropped, flaic_orig_dropped = sum(hierarchy_dropped),
                                                       diff_orig_dropped = sum(hierarchy_dropped) - icomb_dropped, pc_retained = sum(!pc_dropped), pc_retained_pct = 100 * mean(!pc_dropped))
        k2 <- k2 + 1L
        pc_retention_detail[[k2]] <- tibble::tibble(origin = T, p = p, center = C, scale = S, CS = get_cs_label(C, S), pc = seq_len(p), retained = !pc_dropped)
      }
    }
  }
}

flaic_selection_origin <- dplyr::bind_rows(flaic_selection_origin[seq_len(k1)])
pc_retention_detail <- dplyr::bind_rows(pc_retention_detail[seq_len(k2)])

#### Lasso selection summaries and table ---------------------------------------------------------------------------------------------------------

pca_order <- c("Centred, unscaled", "Centred, scaled", "Uncentred, unscaled", "Uncentred, scaled")

lasso_selection_summary <- flaic_selection_origin %>%
  dplyr::mutate(PCA = factor(get_preprocessing_label(center, scale), levels = pca_order),
                delta_hierarchy_retained = icomb_orig_dropped - flaic_orig_dropped) %>%
  dplyr::group_by(PCA, p) %>%
  dplyr::summarise(delta_hier_mean = mean(delta_hierarchy_retained), delta_hier_median = median(delta_hierarchy_retained),
                   delta_hier_min = min(delta_hierarchy_retained), delta_hier_max = max(delta_hierarchy_retained),
                   pc_ret_mean = mean(pc_retained), pc_ret_median = median(pc_retained), pc_ret_min = min(pc_retained), pc_ret_max = max(pc_retained),
                   pc_pct_mean = mean(pc_retained_pct), pc_pct_median = median(pc_retained_pct), .groups = "drop") %>%
  dplyr::arrange(PCA, p)

icomb_retained_median <- ny - median(icomb_drop_origin$icomb_orig_dropped)
fmt_compact <- function(x) ifelse(abs(x - round(x)) < 1e-10, sprintf("%.0f", x), sprintf("%.1f", x))
fmt_pc_main <- function(n, pct) paste0(fmt_compact(n), " (", sprintf("%.1f", pct), "\\%)")

lasso_selection_table <- lasso_selection_summary %>%
  dplyr::mutate(`Delta hierarchy retained` = fmt_compact(delta_hier_median),
                `PC forecasts retained` = fmt_pc_main(pc_ret_median, pc_pct_median)) %>%
  dplyr::select(p, PCA, `Delta hierarchy retained`, `PC forecasts retained`) %>%
  tidyr::pivot_wider(names_from = PCA, values_from = c(`Delta hierarchy retained`, `PC forecasts retained`),
                     names_glue = "{PCA}__{.value}") %>%
  dplyr::arrange(p)

saveRDS(lasso_selection_table, "rds/table_R2_lasso_predictor_selection.rds")
saveRDS(list(icomb_retained_median = icomb_retained_median, hierarchy_total = ny), "rds/table_R2_lasso_predictor_selection_meta.rds")
saveRDS(list(icomb_drop_origin = icomb_drop_origin, flaic_selection_origin = flaic_selection_origin,
             pc_retention_detail = pc_retention_detail, summary = lasso_selection_summary), "rds/lasso_selection_diagnostics.rds")

#### PC forecast retention heatmap: p = 20 

make_pc_retention_heatmap <- function(p_show = 20, show_title = TRUE, label_size = 3.1, label_bold = TRUE) {
  plot_df <- pc_retention_detail %>%
    dplyr::filter(p == p_show) %>%
    dplyr::group_by(center, scale, pc) %>%
    dplyr::summarise(retained_pct = 100 * mean(retained), .groups = "drop") %>%
    dplyr::mutate(PCA = factor(get_preprocessing_label(center, scale), levels = pca_order),
                  PC = factor(paste0("PC", pc), levels = paste0("PC", seq_len(p_show))),
                  label = sprintf("%.0f", retained_pct),
                  text_col = ifelse(retained_pct >= 60, "white", "black"))
  
  ggplot(plot_df, aes(x = PC, y = PCA, fill = retained_pct)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    geom_text(aes(label = label, colour = text_col), size = label_size, fontface = ifelse(label_bold, "bold", "plain")) +
    scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 100), breaks = seq(0, 100, 20),
                        labels = function(x) paste0(x, "%"), name = "Retained\n(% forecast origins)") +
    scale_colour_identity() +
    labs(title = if (show_title) paste0("Lasso retention of principal-component forecasts at p = ", p_show) else NULL,
         x = "Principal-component forecast", y = "PCA specification") +
    theme_bw(base_size = 9.5) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 8.5), axis.title = element_text(size = 9),
          legend.position = "right", plot.title = element_text(face = "plain", size = 12, hjust = 0))
}

pc_retention_p20_preview <- make_pc_retention_heatmap(20, show_title = TRUE)
pc_retention_p20 <- make_pc_retention_heatmap(20, show_title = FALSE)

pc_retention_p20_preview
ggsave("figures/figure_pc_retention_p20.png", pc_retention_p20, width = 11.5, height = 4.8, dpi = 300)

### Analysis 7: FLAP benchmark tables for Quarto -------------------------------------------------------------------------------------------------

# Table R4: FLAP PRIAL for Average performance.
# Table R5: Panel A = MinT_shr - FLAP; Panel B = FLAP - FLAIC-S; Panel C = FLAP - FLAIC-L.
flap_table_p <- dimnames(prial_array_flap)$p
flap_table_cs <- tibble::tibble(center = c(FALSE, FALSE, TRUE, TRUE), scale = c(FALSE, TRUE, FALSE, TRUE))
flap_table_h <- c("h1", "h6", "h12", "1-6", "1-12")
flap_tf <- function(x) ifelse(x, "T", "F")
flap_hkey <- function(x) ifelse(x == "1-6", "h16", ifelse(x == "1-12", "h112", x))

make_flap_table_rows <- function() {
  purrr::map_dfr(seq_len(nrow(flap_table_cs)), function(i) {
    C <- flap_table_cs$center[i]
    S <- flap_table_cs$scale[i]
    tibble::tibble(CS = paste0("C=", flap_tf(C), ", S=", flap_tf(S)), p = flap_table_p, center = C, scale = S)
  })
}

make_flap_prial_table <- function(prial_flap, group = "Average", horizons = flap_table_h) {
  specs <- make_flap_table_rows()
  out <- specs
  for (hor in horizons) {
    key <- flap_hkey(hor)
    out[[key]] <- vapply(seq_len(nrow(specs)), function(i) as.numeric(prial_flap[hor, group, specs$p[i], as.character(specs$center[i]), as.character(specs$scale[i])]), numeric(1))
  }
  out
}

make_flap_comparison_table <- function(prial_flap, prial_flaic, prial_comparison, group = "Average", horizons = flap_table_h) {
  specs <- make_flap_table_rows()
  out <- specs
  for (hor in horizons) {
    key <- flap_hkey(hor)
    flap <- vapply(seq_len(nrow(specs)), function(i) as.numeric(prial_flap[hor, group, specs$p[i], as.character(specs$center[i]), as.character(specs$scale[i])]), numeric(1))
    mint <- as.numeric(prial_comparison[hor, group, "mint_shr"])
    flaic_s <- vapply(seq_len(nrow(specs)), function(i) as.numeric(prial_flaic[hor, group, "ridge", specs$p[i], as.character(specs$center[i]), as.character(specs$scale[i])]), numeric(1))
    flaic_l <- vapply(seq_len(nrow(specs)), function(i) as.numeric(prial_flaic[hor, group, "lasso", specs$p[i], as.character(specs$center[i]), as.character(specs$scale[i])]), numeric(1))
    out[[paste0("A__", key)]] <- mint - flap
    out[[paste0("B__", key)]] <- flap - flaic_s
    out[[paste0("C__", key)]] <- flap - flaic_l
  }
  out
}

table_R4_flap_prial <- make_flap_prial_table(prial_array_flap, "Average", flap_table_h)
table_R5_flap_comparison <- make_flap_comparison_table(prial_array_flap, prial_array_flaic, prial_array, "Average", flap_table_h)
saveRDS(table_R4_flap_prial, "rds/table_R4_flap_prial_average.rds")
saveRDS(table_R5_flap_comparison, "rds/table_R5_flap_method_comparison.rds")
message("FLAP PRIAL and FLAP comparison tables saved for Quarto.")

### Analysis 8: improved explore_flaic wrapper ---------------------------------------------------------------------------------------------------

plot_prial_by_p_levels <- function(center = TRUE, scale = FALSE, horizon = "1-12") {
  groups <- analysis_levels
  plot_df <- tidyr::expand_grid(group = groups, penalty = c("ridge", "lasso"), p = analysis_p) %>% dplyr::rowwise() %>%
    dplyr::mutate(gain_pp = as.numeric(prial_array[horizon, group, get_icomb_method(penalty)] - prial_array_flaic[horizon, group, penalty, as.character(p), as.character(center), as.character(scale)])) %>%
    dplyr::ungroup() %>% dplyr::mutate(p = factor(p, levels = analysis_p), Estimator = factor(ifelse(penalty == "ridge", "FLAIC-S", "FLAIC-L"), levels = c("FLAIC-S", "FLAIC-L")),
                                                 Level = factor(group, levels = analysis_levels))
  g <- ggplot(plot_df, aes(p, gain_pp, colour = Level, group = Level)) + geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45") +
    geom_line(linewidth = 0.65) + geom_point(size = 1.7) + facet_wrap(~Estimator, nrow = 1) + scale_colour_manual(values = level_cols) +
    labs(title = paste0("FLAIC performance by hierarchy level: C=", ifelse(center, "T", "F"), ", S=", ifelse(scale, "T", "F")), subtitle = paste0("Forecast horizon: ", horizon),
         x = "Number of principal components (p)", y = "Incremental PRIAL gain (pp)\n(positive = improvement)", colour = "Hierarchy level") +
    theme_bw(base_size = 9.5) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  list(plot = g, data = plot_df)
}

plot_prial_by_p_preprocessing <- function(group = "Average", horizon = "1-12") {
  plot_df <- tidyr::expand_grid(penalty = c("ridge", "lasso"), p = analysis_p, center = c(FALSE, TRUE), scale = c(FALSE, TRUE)) %>% dplyr::rowwise() %>%
    dplyr::mutate(gain_pp = as.numeric(prial_array[horizon, group, get_icomb_method(penalty)] - prial_array_flaic[horizon, group, penalty, as.character(p), as.character(center), as.character(scale)])) %>%
    dplyr::ungroup() %>% dplyr::mutate(p = factor(p, levels = analysis_p), Estimator = factor(ifelse(penalty == "ridge", "FLAIC-S", "FLAIC-L"), levels = c("FLAIC-S", "FLAIC-L")),
                                                 Scaling = factor(ifelse(scale, "Scaled", "Unscaled"), levels = c("Unscaled", "Scaled")),
                                                 Centring = factor(ifelse(center, "Centred", "Uncentred"), levels = c("Centred", "Uncentred")))
  g <- ggplot(plot_df, aes(p, gain_pp, colour = Scaling, linetype = Centring, group = interaction(Scaling, Centring))) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45") + geom_line(linewidth = 0.65) + geom_point(size = 1.7) + facet_wrap(~Estimator, nrow = 1) +
    scale_colour_manual(values = scale_cols) + scale_linetype_manual(values = center_ltypes) +
    labs(title = paste0("FLAIC sensitivity to PCA preprocessing: ", group), subtitle = paste0("Forecast horizon: ", horizon), x = "Number of principal components (p)",
         y = "Incremental PRIAL gain (pp)\n(positive = improvement)", colour = "Scaling", linetype = "Centring") +
    theme_bw(base_size = 9.5) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  list(plot = g, data = plot_df)
}

explore_flaic <- function(p_retention = 20, center = TRUE, scale = FALSE, horizon = "1-12", group = "Average") {
  if (!as.character(p_retention) %in% as.character(analysis_p)) stop("p_retention must be one of: ", paste(analysis_p, collapse = ", "))
  if (!horizon %in% dimnames(prial_array)$horizon) stop("horizon must be one of: ", paste(dimnames(prial_array)$horizon, collapse = ", "))
  if (!group %in% names(group_idx)) stop("group must be one of: ", paste(names(group_idx), collapse = ", "))
  selected_summary <- delta_prial_heat_all %>% dplyr::filter(Group == group, p == p_retention, center == .env$center, scale == .env$scale) %>%
    dplyr::left_join(origin_heat_all %>% dplyr::filter(Group == group, p == p_retention, center == .env$center, scale == .env$scale) %>%
                       dplyr::select(horizon, penalty, pct_flaic_better), by = c("horizon", "penalty"))
  levels_result <- plot_prial_by_p_levels(center, scale, horizon)
  preprocessing_result <- plot_prial_by_p_preprocessing(group, horizon)
  list(settings = list(p_retention = p_retention, center = center, scale = scale, horizon = horizon, group = group),
       plots = list(retention = make_pc_retention_heatmap(p_retention, show_title = FALSE), prial_levels = levels_result$plot, prial_preprocessing = preprocessing_result$plot,
                    incremental_prial_heatmap = prial_heatmaps[[group]], origin_win_heatmap = origin_heatmaps[[group]], origin_vs_prial_scatter = figure_origin_prial_scatter),
       data = list(retention = pc_retention_detail %>% dplyr::filter(p == p_retention), prial_levels = levels_result$data, prial_preprocessing = preprocessing_result$data,
                   incremental_prial_heatmap = delta_prial_heat_all %>% dplyr::filter(Group == group), origin_win_heatmap = origin_heat_all %>% dplyr::filter(Group == group),
                   selected_summary = selected_summary))
}

# Example:
x <- explore_flaic(p_retention = 5, center = TRUE, scale = FALSE, horizon = "1-12", group = "Average")
# x$plots$incremental_prial_heatmap
# x$plots$origin_win_heatmap
x$plots$prial_preprocessing
x$plots$prial_levels
x$plots$retention
x$data$selected_summary

### Save workbench objects -----------------------------------------------------------------------------------------------------------------------

saveRDS(
  list(
    origin_heat_all = origin_heat_all,
    delta_prial_heat_all = delta_prial_heat_all,
    origin_heatmaps = origin_heatmaps,
    prial_heatmaps = prial_heatmaps,
    hierarchy_origin_heatmaps = hierarchy_origin_heatmaps_,
    hierarchy_prial_heatmaps = hierarchy_prial_heatmaps_,
    case_grid_average = case_grid_average,
    case_grid_quadrants = case_grid_quadrants,
    origin_metric_correlations = origin_metric_correlations,
    mcb_tests = lapply(mcb_figures, function(x) x$tests),
    mcb_base_claim_summary = mcb_base_claim_summary,
    lasso_selection_summary = lasso_selection_summary,
    lasso_selection_table = lasso_selection_table,
    table_R4_flap_prial = table_R4_flap_prial,
    table_R5_flap_comparison = table_R5_flap_comparison
  ),
  "rds/flaic_analysis_workbench_objects.rds"
)

message("07b_analysis_workbench complete.")






