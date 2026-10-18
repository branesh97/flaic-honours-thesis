# Descriptive analysis for FLAIC thesis
# Produces:
#   figures/figure_d1_descriptive_tourism.png
#   figures/figure_d2_selected_pcs.png
#   rds/table_D1_pc_descriptives.rds

library(tidyverse)
library(patchwork)
library(zoo)

dir.create("figures", showWarnings = FALSE)
dir.create("rds", showWarnings = FALSE)

# Load data/settings and final-run PCA objects only.
data_setup <- readRDS("rds/flaic_01_data_setup.rds")
list2env(data_setup, envir = .GlobalEnv)
rm(data_setup)
forecast_results <- readRDS("rds/flaic_03_forecasts_complete.rds")
list2env(forecast_results, envir = .GlobalEnv)
rm(forecast_results)

required_objects <- c("ally", "origins", "centers", "scales", "group_idx",
                      "pc_actual_array", "pc_eigval_array", "pc_eigvec_array")
stopifnot(all(vapply(required_objects, exists, logical(1), inherits = TRUE)))

### Settings -------------------------------------------------------------------

initial_window <- window_size
data_start <- c(1998L, 1L)
descriptive_origin <- as.integer(origins[1])
selected_pcs <- 1:5
n_side <- 2L

state_label_lookup <- c(
  A = "NSW", B = "VIC", C = "QLD", D = "SA",
  E = "WA", F = "TAS", G = "NT"
)
state_full_names <- c(
  NSW = "New South Wales", VIC = "Victoria", QLD = "Queensland",
  SA = "South Australia", WA = "Western Australia",
  TAS = "Tasmania", NT = "Northern Territory"
)

get_array_item <- function(arr, T, center, scale) {
  arr[[as.character(T), as.character(center), as.character(scale)]]
}

zscore <- function(x) {
  x <- as.numeric(x)
  sx <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
  as.numeric(scale(x, center = TRUE, scale = TRUE))
}

window_indices <- function(T = descriptive_origin, window_size = initial_window) {
  seq.int(T - window_size + 1L, T)
}

make_window_ts <- function(T = descriptive_origin, window_size = initial_window) {
  idx <- window_indices(T, window_size)
  first_index <- idx[1] - 1L
  start_year <- data_start[1] + first_index %/% seasonal_period
  start_month <- data_start[2] + first_index %% seasonal_period
  ts(ally[idx, , drop = FALSE], start = c(start_year, start_month), frequency = seasonal_period)
}

### Figure D1: tourism series and seasonal profiles ----------------------------

make_figure_d1 <- function(T = descriptive_origin, window_size = initial_window) {
  idx <- window_indices(T, window_size)
  y_window <- make_window_ts(T, window_size)

  total_code <- if ("Total" %in% colnames(ally)) {
    "Total"
  } else if ("Australia" %in% colnames(ally)) {
    "Australia"
  } else {
    stop("Need a 'Total' or 'Australia' column in ally.")
  }

  state_codes <- names(state_label_lookup)
  stopifnot(all(state_codes %in% colnames(ally)))

  date_vec <- as.Date(zoo::as.yearmon(time(y_window)))
  panel_order <- c("Australia", unname(state_label_lookup))
  state_order <- unname(state_label_lookup)
  label_map <- c(setNames("Australia", total_code), state_label_lookup)

  series_df <- as.data.frame(y_window[, c(total_code, state_codes), drop = FALSE]) %>%
    dplyr::mutate(Date = date_vec) %>%
    tidyr::pivot_longer(-Date, names_to = "series_code", values_to = "value") %>%
    dplyr::mutate(series = factor(unname(label_map[series_code]), levels = panel_order))

  # Standardise each state/territory series over the initial window using
  # (x - sample mean) / sample standard deviation, then average by calendar month.
  state_z <- scale(as.matrix(y_window[, state_codes, drop = FALSE]), center = TRUE, scale = TRUE)
  season_df <- as.data.frame(state_z) %>%
    dplyr::mutate(month = factor(cycle(y_window), levels = 1:12, labels = month.abb)) %>%
    tidyr::pivot_longer(-month, names_to = "state_code", values_to = "value") %>%
    dplyr::group_by(month, state_code) %>%
    dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(`State or territory` = factor(unname(state_label_lookup[state_code]), levels = state_order))

  state_cols <- c(
    NSW = "#D55E00", VIC = "#7570B3", QLD = "#009E73", SA = "#A6761D",
    WA = "#CC79A7", TAS = "#0072B2", NT = "#A9A9A9"
  )

  fig_theme <- theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold", size = 8.5),
      plot.title = element_text(face = "bold", size = 10.5, margin = margin(b = 4)),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )

  p1 <- ggplot(series_df, aes(Date, value)) +
    geom_line(linewidth = 0.35) +
    facet_wrap(~series, ncol = 2, scales = "free_y") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(title = paste0("Panel A. National and state/territory tourism series, ",
                        format(min(date_vec), "%b %Y"), " – ", format(max(date_vec), "%b %Y")),
         x = NULL, y = NULL) +
    fig_theme

  p2 <- ggplot(season_df, aes(month, value, group = `State or territory`, colour = `State or territory`)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.5) +
    scale_colour_manual(values = state_cols) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(title = "Panel B. Average monthly seasonal profiles", x = NULL,
         y = "Average standardised tourism", colour = "State or territory") +
    fig_theme

  p <- p1 / p2 + patchwork::plot_layout(heights = c(2.25, 1.35))

  list(plot = p, series_df = series_df, season_df = season_df,
       window_start = min(date_vec), window_end = max(date_vec))
}

# Main-text initial window plus appendix mid/end rolling windows.
figure_d1 <- make_figure_d1(T = 120L)
figure_d1_mid <- make_figure_d1(T = 192L)
figure_d1_end <- make_figure_d1(T = 263L)

figure_d1$plot

ggsave("figures/figure_d1_descriptive_tourism.png", figure_d1$plot,
       width = 8.3, height = 9.2, dpi = 300)
ggsave("figures/figure_d1_descriptive_tourism_T192.png", figure_d1_mid$plot,
       width = 8.3, height = 9.2, dpi = 300)
ggsave("figures/figure_d1_descriptive_tourism_T263.png", figure_d1_end$plot,
       width = 8.3, height = 9.2, dpi = 300)

### Figure D2: selected PC scores and illustrative observed patterns ------------

make_figure_d2 <- function(T = descriptive_origin, center = TRUE, scale = FALSE,
                           pcs = selected_pcs, n_side = 2L,
                           window_size = initial_window) {
  idx <- window_indices(T, window_size)
  y_window <- make_window_ts(T, window_size)
  time_num <- as.numeric(time(y_window))

  scores <- get_array_item(pc_actual_array, T, center, scale)
  loadings <- get_array_item(pc_eigvec_array, T, center, scale)

  state_idx <- group_idx$States
  total_idx <- group_idx$Australia[1]
  state_codes <- colnames(ally)[state_idx]

  state_z <- apply(as.matrix(y_window[, state_idx, drop = FALSE]), 2, zscore)
  colnames(state_z) <- state_codes

  results <- purrr::map(pcs, function(j) {
    pc_raw <- as.numeric(scores[, j])

    if (j == 1L) {
      observed_raw <- as.numeric(y_window[, total_idx])
      if (cor(pc_raw, observed_raw, use = "complete.obs") < 0) pc_raw <- -pc_raw
      panel_label <- "PC1: common national movement"
      definition <- tibble::tibble(
        PC = "PC1", positive_group = "Australia", negative_group = NA_character_,
        interpretation = "National/common movement"
      )
    } else {
      state_loading <- as.numeric(loadings[state_idx, j])
      names(state_loading) <- state_codes

      # PCA signs are arbitrary. Orient the component so the state with the
      # largest absolute loading lies on the positive side before selecting groups.
      largest_state <- names(which.max(abs(state_loading)))
      if (state_loading[largest_state] < 0) {
        state_loading <- -state_loading
        pc_raw <- -pc_raw
      }

      positive_candidates <- sort(state_loading[state_loading > 0], decreasing = TRUE)
      negative_candidates <- sort(state_loading[state_loading < 0], decreasing = FALSE)

      if (length(positive_candidates) < n_side || length(negative_candidates) < n_side) {
        ordered <- sort(state_loading, decreasing = TRUE)
        pos_states <- names(head(ordered, n_side))
        neg_states <- names(tail(ordered, n_side))
      } else {
        pos_states <- names(head(positive_candidates, n_side))
        neg_states <- names(head(negative_candidates, n_side))
      }

      # Each state series has already been standardised. The observed pattern is
      # the simple mean of the selected positive-side states minus the simple
      # mean of the selected negative-side states.
      observed_raw <- rowMeans(state_z[, pos_states, drop = FALSE]) -
        rowMeans(state_z[, neg_states, drop = FALSE])

      pos_labels <- unname(ifelse(pos_states %in% names(state_label_lookup),
                                  state_label_lookup[pos_states], pos_states))
      neg_labels <- unname(ifelse(neg_states %in% names(state_label_lookup),
                                  state_label_lookup[neg_states], neg_states))

      panel_label <- paste0("PC", j, ": ", paste(pos_labels, collapse = "+"),
                            " vs ", paste(neg_labels, collapse = "+"))
      definition <- tibble::tibble(
        PC = paste0("PC", j),
        positive_group = paste(pos_labels, collapse = " + "),
        negative_group = paste(neg_labels, collapse = " + "),
        interpretation = "Illustrative cross-state contrast"
      )
    }

    # Standardise the plotted PC score and observed pattern separately to mean
    # zero and unit sample standard deviation for visual comparison only.
    pc_plot <- zscore(pc_raw)
    observed_plot <- zscore(observed_raw)

    # Orient the displayed observed pattern to the PC score. This has no
    # substantive effect because the sign of a principal component is arbitrary.
    if (cor(pc_plot, observed_plot, use = "complete.obs") < 0) observed_plot <- -observed_plot

    plot_df <- tibble::tibble(
      time = rep(time_num, 2),
      series = factor(rep(c("PC score", "Observed tourism pattern"), each = length(time_num)),
                      levels = c("PC score", "Observed tourism pattern")),
      value = c(pc_plot, observed_plot),
      panel = panel_label
    )

    list(data = plot_df, definition = definition)
  })

  plot_df <- dplyr::bind_rows(purrr::map(results, "data"))
  definitions <- dplyr::bind_rows(purrr::map(results, "definition"))
  plot_df$panel <- factor(plot_df$panel, levels = unique(plot_df$panel))

  p <- ggplot(plot_df, aes(time, value, colour = series, linetype = series)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~panel, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = c("PC score" = "#FF6258", "Observed tourism pattern" = "#20BFC6")) +
    scale_linetype_manual(values = c("PC score" = "solid", "Observed tourism pattern" = "dashed")) +
    scale_x_continuous(breaks = seq(floor(min(time_num)), ceiling(max(time_num)), by = 2)) +
    labs(x = NULL, y = "Standardised value", colour = NULL, linetype = NULL) +
    theme_bw(base_size = 8.5) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "plain", size = 12, hjust = 0, margin = margin(b = 4)),
      strip.text = element_text(face = "bold", size = 8.2),
      axis.text = element_text(size = 7.2),
      axis.title = element_text(size = 9),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      panel.spacing = unit(0.25, "lines"),
      plot.margin = margin(4, 6, 4, 6)
    )

  list(plot = p, plot_data = plot_df, definitions = definitions)
}

# Main-text initial window plus appendix mid/end rolling windows.
figure_d2 <- make_figure_d2(T = 120L, center = TRUE, scale = FALSE, pcs = 1:5, n_side = 2L)
figure_d2_mid <- make_figure_d2(T = 192L, center = TRUE, scale = FALSE, pcs = 1:5, n_side = 2L)
figure_d2_end <- make_figure_d2(T = 263L, center = TRUE, scale = FALSE, pcs = 1:5, n_side = 2L)

figure_d2$plot
figure_d2$definitions
figure_d2_mid$definitions
figure_d2_end$definitions

ggsave("figures/figure_d2_selected_pcs.png", figure_d2$plot,
       width = 8.3, height = 8.8, dpi = 300)

ggsave("figures/figure_d2_selected_pcs_T192.png", figure_d2_mid$plot,
       width = 8.3, height = 8.8, dpi = 300)

ggsave("figures/figure_d2_selected_pcs_T263.png", figure_d2_end$plot,
       width = 8.3, height = 8.8, dpi = 300)

### Table D1: PC shares and seasonal characteristics across origins -------------

pc_feature_origin_all4 <- purrr::map_dfr(origins, function(T) {
  purrr::map_dfr(centers, function(C) {
    purrr::map_dfr(scales, function(S) {
      scores <- get_array_item(pc_actual_array, T, C, S)
      eigvals <- get_array_item(pc_eigval_array, T, C, S)
      month_factor <- factor((seq_len(nrow(scores)) - 1L) %% 12L + 1L)

      purrr::map_dfr(1:5, function(j) {
        x <- as.numeric(scores[, j])
        tibble::tibble(
          origin = T, center = C, scale = S, pc = j,
          component_share_pct = 100 * eigvals[j] / sum(eigvals),
          seasonal_R2_pct = 100 * summary(lm(x ~ month_factor))$r.squared
        )
      })
    })
  })
})

fmt_share <- function(x) sprintf("%.2f [%.2f, %.2f]", median(x, na.rm = TRUE),
                                 quantile(x, 0.025, na.rm = TRUE),
                                 quantile(x, 0.975, na.rm = TRUE))
fmt_r2 <- function(x) sprintf("%.1f [%.1f, %.1f]", median(x, na.rm = TRUE),
                              quantile(x, 0.025, na.rm = TRUE),
                              quantile(x, 0.975, na.rm = TRUE))

pc_desc_table_all4 <- pc_feature_origin_all4 %>%
  dplyr::mutate(
    Specification = dplyr::case_when(
      center & !scale ~ "Centred, unscaled",
      center & scale ~ "Centred, scaled",
      !center & !scale ~ "Uncentred, unscaled",
      !center & scale ~ "Uncentred, scaled"
    ),
    Specification = factor(Specification,
                           levels = c("Centred, unscaled", "Centred, scaled",
                                      "Uncentred, unscaled", "Uncentred, scaled"))
  ) %>%
  dplyr::group_by(Specification, pc) %>%
  dplyr::summarise(
    `Total sum-of-squares share (%)` = fmt_share(component_share_pct),
    `Seasonal R2 (%)` = fmt_r2(seasonal_R2_pct),
    .groups = "drop"
  ) %>%
  dplyr::mutate(PC = paste0("PC", pc)) %>%
  dplyr::select(Specification, PC, `Total sum-of-squares share (%)`, `Seasonal R2 (%)`)

centred_desc <- pc_desc_table_all4 %>%
  dplyr::filter(grepl("^Centred", Specification)) %>%
  dplyr::mutate(Specification = sub("^Centred, ", "", Specification)) %>%
  dplyr::rename(
    `Centred` = Specification,
    `Share (%)` = `Total sum-of-squares share (%)`,
    `Seasonal R2 (%)` = `Seasonal R2 (%)`
  )

uncentred_desc <- pc_desc_table_all4 %>%
  dplyr::filter(grepl("^Uncentred", Specification)) %>%
  dplyr::mutate(Specification = sub("^Uncentred, ", "", Specification)) %>%
  dplyr::rename(
    `Uncentred` = Specification,
    `PC ` = PC,
    `Share (%) ` = `Total sum-of-squares share (%)`,
    `Seasonal R2 (%) ` = `Seasonal R2 (%)`
  )

pc_desc_table_side_by_side <- dplyr::bind_cols(centred_desc, uncentred_desc)

saveRDS(pc_desc_table_side_by_side, "rds/table_D1_pc_descriptives.rds")
saveRDS(list(pc_feature_origin_all4 = pc_feature_origin_all4,
             figure_d2_definitions_initial = figure_d2$definitions,
             figure_d2_definitions_T192 = figure_d2_mid$definitions,
             figure_d2_definitions_T263 = figure_d2_end$definitions),
        "rds/descriptive_analysis_objects.rds")

pc_desc_table_side_by_side