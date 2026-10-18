# Forecasting with linear augmented information combination (FLAIC)

library(tidyverse)
library(hts)

### Main settings -------------------------------------------------------------------------------------------------------------------------------

# forecast horizon and seasonal frequency are both 12 for monthly data, but are kept as separate parameters
h <- 12L
seasonal_period <- 12L
# rolling estimation window is 120 months (10 years)
window_size <- 120L
# regularised IComb / FLAIC use the last 40 observations of the outer training window for internal validation
validation_size <- 40L
inner_train_size <- window_size - validation_size
# regularisation grid follows the IComb implementation
nlambda <- 200L
# principal-component augmentation dimensions considered in the final empirical investigation
p_values <- c(1, 2, 3, 5, 10, 20, 40, 77)
penalties <- c("ridge", "lasso")
centers <- c(TRUE, FALSE)
scales <- c(TRUE, FALSE)
# loss used for internal lambda selection in regularised IComb / FLAIC
lambda_loss_scale <- "original"
# whether we use PC scores as our constructed component values
use_pc_scores <- TRUE
# parallel workers
workers <- 8L
# checkpoint batch size for the base/IComb/PC stages
checkpoint_every <- 12L

# folders used by the final run
dir.create("checkpoints", showWarnings = FALSE)
dir.create("coefficients", showWarnings = FALSE)
dir.create("rds", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("supplementary", showWarnings = FALSE)

### Getting data ready ---------------------------------------------------------------------------------------------------------------------------

tourism <- read.csv(file.path("data", "TRA-final.csv"), header = TRUE)

# extract the bottom-level time series (bts) with purpose of travel
bts_p <- tourism[, 3:ncol(tourism)]
bts_p <- ts(bts_p, start = c(1998, 1), frequency = seasonal_period)
# there are 308 series which is 77 regions (bts we need + purposes of: holiday, visiting friends & relatives, business, other)
# so we need to aggregate 4 columns by each row to get just region level data in each column

# for obtaining regional level data
# we want to sum col 1:4, 5:8, ..., 301:304, 305:308 by each row
idx_s <- seq(1, ncol(bts_p), by = 4) # gives 1, 5, ..., 301, 305
idx_e <- seq(4, ncol(bts_p), by = 4) # gives 4, 8, ..., 304, 308
colnames_bts <- substr(colnames(bts_p), 1, 3)[idx_s] # extracts the region identifier from the first purpose column for each region
bts <- matrix(, nrow(bts_p), ncol(bts_p)/4) # empty 264 by 77 matrix
# summation of relevant columns by each row happens here
for (i in 1:length(colnames_bts)) {
  bts[, i] <- rowSums(bts_p[, idx_s[i]:idx_e[i]])
}
colnames(bts) <- colnames_bts
bts <- ts(bts, start = c(1998, 1), frequency = seasonal_period)

# create hierarchical time series object from regional level data
tourism_hts <- hts(bts, characters = c(1, 1, 1))

# clean-up
rm(bts, bts_p, tourism, idx_e, idx_s, colnames_bts)

# create multivariate time series (y) object since hts stores data by bottom series
ally <- allts(tourism_hts)
ny <- ncol(ally) # number of series in y = 112
nT <- nrow(ally) # number of observations in time for each series in y = 264 months (22 years)
nb <- ncol(tourism_hts[["bts"]]) # number of bottom series in y = 77
na <- ny - nb # number of non-bottom series in y = 35
S <- smatrix(tourism_hts) # summing matrix
colnames(S) <- colnames(tourism_hts$bts)
rownames(S) <- colnames(ally)
bottom_idx <- (ny - nb + 1):ny
group_idx <- list(Australia = 1, States = 2:8, Zones = 9:35, Regions = 36:112, Average = 1:112)
# labels used for reporting
state_label_lookup <- c(A = "NSW", B = "VIC", C = "QLD", D = "SA", E = "WA", F = "TAS", G = "NT")
series_label_lookup <- c(Total = "Australia", state_label_lookup)
nodes <- get_nodes(tourism_hts) # needed in combinef function of hts for reconciliation
w_struct <- 1 / rowSums(S) # needed for WLS_s reconciliation
# forecast origins are T = 120 to T = 263 for the 264-month sample
origins <- window_size:(nT - 1L)

### Save data and settings ------------------------------------------------------------------------------------------------------------------------

data_setup <- list(tourism_hts = tourism_hts, ally = ally, ny = ny, nT = nT, nb = nb, na = na, S = S, bottom_idx = bottom_idx, group_idx = group_idx,
                   state_label_lookup = state_label_lookup, series_label_lookup = series_label_lookup, nodes = nodes, w_struct = w_struct,
                   h = h, seasonal_period = seasonal_period, window_size = window_size, validation_size = validation_size,
                   inner_train_size = inner_train_size, nlambda = nlambda, origins = origins, p_values = p_values, penalties = penalties,
                   centers = centers, scales = scales, lambda_loss_scale = lambda_loss_scale, use_pc_scores = use_pc_scores,
                   workers = workers, checkpoint_every = checkpoint_every)
saveRDS(data_setup, "rds/flaic_01_data_setup.rds")
message("01_data_setup complete: data and final-run settings saved.")