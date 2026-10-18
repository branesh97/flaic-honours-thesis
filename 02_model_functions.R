# Forecasting with linear augmented information combination (FLAIC)
# Model functions used by the final empirical investigation

library(hts)
library(forecast)
library(glmnet)

### Functions ------------------------------------------------------------------------------------------------------------------------------------

#### Function for base forecasts + BU and conventional reconciliation, using ETS ---------------------------------------------------------------
# base forecasts are generated in the same function as reconciled forecasts for the non-IComb comparison methods

bfcst_ets <- function(T) {
  # start time to compile elapsed minutes per origin(T)
  start_time <- Sys.time()
  data_a <- ts(ally[(T-window_size+1):T,], start = c(1998 + (T-window_size) %/% seasonal_period, 1 + (T-window_size) %% seasonal_period), frequency = seasonal_period)
  # for each col (series), we need to fit ETS model
  resid_fitted <- matrix(NA_real_, nrow = nrow(data_a), ncol = ncol(data_a))
  colnames(resid_fitted) <- colnames(ally)
  bfcst <- matrix(NA_real_, nrow = h, ncol = ncol(data_a))
  colnames(bfcst) <- colnames(ally)
  rownames(bfcst) <- paste0("h", 1:h)
  fitted_mat <- matrix(NA_real_, nrow = nrow(data_a), ncol = ncol(data_a), dimnames = list(NULL, colnames(ally)))
  for (i in 1:ncol(data_a)) {
    fit_m <- ets(data_a[,i], ic = "aicc")
    fitted_mat[, i] <- fitted(fit_m)
    resid_fitted[, i] <- data_a[, i] - fitted(fit_m)
    bfcst[, i] <- forecast(fit_m, h = h)$mean
  }
  # bottom-up reconciliation
  bu_fcst <- bfcst[, bottom_idx] %*% t(S)
  # conventional reconciliation
  # weights accept a numeric vector which would be diagonal entries of the precision matrix (inverse of weighting matrix)
  ols_fcst <- combinef(bfcst, nodes = nodes, weights = NULL, keep = "all", algorithms = "lu") # would be same if we used weights as rep(1,112) instead
  wls_s_fcst <- combinef(bfcst, nodes = nodes, weights = w_struct, keep = "all", algorithms = "lu")
  wls_v_fcst <- combinef(bfcst, nodes = nodes, weights = 1/colMeans(resid_fitted^2), keep = "all", algorithms = "lu")
  mint_shr_fcst <- MinT(bfcst, nodes = nodes, residual = resid_fitted, covariance = "shr", keep = "all", algorithms = "lu")
  # result returned
  list(T = T, bfcst = bfcst, fitted_mat = fitted_mat, bu_fcst = bu_fcst, ols_fcst = ols_fcst, wls_s_fcst = wls_s_fcst,
       wls_v_fcst = wls_v_fcst, mint_shr_fcst = mint_shr_fcst, elapsed_minutes = as.numeric(difftime(Sys.time(), start_time, units = "mins")))
}

#### Functions for IComb ------------------------------------------------------------------------------------------------------------------------
# lambda_path and icomb_reg are very similar to functions in Shanika's GitHub IComb repository

##### Function: lambda_path ----------------------------------------------------------------------------------------------------------------------
# IComb paper had grid of 200, so will do the same
# alpha = 0: ridge and alpha = 1: lasso
# lambda is linear in log scale

lambda_path <- function(x, y, alpha = 1, standardize = FALSE, standardize_response = FALSE, intercept = TRUE, lambda_min_ratio = "expand", nlambda = 200) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  # your x is the X predictor matrix excluding intercept and rows are your time in training sample and columns are your predictors
  dimx <- dim(x) # first value will be no of obs in training sample and second value will be number of predictors (excluding intercept; depends which have been supplied)
  nobs <- dimx[1]
  nvars <- dimx[2]
  # your y is the response variable and it will have dimy[1] = nobs rows and ny columns (all series in ally will be in response)
  dimy <- dim(y)
  # if intercept = TRUE, then we centre the variables; intercept absorbs the column means and in glmnet it is not penalised
  if (intercept) {
    xm <- colMeans(x)
    ym <- colMeans(y)
  } else {
    xm <- rep(0, times = nvars)
    ym <- rep(0, times = dimy[2])
  }
  # if we standardise response (y), calculate response std dev; else set value 1
  if (standardize_response) {
    ysd <- sqrt(colMeans(scale(y, center = TRUE, scale = FALSE)^2))
  } else {
    ysd <- rep(1, times = dimy[2])
  }
  # scaling y; if intercept = TRUE, we subtract means and if standardize_response = TRUE, we divide by std dev
  y <- scale(y, center = ym, scale = ysd)
  # if we standardise predictors (x), calculate predictor std dev; else set value 1
  if (standardize) {
    xsd <- sqrt(colMeans(scale(x, center = TRUE, scale = FALSE)^2))
  } else {
    xsd <- rep(1, times = nvars)
  }
  # scaling x
  x <- scale(x, center = xm, scale = xsd)
  # compute lambda_max
  # we find which predictor has strongest combined relationship with response matrix
  # and then divide it by something related to alpha which is 1 for lasso and 0.001 for ridge (pmax(0,0.001) = 0.001)
  # for all slope coefficients to remain zero, the group-lasso penalty must be at least as large as every predictor's gradient
  lambda_max <- max(sqrt(rowSums(crossprod(x, y)^2))/nobs) / pmax(alpha, 0.001)
  # compute lambda_min
  # if lambda_min_ratio == "expand", we use as per IComb paper
  # if lambda_min_ratio == "glmnet", we use glmnet
  if (lambda_min_ratio == "expand")
    lambda_min_ratio <- 10^(-floor(log10(lambda_max))-2)
  else if (lambda_min_ratio == "glmnet")
    lambda_min_ratio <- ifelse(nobs < nvars, 0.01, 1e-04)
  lambda_min <- lambda_min_ratio * lambda_max
  # return sequence of lambdas depending on min and max, linear on log scale
  exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
}

##### Function: IComb with regularisation --------------------------------------------------------------------------------------------------------
# we have added thresh, but kept same as default in glmnet
# we also have added argument of window being rolling or expanding for user's choice during best lambda computation
# loss_scale option is available as per one version of IComb code which selects which loss we minimise during nested CV to select lambda
# we use "original" as per IComb paper, as noted in the paper

icomb_reg <- function(fitted, actual, forecasts, train_size, window = c("rolling", "expanding"), alpha = 1, standardize = FALSE,
                      standardize_response = FALSE, intercept = TRUE, lambda = NULL, lambda_min_ratio = "expand",
                      nlambda = 200, maxit = 1e+07, thresh = 1e-07, loss_scale = c("original", "variance", "snaive")) {
  window <- match.arg(window)
  loss_scale <- match.arg(loss_scale)
  # FPP seasonal-naive lambda selection can be added later
  if (loss_scale == "snaive") stop("loss_scale = 'snaive' is not implemented yet.")
  # predictor matrix excluding intercept; rows are time and columns are predictors
  dimx <- dim(fitted)
  nobs <- dimx[1]
  nvars <- dimx[2]
  # response matrix; rows are time and columns are hierarchy series
  dimy <- dim(actual)
  # dimension of forecasts; rows can be all h horizons and columns should match predictors
  dimfc <- dim(forecasts)
  # number of iterations over which we fit all lambdas and find lambda with lowest MSE
  niter <- nobs - train_size
  # if lambda was not specified, get lambda sequence from lambda_path
  if (is.null(lambda)) lambda <- lambda_path(fitted, actual, alpha, standardize, standardize_response, intercept, lambda_min_ratio, nlambda)
  # store reconciled forecasts for each internal validation iteration, response series and lambda
  recon <- array(NA_real_, dim = c(niter, dimy[2], length(lambda)))
  # actual validation observations
  test <- actual[(train_size + 1):nobs, , drop = FALSE]
  # response SD used for variance-standardising each internal validation error
  ysd_all <- matrix(NA_real_, nrow = niter, ncol = dimy[2], dimnames = list(NULL, colnames(actual)))
  for (i in 1:niter) {
    # IComb paper uses an expanding window
    if (window == "rolling") {
      train_set <- (1+i-1):(train_size+i-1)
    } else {
      train_set <- 1:(train_size+i-1)
    }
    xdata <- fitted[train_set, ]
    ydata <- actual[train_set, ]
    # response SDs based ONLY on this internal training sample
    ysd <- sqrt(colMeans(scale(ydata, center = TRUE, scale = FALSE)^2))
    if (any(ysd < sqrt(.Machine$double.eps))) stop("At least one response has approximately zero variance in internal CV iteration ", i)
    ysd_all[i, ] <- ysd
    fit <- glmnet(xdata, ydata, family = "mgaussian", standardize = standardize, standardize.response = standardize_response,
                  intercept = intercept, alpha = alpha, lambda = lambda, maxit = maxit, thresh = thresh)
    pred_i <- predict(fit, newx = fitted[train_size+i, , drop = FALSE])
    recon[i, ,] <- pred_i[1, , ]
  }
  # validation errors for every CV iteration, response and lambda
  err <- sweep(recon, MARGIN = c(1, 2), STATS = test, FUN = "-")
  if (loss_scale == "original") {
    loss_array <- err^2
  } else if (loss_scale == "variance") {
    err_scaled <- sweep(err, MARGIN = c(1, 2), STATS = ysd_all, FUN = "/")
    loss_array <- err_scaled^2
  }
  # average across the 40 validation origins and 112 response series, separately for each lambda
  mse_lambda <- apply(loss_array, MARGIN = 3, FUN = mean, na.rm = TRUE)
  idx_min <- which.min(mse_lambda)
  lambda_best <- lambda[idx_min]
  lambda_info <- list(lambda_max = max(lambda), lambda_min = min(lambda), lambda_best = lambda_best, lambda_best_idx = idx_min,
                      path_position = (idx_min - 1) / (length(lambda) - 1), at_lambda_max = idx_min == 1L, at_lambda_min = idx_min == length(lambda),
                      log10_span = log10(max(lambda) / min(lambda)))
  # with lambda going from max to min, fit from lambda_max through to lambda_best
  lambda_subset <- lambda[seq_len(idx_min)]
  fit_final <- glmnet(fitted, actual, family = "mgaussian", standardize = standardize, standardize.response = standardize_response, intercept = intercept,
                      alpha = alpha, lambda = lambda_subset, maxit = maxit, thresh = thresh)
  coefs <- coef(fit_final)
  pred_array <- predict(fit_final, newx = forecasts, s = lambda_best)
  predictions <- matrix(pred_array[, , 1], nrow = nrow(forecasts), ncol = ncol(actual), dimnames = list(rownames(forecasts), colnames(actual)))
  # save final coefficient matrix at lambda_best, including the unregularised intercept
  best_coef <- do.call(cbind, lapply(coefs, function(x) as.numeric(x[, ncol(x)])))
  rownames(best_coef) <- rownames(coefs[[1]])
  colnames(best_coef) <- colnames(actual)
  # check that the saved coefficient matrix reproduces the glmnet prediction at lambda_best
  coef_predictions <- cbind("(Intercept)" = 1, as.matrix(forecasts)) %*% best_coef
  stopifnot(isTRUE(all.equal(unname(predictions), unname(coef_predictions), tolerance = 1e-08, check.attributes = FALSE)))
  list(preds = predictions, info = list(lambda_info = lambda_info, lambda_path = lambda, mse_info = mse_lambda, loss_scale = loss_scale,
                                         nnzeros = fit_final$dfmat[1, ncol(fit_final$dfmat)]), coefs = best_coef)
}

##### Function: Obtaining IComb regularised reconciled forecasts --------------------------------------------------------------------------------
# we only consider the case where X and Y are standardised

icomb_rec <- function(T, loss_scale) {
  # start time to compile elapsed minutes per origin(T)
  start_time <- Sys.time()
  # actual for IComb
  Y <- as.matrix(ally[(T-window_size+1):T, , drop = FALSE])
  # predictors or fitted (no intercept based on github code)
  X <- fitted_array[as.character(T), , ]
  # obtain base forecast computed earlier
  bfcst <- fcst_array[as.character(T), , , "base"]
  # ridge, standardise X and Y, with intercept
  icombs_xy1_res <- icomb_reg(fitted = X, actual = Y, forecasts = bfcst, train_size = inner_train_size, window = "expanding", alpha = 0, standardize = TRUE,
                              standardize_response = TRUE, intercept = TRUE, lambda = NULL, lambda_min_ratio = "expand", nlambda = nlambda, maxit = 1e+07, loss_scale = loss_scale)
  # lasso, standardise X and Y, with intercept
  icombl_xy1_res <- icomb_reg(fitted = X, actual = Y, forecasts = bfcst, train_size = inner_train_size, window = "expanding", alpha = 1, standardize = TRUE,
                              standardize_response = TRUE, intercept = TRUE, lambda = NULL, lambda_min_ratio = "expand", nlambda = nlambda, maxit = 1e+07, loss_scale = loss_scale)
  # result returned
  list(T = T, icombs_xy1_res = icombs_xy1_res, icombl_xy1_res = icombl_xy1_res, loss_scale = loss_scale,
       elapsed_minutes = as.numeric(difftime(Sys.time(), start_time, units = "mins")))
}

#### Functions for FLAIC ------------------------------------------------------------------------------------------------------------------------

##### Function: Generating principal components -------------------------------------------------------------------------------------------------
# function depends on:
# origin (T) as data from each training window to origin used to generate PCs
# center = TRUE or FALSE (whether to de-mean prior to generating PCs)
# scale = TRUE or FALSE (whether to equalise variance prior to generating PCs)
# have set retx = TRUE because the returned scores are used as the constructed component series

gen_pc <- function(T, center = FALSE, scale = FALSE, use_pc_scores = TRUE) {
  # start time to compile elapsed minutes per origin(T)
  start_time <- Sys.time()
  data_a <- ts(ally[(T-window_size+1):T,], start = c(1998 + (T-window_size) %/% seasonal_period, 1 + (T-window_size) %% seasonal_period), frequency = seasonal_period)
  # note: subtract 1 from time periods for rank since centering
  rank_pc <- min(nb, nrow(data_a)-1)
  pca_res <- prcomp(x = data_a, retx = TRUE, center = center, scale. = scale, rank. = rank_pc)
  # actual PC scores (data is demeaned and/or scaled prior to multiplying with rotation matrix)
  pc_scores <- pca_res$x
  # if use_pc_scores = TRUE, then the constructed components are PC scores
  # else, the constructed components are data_a %*% pca_res$rotation
  if (use_pc_scores == TRUE) {
    pca_comp <- ts(pc_scores, start = start(data_a), frequency = frequency(data_a))
  } else {
    pca_comp <- ts(data_a %*% pca_res$rotation, start = start(data_a), frequency = frequency(data_a))
  }
  # correlation matrix with constructed components
  cor_mat_lc <- cor(data_a, pca_comp)
  colnames(cor_mat_lc) <- paste0("PC", 1:ncol(cor_mat_lc))
  rownames(cor_mat_lc) <- colnames(ally)
  # correlation matrix with PC scores
  cor_mat_pc <- cor(data_a, pc_scores)
  colnames(cor_mat_pc) <- paste0("PC", 1:ncol(cor_mat_pc))
  rownames(cor_mat_pc) <- colnames(ally)
  # fit ETS model to each component and save fitted values and base forecasts
  pc_fitted_matrix <- matrix(NA_real_, nrow = nrow(pca_comp), ncol = ncol(pca_comp), dimnames = list(rownames(pca_comp), colnames(pca_comp)))
  pc_fcst_matrix <- matrix(NA_real_, nrow = h, ncol = ncol(pca_comp), dimnames = list(paste0("h", 1:h), colnames(pca_comp)))
  for (i in seq_len(ncol(pca_comp))) {
    fit_m <- ets(pca_comp[,i], ic = "aicc")
    pc_fitted_matrix[, i] <- fitted(fit_m)
    pc_fcst_matrix[, i] <- forecast(fit_m, h = h)$mean
  }
  # eigenvalues, proportional variance and cumulative variance
  eigval <- pca_res$sdev^2
  propvar <- eigval / sum(eigval)
  cumvar <- cumsum(propvar)
  # result returned
  list(T = T, center = center, scale = scale, use_pc_scores = use_pc_scores, eigval = eigval, eigvec = pca_res$rotation, propvar = propvar, cumvar = cumvar,
       cor_mat_lc = cor_mat_lc, cor_mat_pc = cor_mat_pc, pca_comp = pca_comp, pc_scores = pc_scores, pc_fitted_matrix = pc_fitted_matrix,
       pc_fcst_matrix = pc_fcst_matrix, elapsed_minutes = as.numeric(difftime(Sys.time(), start_time, units = "mins")))
}

##### Function: FLAIC reconciliation by origin (T) and number of principal components (p) --------------------------------------------------------
# X and Y are standardised when computing regularised IComb coefficients, the intercept is included but not regularised,
# and lambda is selected over the final 40 observations using an expanding internal training window

flaic_rec <- function(T, p, center, scale, loss_scale) {
  # start time to compile elapsed minutes per origin(T)-p combination
  start_time <- Sys.time()
  # actual for IComb
  Y <- as.matrix(ally[(T-window_size+1):T, , drop = FALSE])
  # predictors or fitted (original hierarchy fitted values augmented with fitted PC values)
  X <- cbind(fitted_array[as.character(T), , ], pc_fitted_array[[as.character(T), as.character(center), as.character(scale)]][, seq_len(p), drop = FALSE])
  # obtain base forecasts (original hierarchy base forecasts augmented with PC forecasts)
  bfcst <- cbind(fcst_array[as.character(T), , , "base"], pc_fcst_array[[as.character(T), as.character(center), as.character(scale)]][, seq_len(p), drop = FALSE])
  # FLAIC: ridge, standardise X and Y, with intercept
  flaics_xy1_res <- icomb_reg(fitted = X, actual = Y, forecasts = bfcst, train_size = inner_train_size, window = "expanding", alpha = 0, standardize = TRUE,
                              standardize_response = TRUE, intercept = TRUE, lambda = NULL, lambda_min_ratio = "expand", nlambda = nlambda, maxit = 1e+07, loss_scale = loss_scale)
  # FLAIC: lasso, standardise X and Y, with intercept
  flaicl_xy1_res <- icomb_reg(fitted = X, actual = Y, forecasts = bfcst, train_size = inner_train_size, window = "expanding", alpha = 1, standardize = TRUE,
                              standardize_response = TRUE, intercept = TRUE, lambda = NULL, lambda_min_ratio = "expand", nlambda = nlambda, maxit = 1e+07, loss_scale = loss_scale)
  # result returned
  list(T = T, p = p, center = center, scale = scale, flaics_xy1_res = flaics_xy1_res, flaicl_xy1_res = flaicl_xy1_res, loss_scale = loss_scale,
       elapsed_minutes = as.numeric(difftime(Sys.time(), start_time, units = "mins")))
}

#### Function for FLAP benchmark -----------------------------------------------------------------------------------------------------------------
# Required existing objects:
# ally, origins, centers, scales, window_size, h, ny, nb, S, bottom_idx,
# fitted_array, fcst_array, pc_scores_array, pc_fitted_array, pc_fcst_array and pc_eigvec_array
# This function produces FLAP forecasts only; forecast evaluation is kept in 05_evaluate.R with the other methods.

run_flap_benchmark <- function(p_flap, verbose = TRUE) {
  # 0. Checks / setup
  p_flap <- as.integer(p_flap)
  stopifnot(length(p_flap) == 1L, p_flap >= 1L, p_flap <= nb)
  n_aug <- ny + p_flap
  if (verbose) message("Starting FLAP benchmark with p = ", p_flap, " | augmented dimension = ", n_aug)
  # 1. Storage
  fcst_array_star <- array(NA_real_, dim = c(origin = length(origins), horizon = h, series = ny, center = length(centers), scale = length(scales)),
                           dimnames = list(origin = as.character(origins), horizon = paste0("h", seq_len(h)), series = unname(colnames(ally)),
                                           center = as.character(centers), scale = as.character(scales)))
  # 2. Generate FLAP forecasts
  for (T in origins) {
    T_chr <- as.character(T)
    # 2.1 Actual hierarchy over rolling window
    Y_T <- as.matrix(ally[(T - window_size + 1):T, , drop = FALSE])
    B_T <- Y_T[ , bottom_idx, drop = FALSE]
    stopifnot(all(dim(Y_T) == c(window_size, ny)), all(dim(B_T) == c(window_size, nb)))
    # 2.2 Original hierarchy fitted values
    y_fitted_T <- as.matrix(fitted_array[T_chr, , , drop = TRUE])
    if (length(dim(y_fitted_T)) == 2L && all(dim(y_fitted_T) == c(ny, window_size))) y_fitted_T <- t(y_fitted_T)
    stopifnot(all(dim(y_fitted_T) == c(window_size, ny)), all(is.finite(y_fitted_T)))
    # 2.3 Original hierarchy base forecasts
    y_hat_T <- as.matrix(fcst_array[T_chr, , , "base", drop = TRUE])
    if (length(dim(y_hat_T)) == 2L && all(dim(y_hat_T) == c(ny, h))) y_hat_T <- t(y_hat_T)
    stopifnot(all(dim(y_hat_T) == c(h, ny)), all(is.finite(y_hat_T)))
    # 2.4 Center-scale combinations
    for (center in centers) {
      center_chr <- as.character(center)
      for (scale in scales) {
        scale_chr <- as.character(scale)
        # exact preprocessing used by prcomp
        preprocessed_T <- scale(Y_T, center = center, scale = scale)
        if (center) center_vec_T <- as.numeric(attr(preprocessed_T, "scaled:center")) else center_vec_T <- rep(0, ny)
        if (scale) scale_vec_T <- as.numeric(attr(preprocessed_T, "scaled:scale")) else scale_vec_T <- rep(1, ny)
        stopifnot(length(center_vec_T) == ny, length(scale_vec_T) == ny, all(is.finite(center_vec_T)), all(is.finite(scale_vec_T)), all(scale_vec_T > sqrt(.Machine$double.eps)))
        # V_T: take only first p_flap PCs; dimensions: 112 x p_flap
        V_T <- pc_eigvec_array[[T_chr, center_chr, scale_chr]][ , seq_len(p_flap), drop = FALSE]
        stopifnot(all(dim(V_T) == c(ny, p_flap)))
        # location adjustment: V' D^-1 mu, where mu = ybar if centered and 0 otherwise
        adjustment_T <- as.numeric(t(V_T) %*% (center_vec_T / scale_vec_T))
        # actual v
        pc_actual_T <- as.matrix(pc_scores_array[[T_chr, center_chr, scale_chr]][ , seq_len(p_flap), drop = FALSE])
        v_actual_T <- sweep(pc_actual_T, MARGIN = 2, STATS = adjustment_T, FUN = "+")
        # fitted v
        pc_fitted_T <- as.matrix(pc_fitted_array[[T_chr, center_chr, scale_chr]][ , seq_len(p_flap), drop = FALSE])
        v_fitted_T <- sweep(pc_fitted_T, MARGIN = 2, STATS = adjustment_T, FUN = "+")
        # forecast v
        pc_hat_T <- as.matrix(pc_fcst_array[[T_chr, center_chr, scale_chr]][ , seq_len(p_flap), drop = FALSE])
        if (length(dim(pc_hat_T)) == 2L && all(dim(pc_hat_T) == c(p_flap, h))) pc_hat_T <- t(pc_hat_T)
        stopifnot(all(dim(pc_hat_T) == c(h, p_flap)))
        v_hat_T <- sweep(pc_hat_T, MARGIN = 2, STATS = adjustment_T, FUN = "+")
        v_names <- paste0("v_PC", seq_len(p_flap))
        colnames(v_actual_T) <- v_names
        colnames(v_fitted_T) <- v_names
        colnames(v_hat_T) <- v_names
        # S*_T = [V' D^-1 S; S]
        D_inv_S_T <- sweep(S, MARGIN = 1, STATS = scale_vec_T, FUN = "/")
        S_v_T <- t(V_T) %*% D_inv_S_T
        S_star_T <- rbind(S_v_T, S)
        rownames(S_star_T) <- c(v_names, unname(colnames(ally)))
        colnames(S_star_T) <- unname(colnames(S))
        stopifnot(all(dim(S_star_T) == c(n_aug, nb)))
        # algebra check 1: v = Y D^-1 V
        v_direct_T <- sweep(Y_T, MARGIN = 2, STATS = scale_vec_T, FUN = "/") %*% V_T
        if (!isTRUE(all.equal(v_actual_T, v_direct_T, tolerance = 1e-7, check.attributes = FALSE))) stop("v = Y D^-1 V failed: T=", T, ", C=", center, ", S=", scale, ", p=", p_flap)
        # algebra check 2: u_actual = B S*'
        u_actual_T <- cbind(v_actual_T, Y_T)
        u_structural_T <- B_T %*% t(S_star_T)
        if (!isTRUE(all.equal(u_actual_T, u_structural_T, tolerance = 1e-7, check.attributes = FALSE))) stop("u = B S*' failed: T=", T, ", C=", center, ", S=", scale, ", p=", p_flap)
        # augmented fitted values and errors
        u_fitted_T <- cbind(v_fitted_T, y_fitted_T)
        u_error_T <- u_actual_T - u_fitted_T
        stopifnot(all(dim(u_error_T) == c(window_size, n_aug)), all(is.finite(u_error_T)))
        # shrinkage covariance consistent with MinT-shrink using hts:::shrink.estim
        tar <- hts:::lowerD(u_error_T)
        shrink_res_T <- hts:::shrink.estim(u_error_T, tar)
        W_star_T <- shrink_res_T[[1]]
        lambda_shr_T <- as.numeric(shrink_res_T[[2]][2])
        stopifnot(is.matrix(W_star_T), is.numeric(W_star_T), all(dim(W_star_T) == c(n_aug, n_aug)), all(is.finite(W_star_T)),
                  is.finite(lambda_shr_T), lambda_shr_T >= 0, lambda_shr_T <= 1)
        # compute G*_T directly -- no need to explicitly form W*^{-1}
        W_inv_S_T <- solve(W_star_T, S_star_T)
        middle_T <- crossprod(S_star_T, W_inv_S_T)
        G_star_T <- solve(middle_T, t(W_inv_S_T))
        stopifnot(all(dim(G_star_T) == c(nb, n_aug)), all(is.finite(G_star_T)))
        if (!isTRUE(all.equal(G_star_T %*% S_star_T, diag(nb), tolerance = 1e-6, check.attributes = FALSE))) stop("G*S*=I failed: T=", T, ", C=", center, ", S=", scale, ", p=", p_flap)
        # augmented forecasts
        u_hat_T <- cbind(v_hat_T, y_hat_T)
        stopifnot(all(dim(u_hat_T) == c(h, n_aug)))
        # reconciled bottom forecasts and rebuild coherent original hierarchy
        b_tilde_T <- u_hat_T %*% t(G_star_T)
        y_tilde_T <- b_tilde_T %*% t(S)
        colnames(y_tilde_T) <- unname(colnames(ally))
        rownames(y_tilde_T) <- paste0("h", seq_len(h))
        stopifnot(all(dim(y_tilde_T) == c(h, ny)), all(is.finite(y_tilde_T)))
        # final coherence
        y_bottom_T <- y_tilde_T[ , bottom_idx, drop = FALSE]
        if (!isTRUE(all.equal(y_tilde_T, y_bottom_T %*% t(S), tolerance = 1e-7, check.attributes = FALSE))) stop("Final coherence failed: T=", T, ", C=", center, ", S=", scale, ", p=", p_flap)
        # store
        fcst_array_star[T_chr, , , center_chr, scale_chr] <- y_tilde_T
      }
    }
    if (verbose && match(T, origins) %% 12L == 0L) message("p = ", p_flap, " completed through origin ", T, " | ", match(T, origins), "/", length(origins))
  }
  # result returned
  list(p_flap = p_flap, n_aug = n_aug, fcst_array_star = fcst_array_star)
}