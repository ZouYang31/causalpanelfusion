#' Prediction intervals for synthetic control data fusion
#'
#' Constructs prediction intervals for the time-averaged effect
#' \eqn{\psi_0} and the ingredients for per-period bands ([pi_bands()]),
#' following Section 4.1 of the paper. The estimation error is decomposed as
#' \deqn{\hat\psi - \psi_0 = \Delta_w + e,}
#' where \eqn{\Delta_w} is the weight-construction (in-sample) error from
#' fitting the weights in the reference domain and \eqn{e} is the transfer
#' (out-of-sample) error in the target domain. The error probability is split
#' as \eqn{\alpha = \alpha_1 + \alpha_2}.
#'
#' @section Weight-construction error (\eqn{M_1}):
#' Residuals of the reference-domain fit are centred at their mean and
#' resampled (moving-block or i.i.d.) to create `B` bootstrap reference
#' trajectories for the treated unit. For each, the full algorithm of
#' [scfusion()] is re-run (including budget re-selection when
#' `reselect_budget = TRUE`) and \eqn{\delta_b = \hat C - \hat C_b} is
#' recorded, where \eqn{\hat C = \sum_j w_j \bar Y_j}. \eqn{M_{1L}, M_{1U}}
#' are the \eqn{\alpha_1/2} and \eqn{1 - \alpha_1/2} quantiles of
#' \eqn{\delta_b}.
#'
#' @section Transfer error (\eqn{M_2}):
#' There is no target-domain pre-period, so the transfer error is calibrated
#' from donor placebos: each calibration donor (by default the donors with
#' positive weight, the feasible analogue of the pseudo-true donor set) is
#' treated as the target and refitted from the remaining donors (the treated
#' unit is always excluded). Their time-averaged target-domain placebo gaps,
#' centred at their mean, give the plug-in scale \eqn{\hat\sigma} (sample
#' standard deviation). Three constructions are reported:
#' \describe{
#'   \item{Approach 2, Gaussian}{\eqn{M_2 = \hat\sigma\,(z_{\alpha_2/2},
#'     z_{1-\alpha_2/2})}.}
#'   \item{Approach 2, empirical}{\eqn{M_2 = \hat\sigma\,(\hat q_{\alpha_2/2},
#'     \hat q_{1-\alpha_2/2})} with quantiles of the standardised placebo
#'     gaps.}
#'   \item{Approach 1, sub-Gaussian}{\eqn{M_2 = \pm\nu\sqrt{2\log(2/\alpha_2)}}
#'     with \eqn{\nu = \hat\sigma} or a user-supplied upper bound.}
#' }
#' The interval is \eqn{[\hat\psi - M_{1U} - M_{2U} - \gamma,\ \hat\psi -
#' M_{1L} - M_{2L} + \gamma]}, where \eqn{\gamma} is `bias_allowance`.
#'
#' @section Interpretation:
#' These are feasible plug-in analogues of the oracle intervals: estimation
#' error in the scale and in the empirical quantiles is not propagated, and
#' Assumption 6 of the paper is maintained rather than verified. With
#' `bias_allowance = 0` the interval targets \eqn{\psi_0 + E[e \mid H]};
#' supply a valid bound on \eqn{|E[e \mid H]|} to target \eqn{\psi_0}.
#'
#' @param fit An `"scfusion"` object from [scfusion()].
#' @param alpha1,alpha2 Error probabilities for the weight-construction and
#'   transfer components. The nominal level is `1 - alpha1 - alpha2`.
#' @param B Number of bootstrap refits.
#' @param resample `"mbb"` (moving-block bootstrap, default) or `"iid"`.
#' @param block_len Block length for the moving-block bootstrap.
#' @param reselect_budget If `TRUE` (default) every refit re-runs the full
#'   budget search; if `FALSE` the budget selected by `fit` is held fixed (a
#'   faster, conditional approximation).
#' @param budget_step Optional coarser budget grid spacing for the refits.
#'   Defaults to the grid used by `fit`.
#' @param placebo_pool Calibration donors: `"positive_weight"` (default) or
#'   `"all"` donors (a robustness check).
#' @param placebo_weight_tolerance Weight above which a donor counts as
#'   positive-weight.
#' @param empirical_quantile_type Quantile type for the empirical method
#'   (default `1`, the inverse empirical CDF); see [stats::quantile()].
#' @param common_scale Optional externally justified common transfer scale
#'   for Approach 2 (replaces the plug-in \eqn{\hat\sigma}).
#' @param scale_upper_bound Optional certified upper bound on the
#'   sub-Gaussian scale for Approach 1.
#' @param bias_allowance Non-negative bound \eqn{\gamma} on the conditional
#'   transfer bias (default `0`).
#' @param seed Optional seed for the bootstrap draws. The caller's random
#'   number stream is left untouched when a seed is given.
#' @param parallel Run the refits through \pkg{future.apply}. Choose a
#'   backend first, e.g. `future::plan(future::multisession)`. Results are
#'   identical to the serial run.
#' @param progress Show progress messages.
#'
#' @return An object of class `"scfusion_pi"` with elements `intervals` (data
#'   frame of average-effect intervals), `estimate`, `M1`, `scale`,
#'   `bootstrap` (draws `delta`, weights `weights` and diagnostics),
#'   `placebo` (donor-by-period placebo gaps), `settings` and `fit`. Use
#'   [pi_bands()] for per-period bands, `plot()` to draw them, and
#'   [sensitivity_pi()] to vary the transfer scale.
#'
#' @references Cattaneo, M. D., Feng, Y. and Titiunik, R. (2021). Prediction
#'   intervals for synthetic control methods. *Journal of the American
#'   Statistical Association*, 116, 1865-1880.
#'
#' @examples
#' \donttest{
#' sim <- simulate_fusion_data(n_units = 16, n_ref = 12, n_target = 4,
#'                             seed = 1)
#' fit <- scfusion(sim$target, sim$reference, X = sim$X, Z = sim$Z,
#'                 budget_step = 0.1)
#' pi <- scfusion_pi(fit, B = 50, seed = 1)
#' pi
#' pi_bands(pi, method = "gaussian")
#' plot(pi)
#' }
#' @export
scfusion_pi <- function(fit, alpha1 = 0.05, alpha2 = 0.05, B = 500L,
                        resample = c("mbb", "iid"), block_len = 3L,
                        reselect_budget = TRUE, budget_step = NULL,
                        placebo_pool = c("positive_weight", "all"),
                        placebo_weight_tolerance = 1e-6,
                        empirical_quantile_type = 1L,
                        common_scale = NULL, scale_upper_bound = NULL,
                        bias_allowance = 0, seed = NULL, parallel = FALSE,
                        progress = interactive()) {
  if (!inherits(fit, "scfusion")) {
    stop("`fit` must be an object returned by scfusion().", call. = FALSE)
  }
  .check_alphas(alpha1, alpha2)
  .check_scalar(B, "B", lower = 2, integer = TRUE)
  resample <- match.arg(resample)
  .check_scalar(block_len, "block_len", lower = 1, integer = TRUE)
  .check_flag(reselect_budget, "reselect_budget")
  placebo_pool <- match.arg(placebo_pool)
  .check_scalar(placebo_weight_tolerance, "placebo_weight_tolerance",
                lower = 0)
  .check_scalar(empirical_quantile_type, "empirical_quantile_type",
                lower = 1, upper = 9, integer = TRUE)
  .check_scalar(bias_allowance, "bias_allowance", lower = 0)
  if (!is.null(common_scale)) {
    .check_scalar(common_scale, "common_scale", lower = 0, strict = TRUE)
  }
  if (!is.null(scale_upper_bound)) {
    .check_scalar(scale_upper_bound, "scale_upper_bound", lower = 0)
  }
  .check_flag(parallel, "parallel")
  .check_flag(progress, "progress")
  if (fit$settings$solver == "gurobi") .require_gurobi()

  data <- fit$data
  settings <- fit$settings
  budgets <- .refit_budgets(fit, reselect_budget, budget_step)
  tr <- fit$treated
  dn <- fit$donors
  w <- as.numeric(fit$weights)
  Y <- data$target
  Fm <- data$reference

  # ---- point estimate and quantities held fixed --------------------------
  Y_treated <- Y[tr, ]
  Y_control <- Y[dn, , drop = FALSE]
  Ybar_control <- rowMeans(Y_control)
  Chat <- sum(w * Ybar_control)
  psi_hat <- mean(Y_treated) - Chat
  Chat_t <- as.numeric(crossprod(Y_control, w))
  synthetic_F <- as.numeric(crossprod(Fm[dn, , drop = FALSE], w))
  u_resid <- Fm[tr, ] - synthetic_F
  u_bar <- mean(u_resid)
  u_c <- u_resid - u_bar

  # ---- bootstrap draws (all randomness happens here) -----------------------
  draws <- .with_seed(seed, lapply(seq_len(B), function(b) {
    synthetic_F + u_bar + .block_resample(u_c, resample, block_len)
  }))

  # ---- weight-construction error -------------------------------------------
  boot_fits <- .cpf_lapply(draws, function(Fb) {
    tryCatch(
      .fit_unit(data, tr, dn, settings, budgets, F_treated = Fb)$weights,
      error = function(e) NULL
    )
  }, parallel = parallel, progress = progress,
  label = "Bootstrapping the weight-construction error")
  ok <- vapply(boot_fits, function(wb) {
    !is.null(wb) && length(wb) == length(dn) && all(is.finite(wb))
  }, logical(1))
  if (sum(ok) < 2L) {
    stop("Fewer than two bootstrap refits succeeded.", call. = FALSE)
  }
  W_boot <- do.call(rbind, boot_fits[ok])
  colnames(W_boot) <- data$unit_names[dn]
  delta <- Chat - as.numeric(W_boot %*% Ybar_control)
  weight_distance <- sqrt(rowSums(sweep(W_boot, 2L, w)^2))
  boot_diag <- list(
    requested = B, usable = sum(ok), failed = sum(!ok),
    distinct_weights = nrow(unique(round(W_boot, 8))),
    median_weight_distance = stats::median(weight_distance),
    sd_delta = stats::sd(delta)
  )
  if (boot_diag$distinct_weights == 1L ||
      boot_diag$sd_delta <= sqrt(.Machine$double.eps)) {
    warning("Bootstrap refits produced essentially no variation in the ",
            "weights; the in-sample component is degenerate.", call. = FALSE)
  }

  # ---- donor placebos for the transfer error -------------------------------
  placebo <- .donor_placebos(fit, budgets, placebo_pool,
                             placebo_weight_tolerance, parallel, progress)

  # ---- assemble ------------------------------------------------------------
  out <- structure(
    list(
      estimate = psi_hat,
      Chat = Chat,
      Chat_t = Chat_t,
      Y_treated = Y_treated,
      Y_control = Y_control,
      bootstrap = c(list(delta = delta, weights = W_boot), boot_diag),
      placebo = placebo,
      settings = list(
        alpha1 = alpha1, alpha2 = alpha2, B = B, resample = resample,
        block_len = block_len, reselect_budget = reselect_budget,
        n_budgets = nrow(budgets), placebo_pool = placebo_pool,
        placebo_weight_tolerance = placebo_weight_tolerance,
        empirical_quantile_type = empirical_quantile_type,
        common_scale = common_scale, scale_upper_bound = scale_upper_bound,
        bias_allowance = bias_allowance, seed = seed
      ),
      fit = fit
    ),
    class = "scfusion_pi"
  )
  avg <- .average_intervals(out, alpha1, alpha2, bias_allowance)
  out$intervals <- avg$intervals
  out$M1 <- avg$M1
  out$scale <- avg$scale
  out
}

.check_alphas <- function(alpha1, alpha2) {
  .check_scalar(alpha1, "alpha1", lower = 0, upper = 1, strict = TRUE)
  .check_scalar(alpha2, "alpha2", lower = 0, upper = 1, strict = TRUE)
  if (alpha1 + alpha2 >= 1) {
    stop("`alpha1 + alpha2` must be below 1.", call. = FALSE)
  }
}

.refit_budgets <- function(fit, reselect_budget, budget_step) {
  if (!reselect_budget) {
    b <- matrix(fit$budget, nrow = 1L)
    colnames(b) <- c("b_F", "b_Z", "b_X")
    return(b)
  }
  if (is.null(budget_step)) return(fit$settings$budgets)
  min_value <- fit$settings$budget_min
  if (is.na(min_value)) min_value <- 0.01
  budget_grid(budget_step, min_value)
}

# Moving-block (or i.i.d.) resample of x with the same length; identical to
# the research implementation (mbb_draw) so seeds reproduce its draws.
.block_resample <- function(x, resample, block_len) {
  n <- length(x)
  L <- as.integer(block_len)
  if (resample == "iid" || n <= L) return(sample(x, n, replace = TRUE))
  nblock <- ceiling(n / L)
  starts <- sample(seq_len(n - L + 1L), nblock, replace = TRUE)
  as.vector(vapply(starts, function(s) x[s:(s + L - 1L)],
                   numeric(L)))[seq_len(n)]
}

# Donor placebo fits in the target domain (treated unit always excluded).
.donor_placebos <- function(fit, budgets, pool, tol, parallel, progress) {
  data <- fit$data
  tr <- fit$treated
  dn <- fit$donors
  pool_idx <- if (pool == "positive_weight") dn[fit$weights > tol] else dn
  if (length(pool_idx) < 2L) {
    stop("Fewer than two calibration donors (positive-weight donors); the ",
         "transfer scale cannot be estimated. Try placebo_pool = \"all\".",
         call. = FALSE)
  }
  fits <- .cpf_lapply(as.list(pool_idx), function(j) {
    others <- setdiff(dn, j)
    res <- tryCatch(.fit_unit(data, j, others, fit$settings, budgets),
                    error = function(e) NULL)
    if (is.null(res) || is.null(res$weights)) return(NULL)
    syn <- as.numeric(crossprod(data$target[others, , drop = FALSE],
                                res$weights))
    data$target[j, ] - syn
  }, parallel = parallel, progress = progress,
  label = "Fitting donor placebos")
  ok <- !vapply(fits, is.null, logical(1))
  if (any(!ok)) {
    warning(sprintf("Donor placebo fits failed for: %s.",
                    paste(data$unit_names[pool_idx[!ok]], collapse = ", ")),
            call. = FALSE)
  }
  if (sum(ok) < 2L) {
    stop("Fewer than two usable donor placebo fits.", call. = FALSE)
  }
  gaps <- do.call(rbind, fits[ok])
  rownames(gaps) <- data$unit_names[pool_idx[ok]]
  colnames(gaps) <- colnames(data$target)
  list(gaps = gaps, donors = pool_idx[ok], requested = pool_idx,
       pool = pool)
}

.average_intervals <- function(x, alpha1, alpha2, gamma) {
  s <- x$settings
  psi <- x$estimate
  delta <- x$bootstrap$delta
  M1 <- c(L = stats::quantile(delta, alpha1 / 2, names = FALSE),
          U = stats::quantile(delta, 1 - alpha1 / 2, names = FALSE))
  avg <- rowMeans(x$placebo$gaps)
  centered <- avg - mean(avg)
  sigma_plugin <- stats::sd(centered)
  if (!is.finite(sigma_plugin) || sigma_plugin <= 0) {
    stop("The donor placebo gaps do not give a positive transfer scale.",
         call. = FALSE)
  }
  sigma_ls <- if (is.null(s$common_scale)) sigma_plugin else s$common_scale
  nu <- if (is.null(s$scale_upper_bound)) sigma_plugin else s$scale_upper_bound
  std <- centered / sigma_ls
  qtype <- s$empirical_quantile_type
  M2 <- rbind(
    gaussian = sigma_ls * stats::qnorm(c(alpha2 / 2, 1 - alpha2 / 2)),
    empirical = sigma_ls * c(
      stats::quantile(std, alpha2 / 2, names = FALSE, type = qtype),
      stats::quantile(std, 1 - alpha2 / 2, names = FALSE, type = qtype)),
    subgaussian = c(-1, 1) * nu * sqrt(2 * log(2 / alpha2))
  )
  level <- 1 - alpha1 - alpha2
  intervals <- data.frame(
    method = c("weight construction only", "location-scale (Gaussian)",
               "location-scale (empirical)", "sub-Gaussian bound"),
    approach = c("diagnostic", "Approach 2", "Approach 2", "Approach 1"),
    level = c(1 - alpha1, level, level, level),
    estimate = psi,
    lower = c(psi - M1[["U"]], psi - M1[["U"]] - M2[, 2] - gamma),
    upper = c(psi - M1[["L"]], psi - M1[["L"]] - M2[, 1] + gamma),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  intervals$width <- intervals$upper - intervals$lower
  list(
    intervals = intervals,
    M1 = M1,
    scale = list(
      plugin = sigma_plugin, location_scale = sigma_ls, subgaussian = nu,
      n_calibration = length(avg), placebo_mean = mean(avg),
      M2 = M2
    )
  )
}

#' @export
print.scfusion_pi <- function(x, digits = 4, ...) {
  s <- x$settings
  cat("Prediction intervals for synthetic control data fusion\n")
  cat(sprintf("  alpha1 = %s (weight construction), alpha2 = %s (transfer)",
              s$alpha1, s$alpha2))
  cat(sprintf(" -> nominal level %s%%\n",
              format(100 * (1 - s$alpha1 - s$alpha2))))
  cat(sprintf("  bootstrap: %d/%d usable refits (%s, block length %d; ",
              x$bootstrap$usable, x$bootstrap$requested, s$resample,
              s$block_len))
  cat(sprintf("budget reselection: %s)\n",
              if (s$reselect_budget) "yes" else "no"))
  cat(sprintf("  calibration donors (%s): %s\n", s$placebo_pool,
              paste(rownames(x$placebo$gaps), collapse = ", ")))
  cat(sprintf("  transfer scale: plug-in sd = %s",
              signif(x$scale$plugin, digits)))
  if (!is.null(s$common_scale)) {
    cat(sprintf(", supplied common scale = %s", s$common_scale))
  }
  if (!is.null(s$scale_upper_bound)) {
    cat(sprintf(", supplied sub-Gaussian bound = %s", s$scale_upper_bound))
  }
  cat(sprintf(" | bias allowance = %s\n\n", s$bias_allowance))
  iv <- x$intervals
  num <- c("estimate", "lower", "upper", "width")
  iv[num] <- lapply(iv[num], round, digits = digits)
  iv$level <- paste0(format(100 * iv$level), "%")
  print(iv, row.names = FALSE)
  if (x$scale$n_calibration < 10L) {
    cat(sprintf(paste0(
      "\nNote: only %d calibration donors; the empirical quantiles are ",
      "determined by the extreme placebo gaps (treat as a sensitivity ",
      "analysis).\n"), x$scale$n_calibration))
  }
  target <- if (s$bias_allowance > 0) {
    "psi_0 (given a valid bias allowance)"
  } else {
    "psi_0 + E[e | H]; read as an interval for psi_0 only if the conditional transfer bias is negligible"
  }
  cat(sprintf("\nFeasible plug-in intervals targeting %s.\n", target))
  invisible(x)
}

#' Per-period prediction bands
#'
#' Pointwise (not simultaneous) prediction bands for the treated unit's
#' untreated target-domain outcomes \eqn{Y_{1s}(0)} and for the per-period
#' effects, built from the bootstrap and placebo draws stored in an
#' `"scfusion_pi"` object. For each period \eqn{s},
#' \eqn{Y_{1s}(0) - \hat C_s = \Delta_{w,s} + e_s}; \eqn{\Delta_{w,s}} is
#' calibrated by the bootstrap weights and \eqn{e_s} by the period-\eqn{s}
#' donor placebo gaps.
#'
#' @param x An `"scfusion_pi"` object.
#' @param method Transfer-error construction: `"gaussian"` or `"empirical"`
#'   (Approach 2) or `"subgaussian"` (Approach 1).
#' @param alpha1,alpha2 Error probabilities; default to those used in `x`.
#' @param bias_allowance Non-negative scalar or per-period vector
#'   \eqn{\gamma_s}; defaults to the value used in `x`.
#' @param period_specific_scale For `"subgaussian"`: use a separate scale per
#'   period (default) or the largest one for all periods.
#' @param period_scale_upper_bound For `"subgaussian"`: optional certified
#'   per-period scale bound(s) replacing the plug-in scales.
#'
#' @return A data frame with one row per target period: the observed and
#'   synthetic outcomes, the weight-construction band (`insample_lo`,
#'   `insample_hi`), the prediction interval for \eqn{Y_{1s}(0)}, which adds
#'   the transfer error (`total_lo`, `total_hi`), the effect
#'   and its band (`effect_lo`, `effect_hi`) and the components `M1_L`,
#'   `M1_U`, `M2_L`, `M2_U`, `sigma`.
#' @examples
#' \donttest{
#' sim <- simulate_fusion_data(n_units = 16, n_ref = 12, n_target = 4,
#'                             seed = 1)
#' fit <- scfusion(sim$target, sim$reference, X = sim$X, Z = sim$Z,
#'                 budget_step = 0.1)
#' pi <- scfusion_pi(fit, B = 50, seed = 1)
#' pi_bands(pi, "empirical")
#' }
#' @export
pi_bands <- function(x, method = c("gaussian", "empirical", "subgaussian"),
                     alpha1 = NULL, alpha2 = NULL, bias_allowance = NULL,
                     period_specific_scale = TRUE,
                     period_scale_upper_bound = NULL) {
  if (!inherits(x, "scfusion_pi")) {
    stop("`x` must be an object returned by scfusion_pi().", call. = FALSE)
  }
  method <- match.arg(method)
  s <- x$settings
  if (is.null(alpha1)) alpha1 <- s$alpha1
  if (is.null(alpha2)) alpha2 <- s$alpha2
  .check_alphas(alpha1, alpha2)
  if (is.null(bias_allowance)) bias_allowance <- s$bias_allowance
  .check_flag(period_specific_scale, "period_specific_scale")
  S <- length(x$Y_treated)
  gamma <- .expand_nonneg(bias_allowance, S, "bias_allowance")

  Chat_t <- x$Chat_t
  Chat_boot <- x$bootstrap$weights %*% x$Y_control
  Delta <- matrix(Chat_t, nrow(Chat_boot), S, byrow = TRUE) - Chat_boot
  M1_L <- apply(Delta, 2L, stats::quantile, probs = alpha1 / 2,
                names = FALSE)
  M1_U <- apply(Delta, 2L, stats::quantile, probs = 1 - alpha1 / 2,
                names = FALSE)

  R <- x$placebo$gaps
  centered <- sweep(R, 2L, colMeans(R))
  sigma <- apply(centered, 2L, stats::sd)
  if (any(!is.finite(sigma))) {
    stop("Invalid per-period placebo scale.", call. = FALSE)
  }

  if (method == "subgaussian") {
    nu <- if (is.null(period_scale_upper_bound)) {
      if (period_specific_scale) sigma else rep(max(sigma), S)
    } else {
      .expand_nonneg(period_scale_upper_bound, S, "period_scale_upper_bound")
    }
    M2_U <- nu * sqrt(2 * log(2 / alpha2))
    M2_L <- -M2_U
    sigma_used <- nu
  } else {
    if (any(sigma <= 0)) {
      stop("Every period needs a positive placebo scale for the ",
           "location-scale bands.", call. = FALSE)
    }
    if (method == "gaussian") {
      M2_L <- sigma * stats::qnorm(alpha2 / 2)
      M2_U <- sigma * stats::qnorm(1 - alpha2 / 2)
    } else {
      qtype <- s$empirical_quantile_type
      M2_L <- M2_U <- numeric(S)
      for (k in seq_len(S)) {
        z <- centered[, k] / sigma[k]
        M2_L[k] <- sigma[k] * stats::quantile(z, alpha2 / 2, names = FALSE,
                                              type = qtype)
        M2_U[k] <- sigma[k] * stats::quantile(z, 1 - alpha2 / 2,
                                              names = FALSE, type = qtype)
      }
    }
    sigma_used <- sigma
  }

  total_lo <- Chat_t + M1_L + M2_L - gamma
  total_hi <- Chat_t + M1_U + M2_U + gamma
  bands <- data.frame(
    period = colnames(x$fit$data$target),
    time = seq_len(S),
    treated = as.numeric(x$Y_treated),
    synthetic = Chat_t,
    insample_lo = Chat_t + M1_L,
    insample_hi = Chat_t + M1_U,
    total_lo = total_lo,
    total_hi = total_hi,
    effect = as.numeric(x$Y_treated) - Chat_t,
    effect_lo = as.numeric(x$Y_treated) - total_hi,
    effect_hi = as.numeric(x$Y_treated) - total_lo,
    M1_L = M1_L, M1_U = M1_U, M2_L = M2_L, M2_U = M2_U,
    sigma = sigma_used,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  attr(bands, "method") <- method
  attr(bands, "level") <- 1 - alpha1 - alpha2
  attr(bands, "alpha") <- c(alpha1 = alpha1, alpha2 = alpha2)
  bands
}

.expand_nonneg <- function(x, S, name) {
  if (length(x) == 1L) x <- rep(x, S)
  if (!is.numeric(x) || length(x) != S || any(!is.finite(x)) || any(x < 0)) {
    stop(sprintf("`%s` must be a non-negative scalar or a vector with one ",
                 name), "value per target period.", call. = FALSE)
  }
  as.numeric(x)
}
