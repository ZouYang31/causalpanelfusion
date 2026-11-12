#' Sensitivity of the estimate to the balance slack \eqn{\eta}
#'
#' Refits the weights over a grid of `(eta_Z, eta_X)` values and reports the
#' estimated effect, the fit criteria and the covariate balance for each
#' combination. The baselines \eqn{NSE_Z^*}, \eqn{NSE_X^*} do not depend on
#' \eqn{\eta} and are reused from `fit`.
#'
#' @param fit An `"scfusion"` object.
#' @param eta_Z,eta_X Numeric vectors of slack values; all combinations are
#'   fitted.
#' @param parallel,progress As in [scfusion_pi()].
#'
#' @return A data frame of class `"scfusion_eta"` with one row per
#'   combination: `eta_Z`, `eta_X`, `estimate`, `NSE_F`, `NSE_Z`, `NSE_X`,
#'   the selected budget, the number of positive weights and the absolute
#'   balance gap for every covariate (`absdiff_Z_<name>`,
#'   `absdiff_X_<name>`). Rows with no feasible fit have `NA` values. Use
#'   `plot()` for a heat map.
#' @examples
#' \donttest{
#' sim <- simulate_fusion_data(n_units = 12, n_ref = 10, n_target = 4,
#'                             seed = 2)
#' fit <- scfusion(sim$target, sim$reference, X = sim$X, Z = sim$Z,
#'                 budget_step = 0.1)
#' sens <- eta_sensitivity(fit, eta_Z = c(0.05, 0.2), eta_X = c(0.05, 0.2))
#' sens
#' plot(sens)
#' }
#' @export
eta_sensitivity <- function(fit, eta_Z = c(0.05, 0.1, 0.15, 0.2),
                            eta_X = c(0.05, 0.1, 0.15, 0.2),
                            parallel = FALSE, progress = interactive()) {
  if (!inherits(fit, "scfusion")) {
    stop("`fit` must be an object returned by scfusion().", call. = FALSE)
  }
  if (!is.numeric(eta_Z) || !length(eta_Z) || any(!is.finite(eta_Z)) ||
      any(eta_Z < 0) || !is.numeric(eta_X) || !length(eta_X) ||
      any(!is.finite(eta_X)) || any(eta_X < 0)) {
    stop("`eta_Z` and `eta_X` must be non-negative numeric vectors.",
         call. = FALSE)
  }
  .check_flag(parallel, "parallel")
  .check_flag(progress, "progress")
  if (fit$settings$solver == "gurobi") .require_gurobi()

  data <- fit$data
  tr <- fit$treated
  dn <- fit$donors
  grid <- expand.grid(eta_X = eta_X, eta_Z = eta_Z)[, c("eta_Z", "eta_X")]
  prob <- .fusion_problem(
    data$reference[tr, ], data$reference[dn, , drop = FALSE],
    data$Z[tr, ], data$Z[dn, , drop = FALSE],
    data$X[tr, ], data$X[dn, , drop = FALSE]
  )
  s <- fit$settings
  rows <- .cpf_lapply(seq_len(nrow(grid)), function(k) {
    res <- tryCatch(
      .fit_fusion_weights(prob, s$budgets, grid$eta_Z[k], grid$eta_X[k],
                          s$solver, s$weight_threshold,
                          baseline = fit$baseline),
      error = function(e) NULL)
    if (is.null(res) || is.null(res$weights)) return(NULL)
    res
  }, parallel = parallel, progress = progress,
  label = "Refitting over the eta grid")

  zn <- paste0("absdiff_Z_", colnames(data$Z))
  xn <- paste0("absdiff_X_", colnames(data$X))
  out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(k) {
    res <- rows[[k]]
    if (is.null(res)) {
      vals <- c(NA_real_, rep(NA_real_, 3L), rep(NA_real_, 3L), NA_real_,
                rep(NA_real_, length(zn) + length(xn)))
    } else {
      w <- res$weights
      syn_Y <- as.numeric(crossprod(data$target[dn, , drop = FALSE], w))
      vals <- c(
        mean(data$target[tr, ] - syn_Y),
        .nse(w, prob$F1, prob$F0), .nse(w, prob$Z1, prob$Z0),
        .nse(w, prob$X1, prob$X0),
        res$budget, sum(w > 0),
        abs(prob$Z1 - as.numeric(crossprod(prob$Z0, w))),
        abs(prob$X1 - as.numeric(crossprod(prob$X0, w)))
      )
    }
    vals
  }))
  out <- as.data.frame(out)
  names(out) <- c("estimate", "NSE_F", "NSE_Z", "NSE_X", "b_F", "b_Z", "b_X",
                  "n_donors", zn, xn)
  out <- cbind(grid, out)
  if (any(is.na(out$estimate))) {
    warning("No feasible fit for some (eta_Z, eta_X) combinations.",
            call. = FALSE)
  }
  class(out) <- c("scfusion_eta", "data.frame")
  out
}

#' Sensitivity of prediction intervals to the transfer-error scale
#'
#' The transfer error is not identified nonparametrically, so this traces
#' the prediction interval as the transfer scale is multiplied by
#' \eqn{\rho}: \eqn{\rho = 0} keeps only the weight-construction error,
#' \eqn{\rho = 1} is the estimated scale, and \eqn{\rho > 1} stress-tests
#' stronger target-only factors or larger target noise.
#'
#' @param x An `"scfusion_pi"` object.
#' @param rho Non-negative multipliers of the transfer scale.
#' @param method `"gaussian"`, `"empirical"` or `"subgaussian"`.
#' @param period `NULL` (default) for the time-averaged effect, or the index
#'   of a target period for that period's effect.
#'
#' @return A data frame with columns `rho`, `scale`, `estimate`, `lower`,
#'   `upper` (class `"scfusion_sensitivity"`, with a `plot()` method).
#' @examples
#' \donttest{
#' sim <- simulate_fusion_data(n_units = 16, n_ref = 12, n_target = 4,
#'                             seed = 1)
#' fit <- scfusion(sim$target, sim$reference, X = sim$X, Z = sim$Z,
#'                 budget_step = 0.1)
#' pi <- scfusion_pi(fit, B = 50, seed = 1)
#' sensitivity_pi(pi)
#' }
#' @export
sensitivity_pi <- function(x, rho = c(0, 0.5, 1, 1.5, 2, 3),
                           method = c("gaussian", "empirical", "subgaussian"),
                           period = NULL) {
  if (!inherits(x, "scfusion_pi")) {
    stop("`x` must be an object returned by scfusion_pi().", call. = FALSE)
  }
  method <- match.arg(method)
  if (!is.numeric(rho) || !length(rho) || any(!is.finite(rho)) ||
      any(rho < 0)) {
    stop("`rho` must be a vector of non-negative numbers.", call. = FALSE)
  }
  s <- x$settings
  gamma <- s$bias_allowance
  if (is.null(period)) {
    est <- x$estimate
    M1 <- x$M1
    M2 <- x$scale$M2[method, ]
    scale <- switch(method, subgaussian = x$scale$subgaussian,
                    x$scale$location_scale)
  } else {
    S <- length(x$Y_treated)
    .check_scalar(period, "period", lower = 1, upper = S, integer = TRUE)
    b <- pi_bands(x, method = method)[period, ]
    est <- b$effect
    M1 <- c(L = b$M1_L, U = b$M1_U)
    M2 <- c(b$M2_L, b$M2_U)
    scale <- b$sigma
    gamma <- .expand_nonneg(gamma, S, "bias_allowance")[period]
  }
  out <- data.frame(
    rho = rho,
    scale = rho * scale,
    estimate = est,
    lower = est - M1[["U"]] - rho * M2[2] - gamma,
    upper = est - M1[["L"]] - rho * M2[1] + gamma
  )
  attr(out, "period") <- period
  attr(out, "method") <- method
  attr(out, "level") <- 1 - s$alpha1 - s$alpha2
  class(out) <- c("scfusion_sensitivity", "data.frame")
  out
}
