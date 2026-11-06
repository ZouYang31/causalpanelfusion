#' Synthetic control data fusion
#'
#' Estimates the effect of an intervention on a single treated unit when the
#' target domain has **no pre-intervention period**, by borrowing a reference
#' domain observed for the same units (Yang, Lee, Kohler and Ghassami). Donor
#' weights are learned in the reference domain and transported to the target
#' domain.
#'
#' @section Algorithm:
#' Write \eqn{F} for reference-domain outcomes, \eqn{Z} and \eqn{X} for the
#' reference- and target-domain covariates, and \eqn{NSE_F(w)},
#' \eqn{NSE_Z(w)}, \eqn{NSE_X(w)} for the normalised squared discrepancies
#' between the treated unit and the weighted donors (squared error divided by
#' the number of columns). The weights are chosen by the data-driven
#' constrained optimisation of the paper:
#' \enumerate{
#'   \item Baselines \eqn{NSE_Z^*} and \eqn{NSE_X^*}: the smallest covariate
#'     discrepancies attainable with simplex weights.
#'   \item For every budget \eqn{b = (b_F, b_Z, b_X)} in [budget_grid()],
#'     minimise \eqn{b_F NSE_F(w) + b_Z NSE_Z(w) + b_X NSE_X(w)} over the
#'     simplex subject to \eqn{NSE_Z(w) \le (1 + \eta_Z)(1 + NSE_Z^*) - 1}
#'     and \eqn{NSE_X(w) \le (1 + \eta_X)(1 + NSE_X^*) - 1}.
#'   \item Keep the budget whose weights give the smallest \eqn{NSE_F}.
#' }
#' Weights below `weight_threshold` are set to zero (without renormalising),
#' as in the original implementation.
#'
#' The estimate of the time-averaged effect on the treated unit is
#' \eqn{\hat\psi = \bar Y_1 - \sum_j w_j \bar Y_j}, and the per-period gaps
#' are \eqn{Y_{1s} - \sum_j w_j Y_{js}}.
#'
#' @section Solvers:
#' `solver = "ecos"` (default) solves every problem as a second-order cone
#' program with the free ECOS solver (package \pkg{ECOSolveR}). `solver =
#' "gurobi"` reproduces the original research implementation: baselines from
#' `kernlab::ipop()` and the constrained problems with Gurobi. It needs the
#' \pkg{gurobi} R package, a Gurobi licence and \pkg{kernlab}. Both solve the
#' same convex problems, so estimates agree up to solver tolerance.
#'
#' @param target Target-domain outcomes: numeric matrix (or data frame) with
#'   one row per unit and one column per post-intervention period (the
#'   paper's \eqn{Y}). The treated unit's row holds its observed (treated)
#'   outcomes.
#' @param reference Reference-domain outcomes, units by periods (the paper's
#'   \eqn{F}). Unaffected by the intervention for every unit.
#' @param X Target-domain covariates, units by covariates.
#' @param Z Reference-domain covariates, units by covariates.
#' @param treated The treated unit: a row index or a row name.
#' @param eta_Z,eta_X Non-negative slack for the covariate-balance constraints
#'   (see Algorithm). Defaults `0.1`.
#' @param budget_step,budget_min Spacing and lower bound of the budget grid;
#'   see [budget_grid()]. Ignored when `budgets` is supplied.
#' @param budgets Optional matrix of candidate budgets with columns `b_F`,
#'   `b_Z`, `b_X` (rows summing to one), replacing the default grid.
#' @param solver `"ecos"` (default) or `"gurobi"`; see Solvers.
#' @param weight_threshold Weights below this value are set to zero.
#'   Default `0.005`.
#' @param keep_path If `TRUE`, store \eqn{NSE_F} for every candidate budget.
#'
#' @details If all four inputs have row names, rows are matched by name (the
#'   order of `target` is used); otherwise rows must already be aligned.
#'   Covariates should be on comparable scales (the paper min-max normalises
#'   them; see [normalize_minmax()]).
#'
#' @return An object of class `"scfusion"`, a list with elements
#' \describe{
#'   \item{`estimate`}{Time-averaged effect \eqn{\hat\psi}.}
#'   \item{`effects`}{Per-period gaps in the target domain.}
#'   \item{`weights`}{Named donor weights.}
#'   \item{`budget`}{Selected budget \eqn{(b_F, b_Z, b_X)}.}
#'   \item{`nse`}{\eqn{NSE_F}, \eqn{NSE_Z}, \eqn{NSE_X} at the weights.}
#'   \item{`baseline`}{The baselines \eqn{NSE_Z^*} and \eqn{NSE_X^*}.}
#'   \item{`synthetic_target`, `synthetic_reference`}{Synthetic trajectories.}
#'   \item{`balance`}{Covariate balance table.}
#'   \item{`data`, `treated`, `donors`, `settings`, `call`}{Inputs and options
#'     used by [scfusion_pi()], [placebo_test()] and the plot methods.}
#' }
#'
#' @references Yang, Z., Lee, S. H., Kohler, J. R. and Ghassami, A. Causal
#'   data fusion for panel data without a pre-intervention period.
#'   arXiv:2410.16391.
#'
#' @seealso [scfusion_pi()] for prediction intervals, [placebo_test()],
#'   [eta_sensitivity()], [equi_confounding()].
#'
#' @examples
#' sim <- simulate_fusion_data(n_units = 16, n_ref = 12, n_target = 4,
#'                             seed = 1)
#' # a coarse budget grid keeps the example fast
#' fit <- scfusion(sim$target, sim$reference, X = sim$X, Z = sim$Z,
#'                 budget_step = 0.05)
#' fit
#' summary(fit)
#' sim$effect            # true per-period effects
#' plot(fit)
#' @export
scfusion <- function(target, reference, X, Z, treated = 1L,
                     eta_Z = 0.1, eta_X = 0.1,
                     budget_step = 0.01, budget_min = 0.01, budgets = NULL,
                     solver = c("ecos", "gurobi"),
                     weight_threshold = 0.005, keep_path = FALSE) {
  cl <- match.call()
  solver <- match.arg(solver)
  .check_scalar(eta_Z, "eta_Z", lower = 0)
  .check_scalar(eta_X, "eta_X", lower = 0)
  .check_scalar(weight_threshold, "weight_threshold", lower = 0, upper = 1)
  .check_flag(keep_path, "keep_path")
  if (solver == "gurobi") .require_gurobi()

  data <- .prepare_fusion_data(target, reference, X, Z, treated)
  budgets <- .resolve_budgets(budgets, budget_step, budget_min)
  settings <- list(
    eta_Z = eta_Z, eta_X = eta_X, solver = solver,
    weight_threshold = weight_threshold,
    budget_step = if (is.null(attr(budgets, "custom"))) budget_step else NA,
    budget_min = if (is.null(attr(budgets, "custom"))) budget_min else NA,
    budgets = budgets
  )

  res <- .fit_unit(data, data$treated, data$donors, settings, budgets,
                   keep_path = keep_path)
  if (is.null(res$weights)) {
    stop("No candidate budget produced a feasible solution. Try larger ",
         "`eta_Z`/`eta_X` or check the covariates.", call. = FALSE)
  }
  out <- .assemble_fit(data, settings, res)
  out$call <- cl
  if (keep_path) out$path <- res$path
  out
}

# Build the object returned by scfusion() from a weight fit.
.assemble_fit <- function(data, settings, res) {
  tr <- data$treated
  dn <- data$donors
  w <- res$weights
  names(w) <- data$unit_names[dn]
  Y <- data$target
  Fm <- data$reference
  syn_Y <- as.numeric(crossprod(Y[dn, , drop = FALSE], w))
  syn_F <- as.numeric(crossprod(Fm[dn, , drop = FALSE], w))
  names(syn_Y) <- colnames(Y)
  names(syn_F) <- colnames(Fm)
  effects <- Y[tr, ] - syn_Y
  nse <- c(
    F = .nse(w, Fm[tr, ], Fm[dn, , drop = FALSE]),
    Z = .nse(w, data$Z[tr, ], data$Z[dn, , drop = FALSE]),
    X = .nse(w, data$X[tr, ], data$X[dn, , drop = FALSE])
  )
  structure(
    list(
      estimate = mean(effects),
      effects = effects,
      weights = w,
      budget = res$budget,
      nse = nse,
      baseline = res$baseline,
      synthetic_target = syn_Y,
      synthetic_reference = syn_F,
      balance = .balance_table(data, w),
      n_feasible = res$n_feasible,
      data = data,
      treated = tr,
      donors = dn,
      settings = settings
    ),
    class = "scfusion"
  )
}

# Fit weights for `target_idx` using donors `donor_idx`.
# `F_treated` optionally replaces the target unit's reference outcomes
# (used by the bootstrap).
.fit_unit <- function(data, target_idx, donor_idx, settings, budgets,
                      F_treated = NULL, keep_path = FALSE) {
  Fm <- data$reference
  prob <- .fusion_problem(
    F_treated = if (is.null(F_treated)) Fm[target_idx, ] else F_treated,
    F_control = Fm[donor_idx, , drop = FALSE],
    Z_treated = data$Z[target_idx, ],
    Z_control = data$Z[donor_idx, , drop = FALSE],
    X_treated = data$X[target_idx, ],
    X_control = data$X[donor_idx, , drop = FALSE]
  )
  .fit_fusion_weights(prob, budgets, settings$eta_Z, settings$eta_X,
                      settings$solver, settings$weight_threshold,
                      keep_path = keep_path)
}

.resolve_budgets <- function(budgets, budget_step, budget_min) {
  if (is.null(budgets)) return(budget_grid(budget_step, budget_min))
  if (is.data.frame(budgets)) budgets <- as.matrix(budgets)
  if (is.list(budgets) && !is.matrix(budgets)) {
    budgets <- do.call(rbind, lapply(budgets, function(b) {
      unlist(b)[c("b_F", "b_Z", "b_X")]
    }))
  }
  if (!is.matrix(budgets) || !is.numeric(budgets) || ncol(budgets) != 3L) {
    stop("`budgets` must be a numeric matrix with three columns ",
         "(b_F, b_Z, b_X).", call. = FALSE)
  }
  if (!is.null(colnames(budgets))) {
    if (!all(c("b_F", "b_Z", "b_X") %in% colnames(budgets))) {
      stop("`budgets` columns must be named b_F, b_Z and b_X.", call. = FALSE)
    }
    budgets <- budgets[, c("b_F", "b_Z", "b_X"), drop = FALSE]
  }
  colnames(budgets) <- c("b_F", "b_Z", "b_X")
  if (any(!is.finite(budgets)) || any(budgets < 0) ||
      any(abs(rowSums(budgets) - 1) > 1e-6)) {
    stop("Each budget must be non-negative and sum to one.", call. = FALSE)
  }
  attr(budgets, "custom") <- TRUE
  budgets
}

# Validate and align the four input blocks.
.prepare_fusion_data <- function(target, reference, X, Z, treated) {
  Y <- .as_numeric_matrix(target, "target")
  Fm <- .as_numeric_matrix(reference, "reference")
  X <- .as_numeric_matrix(X, "X")
  Z <- .as_numeric_matrix(Z, "Z")
  rn <- list(rownames(Y), rownames(Fm), rownames(X), rownames(Z))
  if (all(!vapply(rn, is.null, logical(1)))) {
    units <- rn[[1L]]
    if (anyDuplicated(units)) {
      stop("Row names of `target` must be unique.", call. = FALSE)
    }
    for (k in 2:4) {
      if (!setequal(units, rn[[k]]) || length(rn[[k]]) != length(units)) {
        stop("`target`, `reference`, `X` and `Z` have row names that do not ",
             "match.", call. = FALSE)
      }
    }
    Fm <- Fm[units, , drop = FALSE]
    X <- X[units, , drop = FALSE]
    Z <- Z[units, , drop = FALSE]
  } else {
    n <- c(nrow(Y), nrow(Fm), nrow(X), nrow(Z))
    if (length(unique(n)) != 1L) {
      stop("`target`, `reference`, `X` and `Z` must have the same number of ",
           "rows (one per unit).", call. = FALSE)
    }
    units <- rn[!vapply(rn, is.null, logical(1))]
    units <- if (length(units)) units[[1L]] else paste0("unit_", seq_len(n[1]))
  }
  N <- nrow(Y)
  if (N < 3L) stop("At least three units are required.", call. = FALSE)
  if (is.character(treated)) {
    idx <- match(treated, units)
    if (length(treated) != 1L || is.na(idx)) {
      stop("`treated` must be one of the unit (row) names.", call. = FALSE)
    }
    treated <- idx
  }
  .check_scalar(treated, "treated", lower = 1, upper = N, integer = TRUE)
  treated <- as.integer(treated)
  rownames(Y) <- rownames(Fm) <- rownames(X) <- rownames(Z) <- units
  if (is.null(colnames(Y))) colnames(Y) <- paste0("s", seq_len(ncol(Y)))
  if (is.null(colnames(Fm))) colnames(Fm) <- paste0("t", seq_len(ncol(Fm)))
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  if (is.null(colnames(Z))) colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  list(target = Y, reference = Fm, X = X, Z = Z, treated = treated,
       donors = setdiff(seq_len(N), treated), unit_names = units)
}

.balance_table <- function(data, w) {
  one <- function(M, domain) {
    tr <- M[data$treated, ]
    dn <- M[data$donors, , drop = FALSE]
    data.frame(
      domain = domain,
      covariate = colnames(M),
      treated = as.numeric(tr),
      synthetic = as.numeric(crossprod(dn, w)),
      donor_mean = colMeans(dn),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }
  rbind(one(data$Z, "reference (Z)"), one(data$X, "target (X)"))
}

#' @export
print.scfusion <- function(x, digits = 4, ...) {
  cat("Synthetic control data fusion\n")
  cat(sprintf("  treated unit: %s | donors: %d | target periods: %d | ",
              x$data$unit_names[x$treated], length(x$donors),
              ncol(x$data$target)))
  cat(sprintf("reference periods: %d\n", ncol(x$data$reference)))
  cat(sprintf("  estimated average effect (psi_hat): %s\n",
              .fmt(x$estimate, digits)))
  pos <- x$weights[x$weights > 0]
  pos <- sort(pos, decreasing = TRUE)
  cat(sprintf("  donors with positive weight: %d\n", length(pos)))
  print(round(pos, digits))
  invisible(x)
}

#' Summarise a synthetic control data-fusion fit
#'
#' @param object An `"scfusion"` object.
#' @param x A `"summary.scfusion"` object.
#' @param digits Number of digits to print.
#' @param ... Unused.
#' @return `summary()` returns an object of class `"summary.scfusion"`.
#' @export
summary.scfusion <- function(object, ...) {
  x <- object
  structure(
    list(
      estimate = x$estimate,
      effects = data.frame(
        period = names(x$effects),
        treated = as.numeric(x$data$target[x$treated, ]),
        synthetic = as.numeric(x$synthetic_target),
        effect = as.numeric(x$effects),
        row.names = NULL
      ),
      weights = x$weights[x$weights > 0],
      budget = x$budget,
      nse = x$nse,
      baseline = x$baseline,
      balance = x$balance,
      settings = x$settings[c("eta_Z", "eta_X", "solver",
                              "weight_threshold")],
      n_budgets = nrow(x$settings$budgets),
      n_feasible = x$n_feasible,
      treated = x$data$unit_names[x$treated]
    ),
    class = "summary.scfusion"
  )
}

#' @rdname summary.scfusion
#' @export
print.summary.scfusion <- function(x, digits = 4, ...) {
  cat("Synthetic control data fusion\n\n")
  cat(sprintf("Treated unit: %s\n", x$treated))
  cat(sprintf("Estimated average effect (psi_hat): %s\n\n",
              .fmt(x$estimate, digits)))
  cat("Per-period effects (target domain):\n")
  print(x$effects, digits = digits, row.names = FALSE)
  cat("\nDonor weights (positive only):\n")
  print(round(sort(x$weights, decreasing = TRUE), digits))
  cat(sprintf(
    "\nSelected budget: b_F = %s, b_Z = %s, b_X = %s (%d of %d candidates feasible)\n",
    format(x$budget[["b_F"]]), format(x$budget[["b_Z"]]),
    format(x$budget[["b_X"]]), x$n_feasible, x$n_budgets))
  cat(sprintf("Fit: NSE_F = %s | NSE_Z = %s (baseline %s) | NSE_X = %s (baseline %s)\n",
              signif(x$nse[["F"]], digits), signif(x$nse[["Z"]], digits),
              signif(x$baseline[["Z"]], digits), signif(x$nse[["X"]], digits),
              signif(x$baseline[["X"]], digits)))
  cat(sprintf("Settings: eta_Z = %s, eta_X = %s, solver = %s\n\n",
              x$settings$eta_Z, x$settings$eta_X, x$settings$solver))
  cat("Covariate balance:\n")
  print(x$balance, digits = digits, row.names = FALSE)
  invisible(x)
}

#' @export
coef.scfusion <- function(object, ...) object$weights
