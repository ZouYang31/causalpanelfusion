#' Covariate balance table
#'
#' Compares the treated unit with its synthetic control (and with the donor
#' pool) on every covariate in both domains, the data-fusion analogue of the
#' covariate balance table reported in synthetic control studies.
#'
#' The weights are usually fitted on rescaled covariates (for example
#' min-max normalised with [normalize_minmax()]), which are hard to read.
#' Fits built from [fusion_data()] keep the original covariate values and
#' report them automatically; otherwise supply the covariates in their
#' original units through `X` and `Z`. The discrepancy on the scale the
#' balance constraints act on is always given in `gap_fit_scale`.
#'
#' @param fit An `"scfusion"` object.
#' @param X,Z Optional target- and reference-domain covariates in the units
#'   to report (one row per unit, same columns as used in `fit`). Rows are
#'   matched by row name when present, otherwise they must be in the same
#'   order as in `fit`. By default the original covariate values kept by
#'   [fusion_data()] are used if available, otherwise the covariates used
#'   for fitting.
#' @param compare Optional named list of further `"scfusion"` fits on the
#'   same units (for example with a different `eta_X`); their synthetic
#'   values are added as columns `synthetic_<name>`.
#'
#' @return A data frame of class `"scfusion_balance"` with one row per
#'   covariate and columns
#' \describe{
#'   \item{`domain`, `covariate`}{Domain and covariate name.}
#'   \item{`treated`}{Value of the treated unit.}
#'   \item{`synthetic`}{Weighted donor value \eqn{\sum_j w_j x_j} with the
#'     weights of the estimate.}
#'   \item{`donor_mean`, `donor_min`, `donor_max`}{Donor pool summaries.}
#'   \item{`gap`}{`treated - synthetic` in the reported units.}
#'   \item{`gap_fit_scale`}{`treated - synthetic` on the scale used for
#'     fitting.}
#'   \item{`treated_rank`}{Rank of the treated unit among all units
#'     (1 = largest). A rank near 1 or near the number of units means the
#'     treated unit lies at the edge of the donor pool and cannot be matched
#'     closely by a convex combination.}
#' }
#' The fitted discrepancies \eqn{NSE_Z}, \eqn{NSE_X}, their best attainable
#' values and the limits allowed by `eta_Z`, `eta_X` are stored as
#' attributes and shown by `print()`.
#'
#' @details Weights below the threshold of [scfusion()] are set to zero
#'   without rescaling, so they can sum to slightly less than one; the
#'   synthetic values use them exactly as the estimate does.
#'
#' @examples
#' sim <- simulate_fusion_data(n_units = 16, n_ref = 12, n_target = 4,
#'                             seed = 1)
#' # covariates in "original units" (rescaled here for illustration) are
#' # min-max normalised for fitting and passed back for the table
#' X_raw <- sim$X * 1000
#' Z_raw <- sim$Z * 50
#' fit <- scfusion(sim$target, sim$reference,
#'                 X = normalize_minmax(X_raw), Z = normalize_minmax(Z_raw),
#'                 budget_step = 0.1)
#' covariate_balance(fit, X = X_raw, Z = Z_raw)
#'
#' # side by side with a fit that enforces tighter target-domain balance
#' fit_tight <- scfusion(sim$target, sim$reference,
#'                       X = normalize_minmax(X_raw),
#'                       Z = normalize_minmax(Z_raw),
#'                       eta_X = 0.01, budget_step = 0.1)
#' covariate_balance(fit, X = X_raw, Z = Z_raw,
#'                   compare = list(eta_X_0.01 = fit_tight))
#' @export
covariate_balance <- function(fit, X = NULL, Z = NULL, compare = NULL) {
  if (!inherits(fit, "scfusion")) {
    stop("`fit` must be an object returned by scfusion().", call. = FALSE)
  }
  data <- fit$data
  compare <- .check_compare(compare, fit)
  default_X <- if (is.null(data$X_raw)) data$X else data$X_raw
  default_Z <- if (is.null(data$Z_raw)) data$Z else data$Z_raw
  Zr <- .reporting_covariates(Z, default_Z, data$unit_names, "Z")
  Xr <- .reporting_covariates(X, default_X, data$unit_names, "X")

  one <- function(raw, fitm, domain) {
    tr <- fit$treated
    dn <- fit$donors
    w <- as.numeric(fit$weights)
    donors <- raw[dn, , drop = FALSE]
    synthetic <- as.numeric(crossprod(donors, w))
    out <- data.frame(
      domain = domain,
      covariate = colnames(raw),
      treated = as.numeric(raw[tr, ]),
      synthetic = synthetic,
      donor_mean = colMeans(donors),
      donor_min = apply(donors, 2L, min),
      donor_max = apply(donors, 2L, max),
      gap = as.numeric(raw[tr, ]) - synthetic,
      gap_fit_scale = as.numeric(fitm[tr, ]) -
        as.numeric(crossprod(fitm[dn, , drop = FALSE], w)),
      treated_rank = vapply(seq_len(ncol(raw)), function(k) {
        as.integer(rank(-raw[, k], ties.method = "min")[tr])
      }, integer(1)),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    for (nm in names(compare)) {
      wc <- as.numeric(compare[[nm]]$weights)
      out[[paste0("synthetic_", nm)]] <- as.numeric(crossprod(donors, wc))
    }
    out
  }
  dn <- .domain_names(data)
  out <- rbind(one(Zr, data$Z, dn[["reference"]]),
               one(Xr, data$X, dn[["target"]]))

  s <- fit$settings
  allowed <- c(Z = (1 + s$eta_Z) * (1 + fit$baseline[["Z"]]) - 1,
               X = (1 + s$eta_X) * (1 + fit$baseline[["X"]]) - 1)
  attr(out, "nse") <- fit$nse[c("Z", "X")]
  attr(out, "baseline") <- fit$baseline
  attr(out, "allowed") <- allowed
  attr(out, "weight_sum") <- sum(fit$weights)
  attr(out, "estimates") <- c(main = fit$estimate,
                              vapply(compare, `[[`, numeric(1), "estimate"))
  attr(out, "treated") <- data$unit_names[fit$treated]
  attr(out, "n_units") <- length(data$unit_names)
  attr(out, "original_units") <- !is.null(X) || !is.null(Z) ||
    isTRUE(data$normalized)
  class(out) <- c("scfusion_balance", "data.frame")
  out
}

.check_compare <- function(compare, fit) {
  if (is.null(compare)) return(list())
  if (inherits(compare, "scfusion")) compare <- list(compare)
  if (!is.list(compare) ||
      !all(vapply(compare, inherits, logical(1), "scfusion"))) {
    stop("`compare` must be a list of scfusion() fits.", call. = FALSE)
  }
  nms <- names(compare)
  if (is.null(nms)) nms <- rep("", length(compare))
  nms[nms == ""] <- paste0("fit", which(nms == "") + 1L)
  names(compare) <- make.names(nms, unique = TRUE)
  for (f in compare) {
    if (!identical(f$data$unit_names, fit$data$unit_names) ||
        f$treated != fit$treated) {
      stop("Fits in `compare` must use the same units and treated unit as ",
           "`fit`.", call. = FALSE)
    }
  }
  compare
}

# Align user-supplied reporting covariates with the fitted ones.
.reporting_covariates <- function(raw, fitm, units, name) {
  if (is.null(raw)) return(fitm)
  m <- .as_numeric_matrix(raw, name)
  if (!is.null(rownames(m))) {
    if (length(rownames(m)) != length(units) ||
        !setequal(rownames(m), units)) {
      stop(sprintf("Row names of `%s` do not match the units of `fit`.",
                   name), call. = FALSE)
    }
    m <- m[units, , drop = FALSE]
  } else if (nrow(m) != length(units)) {
    stop(sprintf("`%s` must have one row per unit (%d rows).", name,
                 length(units)), call. = FALSE)
  }
  if (ncol(m) != ncol(fitm)) {
    stop(sprintf("`%s` must have the same %d columns as used in `fit`.",
                 name, ncol(fitm)), call. = FALSE)
  }
  if (is.null(colnames(m))) {
    colnames(m) <- colnames(fitm)
  } else if (setequal(colnames(m), colnames(fitm))) {
    m <- m[, colnames(fitm), drop = FALSE]
  }
  rownames(m) <- units
  m
}

#' @rdname covariate_balance
#' @param x An `"scfusion_balance"` object.
#' @param digits Significant digits to print.
#' @param ... Unused.
#' @export
print.scfusion_balance <- function(x, digits = 4, ...) {
  if (is.null(attr(x, "nse"))) return(NextMethod())   # e.g. after subsetting
  cat(sprintf("Covariate balance: %s vs its synthetic control (%s units)\n\n",
              attr(x, "treated"),
              if (isTRUE(attr(x, "original_units"))) "reported" else
                "fitting-scale"))
  tab <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  attributes(tab)[setdiff(names(attributes(tab)),
                          c("names", "row.names", "class"))] <- NULL
  if (!isTRUE(attr(x, "original_units"))) tab$gap_fit_scale <- NULL  # = gap
  num <- vapply(tab, is.numeric, logical(1)) & names(tab) != "treated_rank"
  tab[num] <- lapply(tab[num], function(v) {
    formatC(v, digits = digits, format = "fg", big.mark = ",")
  })
  print(tab, row.names = FALSE, right = TRUE)
  nse <- attr(x, "nse")
  base <- attr(x, "baseline")
  allowed <- attr(x, "allowed")
  cat(sprintf("\ntreated_rank: 1 = largest of %d units.\n", attr(x, "n_units")))
  cat("Fitting-scale discrepancy NSE (fitted | best attainable | allowed):\n")
  cat(sprintf("  reference: %s | %s | %s\n", signif(nse[["Z"]], 3),
              signif(base[["Z"]], 3), signif(allowed[["Z"]], 3)))
  cat(sprintf("  target:    %s | %s | %s\n", signif(nse[["X"]], 3),
              signif(base[["X"]], 3), signif(allowed[["X"]], 3)))
  ws <- attr(x, "weight_sum")
  if (abs(ws - 1) > 1e-6) {
    cat(sprintf(paste0("Weights sum to %s (weights below the threshold are ",
                       "set to zero without rescaling).\n"),
                format(round(ws, 4))))
  }
  est <- attr(x, "estimates")
  if (length(est) > 1L) {
    cat("Estimates: ", paste(sprintf("%s = %s", names(est),
                                     .fmt(est, digits)), collapse = " | "),
        "\n", sep = "")
  }
  invisible(x)
}
