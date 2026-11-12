#' Equi-confounding data fusion
#'
#' Estimates the average effect on the treated unit under the
#' equi-confounding assumptions of the paper, which anchor the target
#' domain's unobserved confounding bias at the bias observed in the reference
#' domain:
#' \describe{
#'   \item{Linear}{\eqn{\hat\psi = (\bar Y_1 - \bar Y_0) - (\bar F_1 - \bar
#'     F_0)}, a difference-in-differences across domains.}
#'   \item{Logarithmic}{\eqn{\hat\psi = \bar Y_1 - (\bar F_1 / \bar F_0)\,\bar
#'     Y_0}, i.e. equal ratios across domains (average outcomes must be positive).}
#' }
#' Here \eqn{\bar Y_1}, \eqn{\bar F_1} average the treated unit's target and
#' reference outcomes over time and \eqn{\bar Y_0}, \eqn{\bar F_0} average
#' all donors over units and time.
#'
#' @param target Target-domain outcomes (units by periods), or a
#'   [fusion_data()] object (then `reference` and `treated` are taken from
#'   it).
#' @param reference Reference-domain outcomes (units by periods).
#' @param treated The treated unit: a row index or a row name.
#' @param type Which estimator(s) to compute: `"both"` (default),
#'   `"linear"` or `"log"`.
#'
#' @return An object of class `"equi_confounding"`: a list with the
#'   estimate(s) and the four means used, with a `plot()` method that draws
#'   the equi-confounding diagram.
#' @examples
#' sim <- simulate_fusion_data(seed = 3)
#' ec <- equi_confounding(sim$target, sim$reference)
#' ec
#' plot(ec)
#' @export
equi_confounding <- function(target, reference, treated = 1L,
                             type = c("both", "linear", "log")) {
  type <- match.arg(type)
  if (inherits(target, "fusion_data")) {
    if (!missing(reference) || !missing(treated)) {
      stop("When `target` is a fusion_data() object, `reference` and ",
           "`treated` are taken from it; do not supply them.", call. = FALSE)
    }
    reference <- target$reference
    treated <- target$treated
    target <- target$target
  }
  Y <- .as_numeric_matrix(target, "target")
  Fm <- .as_numeric_matrix(reference, "reference")
  if (nrow(Y) != nrow(Fm)) {
    stop("`target` and `reference` must have the same rows (units).",
         call. = FALSE)
  }
  if (!is.null(rownames(Y)) && !is.null(rownames(Fm))) {
    if (!setequal(rownames(Y), rownames(Fm))) {
      stop("Row names of `target` and `reference` do not match.",
           call. = FALSE)
    }
    Fm <- Fm[rownames(Y), , drop = FALSE]
  }
  units <- rownames(Y)
  if (is.character(treated)) {
    idx <- match(treated, units)
    if (length(treated) != 1L || is.na(idx)) {
      stop("`treated` must be one of the unit (row) names.", call. = FALSE)
    }
    treated <- idx
  }
  .check_scalar(treated, "treated", lower = 1, upper = nrow(Y), integer = TRUE)
  means <- c(
    Y1 = mean(Y[treated, ]), Y0 = mean(Y[-treated, , drop = FALSE]),
    F1 = mean(Fm[treated, ]), F0 = mean(Fm[-treated, , drop = FALSE])
  )
  est <- c(
    linear = (means[["Y1"]] - means[["Y0"]]) - (means[["F1"]] - means[["F0"]]),
    log = means[["Y1"]] - (means[["F1"]] / means[["F0"]]) * means[["Y0"]]
  )
  if ((type == "log" || type == "both") && any(means <= 0)) {
    warning("The logarithmic estimator needs positive average outcomes.",
            call. = FALSE)
  }
  if (type != "both") est <- est[type]
  structure(list(estimate = est, means = means,
                 treated = if (is.null(units)) treated else units[treated]),
            class = "equi_confounding")
}

#' @export
print.equi_confounding <- function(x, digits = 4, ...) {
  cat("Equi-confounding data fusion\n")
  for (nm in names(x$estimate)) {
    cat(sprintf("  %-11s psi_hat = %s\n",
                if (nm == "linear") "linear:" else "logarithmic:",
                .fmt(x$estimate[[nm]], digits)))
  }
  invisible(x)
}
