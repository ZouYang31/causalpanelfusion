#' Placebo tests for synthetic control data fusion
#'
#' Two diagnostics that re-run the full fitting algorithm of [scfusion()]:
#' \describe{
#'   \item{`"in_space"`}{Each donor in turn is treated as if it were the
#'     treated unit and fitted from the remaining units. The treated unit's
#'     gaps should stand out from the placebo gaps in the target domain while
#'     all units track their synthetic controls in the reference domain. The
#'     result also carries the RMSPE-ratio table used by [rmspe_test()].}
#'   \item{`"leave_one_out"`}{Each donor with positive weight is dropped in
#'     turn and the treated unit is refitted, showing how much the estimate
#'     depends on any single donor.}
#' }
#'
#' @param fit An `"scfusion"` object.
#' @param type `"in_space"` (default) or `"leave_one_out"`.
#' @param exclude_treated For `"in_space"`: drop the treated unit from every
#'   placebo donor pool. The default `FALSE` follows Abadie, Diamond and
#'   Hainmueller (2010), who move the treated unit into the donor pool; use
#'   `TRUE` to keep its treated target outcomes out of the placebo fits.
#' @param reselect_budget Re-run the budget search for each placebo fit
#'   (default) or hold the budget selected by `fit` fixed.
#' @param parallel,progress As in [scfusion_pi()].
#'
#' @return An object of class `"scfusion_placebo"` with
#' \describe{
#'   \item{`gaps`}{Long data frame of gaps with columns `domain`
#'     (`"reference"`/`"target"`), `unit`, `role` (`"treated"` or
#'     `"placebo"`), `time` and `gap`. For leave-one-out, `unit` names the
#'     donor left out and the full-donor fit has role `"treated"`.}
#'   \item{`estimates`}{One row per fit with the average target gap.}
#'   \item{`rmspe`}{(`"in_space"` only) RMSPE ratio table.}
#'   \item{`weights`}{List of fitted weight vectors.}
#' }
#' @references Abadie, A., Diamond, A. and Hainmueller, J. (2010). Synthetic
#'   control methods for comparative case studies. *Journal of the American
#'   Statistical Association*, 105, 493-505.
#' @examples
#' \donttest{
#' sim <- simulate_fusion_data(n_units = 12, n_ref = 10, n_target = 4,
#'                             seed = 2)
#' fit <- scfusion(sim$target, sim$reference, X = sim$X, Z = sim$Z,
#'                 budget_step = 0.1)
#' pl <- placebo_test(fit)
#' pl
#' plot(pl)
#' loo <- placebo_test(fit, type = "leave_one_out")
#' plot(loo, domain = "target")
#' }
#' @export
placebo_test <- function(fit, type = c("in_space", "leave_one_out"),
                         exclude_treated = FALSE, reselect_budget = TRUE,
                         parallel = FALSE, progress = interactive()) {
  if (!inherits(fit, "scfusion")) {
    stop("`fit` must be an object returned by scfusion().", call. = FALSE)
  }
  type <- match.arg(type)
  .check_flag(exclude_treated, "exclude_treated")
  .check_flag(reselect_budget, "reselect_budget")
  .check_flag(parallel, "parallel")
  .check_flag(progress, "progress")
  if (fit$settings$solver == "gurobi") .require_gurobi()

  data <- fit$data
  tr <- fit$treated
  dn <- fit$donors
  N <- nrow(data$target)
  budgets <- .refit_budgets(fit, reselect_budget, NULL)

  if (type == "in_space") {
    jobs <- lapply(dn, function(j) {
      pool <- setdiff(seq_len(N), j)
      if (exclude_treated) pool <- setdiff(pool, tr)
      list(target = j, pool = pool, label = data$unit_names[j])
    })
  } else {
    left_out <- dn[fit$weights > 0]
    jobs <- lapply(left_out, function(j) {
      list(target = tr, pool = setdiff(dn, j), label = data$unit_names[j])
    })
  }

  fits <- .cpf_lapply(jobs, function(job) {
    res <- tryCatch(.fit_unit(data, job$target, job$pool, fit$settings,
                              budgets),
                    error = function(e) NULL)
    if (is.null(res) || is.null(res$weights)) return(NULL)
    w <- res$weights
    names(w) <- data$unit_names[job$pool]
    list(
      weights = w,
      gap_F = data$reference[job$target, ] -
        as.numeric(crossprod(data$reference[job$pool, , drop = FALSE], w)),
      gap_Y = data$target[job$target, ] -
        as.numeric(crossprod(data$target[job$pool, , drop = FALSE], w))
    )
  }, parallel = parallel, progress = progress,
  label = if (type == "in_space") "Fitting in-space placebos" else
    "Fitting leave-one-out fits")

  ok <- !vapply(fits, is.null, logical(1))
  if (any(!ok)) {
    warning(sprintf("No feasible fit for: %s (skipped).",
                    paste(vapply(jobs[!ok], `[[`, "", "label"),
                          collapse = ", ")), call. = FALSE)
  }
  jobs <- jobs[ok]
  fits <- fits[ok]
  labels <- vapply(jobs, `[[`, "", "label")

  main <- list(
    weights = fit$weights,
    gap_F = data$reference[tr, ] - fit$synthetic_reference,
    gap_Y = fit$effects
  )
  treated_label <- data$unit_names[tr]
  all_fits <- c(list(main), fits)
  all_labels <- c(treated_label, labels)
  roles <- c("treated", rep("placebo", length(fits)))

  gaps <- do.call(rbind, lapply(seq_along(all_fits), function(k) {
    f <- all_fits[[k]]
    rbind(
      data.frame(domain = "reference", unit = all_labels[k], role = roles[k],
                 time = seq_along(f$gap_F), gap = as.numeric(f$gap_F),
                 stringsAsFactors = FALSE),
      data.frame(domain = "target", unit = all_labels[k], role = roles[k],
                 time = seq_along(f$gap_Y), gap = as.numeric(f$gap_Y),
                 stringsAsFactors = FALSE)
    )
  }))

  estimates <- data.frame(
    unit = all_labels, role = roles,
    average_gap = vapply(all_fits, function(f) mean(f$gap_Y), numeric(1)),
    rmspe_reference = vapply(all_fits, function(f) sqrt(mean(f$gap_F^2)),
                             numeric(1)),
    rmspe_target = vapply(all_fits, function(f) sqrt(mean(f$gap_Y^2)),
                          numeric(1)),
    stringsAsFactors = FALSE
  )
  if (type == "leave_one_out") {
    estimates$unit[1] <- "(all donors)"
    gaps$unit[gaps$role == "treated"] <- "(all donors)"
  }

  out <- structure(
    list(
      type = type,
      gaps = gaps,
      estimates = estimates,
      weights = stats::setNames(lapply(all_fits, `[[`, "weights"),
                                estimates$unit),
      exclude_treated = exclude_treated,
      reselect_budget = reselect_budget,
      treated = treated_label,
      fit = fit
    ),
    class = "scfusion_placebo"
  )
  if (type == "in_space") out$rmspe <- .rmspe_table(estimates)
  out
}

.rmspe_table <- function(estimates) {
  tab <- estimates[, c("unit", "role", "rmspe_reference", "rmspe_target",
                       "average_gap")]
  tab$ratio <- tab$rmspe_target / tab$rmspe_reference
  tab <- tab[order(-tab$ratio), , drop = FALSE]
  tab$rank <- seq_len(nrow(tab))
  rownames(tab) <- NULL
  tab
}

#' RMSPE-ratio placebo test
#'
#' Permutation-style test of Abadie, Diamond and Hainmueller (2015), adapted
#' to data fusion: the reference-domain fit plays the role of the
#' pre-intervention period,
#' \deqn{ratio_j = RMSPE_{target}(j) / RMSPE_{reference}(j),}
#' and the exact placebo p-value is the rank of the treated unit's ratio
#' (largest first) divided by the number of units.
#'
#' @param x An `"scfusion"` fit (in-space placebos are then computed with
#'   `exclude_treated = TRUE`) or an in-space `"scfusion_placebo"` object.
#' @param ... Passed to [placebo_test()] when `x` is a fit.
#'
#' @return An object of class `"scfusion_rmspe"` with the ratio table
#'   (`table`), the treated unit's `rank`, `n_units` and `p_value`, and a
#'   robustness check (`p_value_robust`) dropping units whose reference RMSPE
#'   is below a quarter of the median (tiny denominators inflate ratios).
#' @references Abadie, A., Diamond, A. and Hainmueller, J. (2015).
#'   Comparative politics and the synthetic control method. *American
#'   Journal of Political Science*, 59, 495-510.
#' @examples
#' \donttest{
#' sim <- simulate_fusion_data(n_units = 12, n_ref = 10, n_target = 4,
#'                             seed = 2)
#' fit <- scfusion(sim$target, sim$reference, X = sim$X, Z = sim$Z,
#'                 budget_step = 0.1)
#' rmspe_test(fit)
#' }
#' @export
rmspe_test <- function(x, ...) {
  if (inherits(x, "scfusion")) {
    args <- list(...)
    if (is.null(args$exclude_treated)) args$exclude_treated <- TRUE
    x <- do.call(placebo_test, c(list(fit = x, type = "in_space"), args))
  }
  if (!inherits(x, "scfusion_placebo") || x$type != "in_space") {
    stop("`x` must be an scfusion fit or an in-space placebo_test() result.",
         call. = FALSE)
  }
  tab <- x$rmspe
  n <- nrow(tab)
  rank <- tab$rank[tab$role == "treated"]
  tiny <- tab$rmspe_reference < 0.25 * stats::median(tab$rmspe_reference)
  keep <- !tiny | tab$role == "treated"
  sub <- tab[keep, , drop = FALSE]
  rank_sub <- rank(-sub$ratio, ties.method = "min")[sub$role == "treated"]
  structure(
    list(
      table = tab, rank = rank, n_units = n, p_value = rank / n,
      n_tiny = sum(tiny & tab$role != "treated"),
      rank_robust = rank_sub, n_robust = nrow(sub),
      p_value_robust = rank_sub / nrow(sub),
      exclude_treated = x$exclude_treated, treated = x$treated
    ),
    class = "scfusion_rmspe"
  )
}

#' @export
print.scfusion_rmspe <- function(x, digits = 4, ...) {
  cat("RMSPE-ratio placebo test (target RMSPE / reference RMSPE)\n")
  cat(sprintf("  treated unit excluded from placebo pools: %s\n\n",
              if (x$exclude_treated) "yes" else "no"))
  tab <- x$table
  num <- c("rmspe_reference", "rmspe_target", "average_gap", "ratio")
  tab[num] <- lapply(tab[num], signif, digits = digits)
  print(tab, row.names = FALSE)
  cat(sprintf("\nTreated unit (%s) rank: %d of %d -> exact placebo p-value = %s\n",
              x$treated, x$rank, x$n_units, .fmt(x$p_value, 3)))
  if (x$n_tiny > 0L) {
    cat(sprintf(paste0(
      "%d unit(s) have reference RMSPE below 25%% of the median; without ",
      "them the rank is %d of %d (p = %s).\n"),
      x$n_tiny, x$rank_robust, x$n_robust, .fmt(x$p_value_robust, 3)))
  }
  invisible(x)
}

#' @export
print.scfusion_placebo <- function(x, digits = 4, ...) {
  if (x$type == "in_space") {
    cat("In-space placebo test\n")
    cat(sprintf("  %d placebo fits | treated unit in placebo donor pools: %s\n\n",
                nrow(x$estimates) - 1L, if (x$exclude_treated) "no" else "yes"))
  } else {
    cat("Leave-one-out test (each positive-weight donor dropped in turn)\n\n")
  }
  est <- x$estimates
  num <- c("average_gap", "rmspe_reference", "rmspe_target")
  est[num] <- lapply(est[num], signif, digits = digits)
  print(est, row.names = FALSE)
  invisible(x)
}
