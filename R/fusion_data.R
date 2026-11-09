#' Prepare data for synthetic control data fusion
#'
#' Builds the inputs of [scfusion()] from data frames, so that panel data can
#' be used in the layout it comes in. The result is passed to [scfusion()]
#' (and [equi_confounding()]) in place of the four matrices.
#'
#' @section Data layouts:
#' \describe{
#'   \item{Domain as a column (`domain` given)}{One row per unit, domain and
#'     period, with the same variable names in both domains, e.g. columns
#'     `city`, `race`, `month`, `rate`, `income`. `target` and `reference`
#'     are the values of the `domain` column that identify the two domains
#'     (e.g. `"Hispanic"` and `"Black"`), `outcome` is a single column
#'     name.}
#'   \item{Domains as columns (`domain` not given)}{One row per unit and
#'     period, with one outcome column per domain, e.g. `vac_black` and
#'     `vac_hisp`. `outcome` names both columns:
#'     `c(target = "vac_hisp", reference = "vac_black")`. `target` and
#'     `reference` are then optional labels for the domains.}
#' }
#' Covariates can be given in either way:
#' \itemize{
#'   \item a character vector: the same column names are used in both
#'     domains. With a `domain` column in the covariate source, each domain
#'     reads its own rows; otherwise the values are shared by both domains.
#'   \item a list `list(target = ..., reference = ...)` of domain-specific
#'     column names, e.g. `income_hisp` and `income_black`.
#' }
#' Covariates are read from `data` or, if supplied, from `covariate_data`
#' (for example a table with one row per unit, or per unit and domain). It
#' must use the same `unit` column name, and the same `domain` column name
#' if it has one.
#'
#' @section Periods:
#' Each domain uses its own periods: by default every value of `time` at
#' which its outcome is observed for at least one unit (so periods in which
#' a domain is not observed at all are dropped). `target_periods` and
#' `reference_periods` select periods explicitly. Periods are ordered with
#' [sort()], so the `time` column must sort chronologically (numbers,
#' dates, or strings such as `"2021-03"`).
#'
#' @section Covariates that vary over rows:
#' The method uses one value per unit and domain. By default each covariate
#' must be constant across a unit's rows (as when unit-level covariates are
#' repeated in every period). Otherwise choose the rows with
#' `covariate_time` (values of the `time` column of the covariate source)
#' and/or combine them with `covariate_summary` (e.g. `mean`).
#'
#' @param data A data frame in either layout (see Data layouts).
#' @param unit,time Names of the unit and time columns.
#' @param treated The treated unit (a value of the `unit` column).
#' @param outcome Outcome column: one name with a `domain` column, or a
#'   named pair `c(target = , reference = )` without one.
#' @param covariates Covariate columns: a character vector of names used in
#'   both domains, or `list(target = , reference = )` with domain-specific
#'   names.
#' @param domain Optional name of the column identifying the domain.
#' @param target,reference With `domain`: the values of the `domain`
#'   column for the target and reference domain. Without `domain`: optional
#'   labels used in printed output and plots.
#' @param covariate_data Optional data frame holding the covariates.
#' @param covariate_time Optional values of the `time` column; only these
#'   rows of the covariate source are used.
#' @param covariate_summary Optional function combining a unit's covariate
#'   values over rows, e.g. `mean`.
#' @param target_periods,reference_periods Optional periods to use in each
#'   domain.
#' @param donors Optional donor units (values of `unit`). By default every
#'   unit other than the treated one.
#' @param normalize Min-max normalise each covariate to `[0, 1]` across the
#'   units in the analysis, separately in each domain (default `TRUE`, as in
#'   the paper). The original values are kept, and [covariate_balance()]
#'   reports them.
#' @param drop_incomplete Units with missing outcomes or covariates cause an
#'   error by default. With `TRUE` such donors are dropped with a message
#'   (the treated unit must be complete).
#'
#' @return An object of class `"fusion_data"` holding the outcome and
#'   covariate matrices (units in rows, treated unit first), the original
#'   covariate values, domain labels and the settings used. `print()`
#'   summarises it.
#' @seealso [scfusion()], [covariate_balance()]; [panel_to_matrix()] and
#'   [covariate_matrix()] for building matrices by hand.
#' @examples
#' sim <- simulate_fusion_data(n_units = 16, n_ref = 12, n_target = 4,
#'                             seed = 1)
#'
#' # Domain as a column: one row per unit, domain and period
#' head(sim$long)
#' fd <- fusion_data(sim$long, unit = "unit", time = "time",
#'                   treated = "treated", domain = "domain",
#'                   target = "target", reference = "reference",
#'                   outcome = "outcome",
#'                   covariates = c("cov1", "cov2", "cov3"))
#' fd
#'
#' # Domains as columns: one row per unit and period
#' head(sim$wide)
#' fd2 <- fusion_data(sim$wide, unit = "unit", time = "time",
#'                    treated = "treated",
#'                    outcome = c(target = "y_target",
#'                                reference = "y_reference"),
#'                    covariates = list(target = c("X1", "X2", "X3"),
#'                                      reference = c("Z1", "Z2", "Z3")))
#'
#' fit <- scfusion(fd, budget_step = 0.1)
#' fit
#' covariate_balance(fit)   # covariates in their original units
#' @export
fusion_data <- function(data, unit, time, treated, outcome, covariates,
                        domain = NULL, target = NULL, reference = NULL,
                        covariate_data = NULL, covariate_time = NULL,
                        covariate_summary = NULL,
                        target_periods = NULL, reference_periods = NULL,
                        donors = NULL, normalize = TRUE,
                        drop_incomplete = FALSE) {
  cl <- match.call()
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  .check_name(unit, "unit")
  .check_name(time, "time")
  .check_cols(data, c(unit, time), "data")
  .check_flag(normalize, "normalize")
  .check_flag(drop_incomplete, "drop_incomplete")
  if (!is.null(covariate_summary) && !is.function(covariate_summary)) {
    stop("`covariate_summary` must be a function, e.g. mean.", call. = FALSE)
  }
  if (missing(treated) || length(treated) != 1L || is.na(treated)) {
    stop("`treated` must be a single unit identifier.", call. = FALSE)
  }
  treated <- as.character(treated)

  # ---- layout ---------------------------------------------------------------
  long <- !is.null(domain)
  if (long) {
    .check_name(domain, "domain")
    .check_cols(data, domain, "data")
    if (is.null(target) || is.null(reference) || length(target) != 1L ||
        length(reference) != 1L) {
      stop("With `domain`, give the values identifying the two domains in ",
           "`target` and `reference` (e.g. target = \"Hispanic\").",
           call. = FALSE)
    }
    target <- as.character(target)
    reference <- as.character(reference)
    if (identical(target, reference)) {
      stop("`target` and `reference` must differ.", call. = FALSE)
    }
    dom <- as.character(data[[domain]])
    for (v in c(target, reference)) {
      if (!v %in% dom) {
        stop(sprintf("Value '%s' not found in column '%s'.", v, domain),
             call. = FALSE)
      }
    }
    .check_name(outcome, "outcome")
    .check_cols(data, outcome, "data")
    outcome <- c(target = outcome, reference = outcome)
    labels <- c(target = target, reference = reference)
    axis <- c(target = sprintf("%s, %s", outcome[[1]], target),
              reference = sprintf("%s, %s", outcome[[1]], reference))
  } else {
    if (!is.character(outcome) || length(outcome) != 2L ||
        !setequal(names(outcome), c("target", "reference"))) {
      stop("Without `domain`, `outcome` must name one column per domain: ",
           "c(target = \"...\", reference = \"...\").", call. = FALSE)
    }
    outcome <- outcome[c("target", "reference")]
    .check_cols(data, outcome, "data")
    labels <- c(
      target = if (is.null(target)) outcome[["target"]] else
        as.character(target)[1],
      reference = if (is.null(reference)) outcome[["reference"]] else
        as.character(reference)[1])
    axis <- labels
  }

  # ---- covariate specification ------------------------------------------------
  if (is.list(covariates)) {
    if (!setequal(names(covariates), c("target", "reference"))) {
      stop("A list of covariates must have elements `target` and ",
           "`reference`.", call. = FALSE)
    }
    cov_names <- lapply(covariates[c("target", "reference")], as.character)
  } else if (is.character(covariates)) {
    cov_names <- list(target = covariates, reference = covariates)
  } else {
    stop("`covariates` must be a character vector or a list with ",
         "`target` and `reference`.", call. = FALSE)
  }
  if (any(lengths(cov_names) == 0L)) {
    stop("Each domain needs at least one covariate.", call. = FALSE)
  }
  cov_src <- if (is.null(covariate_data)) data else covariate_data
  if (!is.data.frame(cov_src)) {
    stop("`covariate_data` must be a data frame.", call. = FALSE)
  }
  .check_cols(cov_src, unit, "covariate_data")
  .check_cols(cov_src, unique(unlist(cov_names)),
              if (is.null(covariate_data)) "data" else "covariate_data")
  if (!is.null(covariate_time)) {
    .check_cols(cov_src, time,
                if (is.null(covariate_data)) "data" else "covariate_data")
  }
  cov_by_domain <- long && domain %in% names(cov_src)

  # ---- units -------------------------------------------------------------------
  all_units <- unique(as.character(data[[unit]]))
  if (!treated %in% all_units) {
    stop(sprintf("Treated unit '%s' not found in column '%s'.", treated, unit),
         call. = FALSE)
  }
  if (is.null(donors)) {
    donors <- setdiff(all_units, treated)
  } else {
    donors <- unique(as.character(donors))
    unknown <- setdiff(donors, all_units)
    if (length(unknown)) {
      stop(sprintf("Donor(s) not found in column '%s': %s.", unit,
                   .collapse(unknown)), call. = FALSE)
    }
    donors <- setdiff(donors, treated)
  }
  units <- c(treated, donors)

  # ---- outcomes --------------------------------------------------------------------
  u_data <- as.character(data[[unit]])
  keep <- u_data %in% units
  rows_t <- data[keep & (if (long) dom == target else TRUE), , drop = FALSE]
  rows_r <- data[keep & (if (long) dom == reference else TRUE), , drop = FALSE]
  Y <- .outcome_grid(rows_t, unit, time, outcome[["target"]], units,
                     "target")
  Fm <- .outcome_grid(rows_r, unit, time, outcome[["reference"]], units,
                      "reference")
  Y <- .select_periods(Y, target_periods, "target")
  Fm <- .select_periods(Fm, reference_periods, "reference")

  # ---- covariates -----------------------------------------------------------------
  X_raw <- .covariate_block(cov_src, unit, time, units, cov_names$target,
                            if (cov_by_domain) domain else NULL, target,
                            covariate_time, covariate_summary)
  Z_raw <- .covariate_block(cov_src, unit, time, units, cov_names$reference,
                            if (cov_by_domain) domain else NULL, reference,
                            covariate_time, covariate_summary)

  # ---- missing data -------------------------------------------------------------
  issues <- .missing_issues(units, list(
    "target outcome" = Y, "reference outcome" = Fm,
    "target covariates" = X_raw, "reference covariates" = Z_raw))
  dropped <- character()
  if (length(issues)) {
    if (treated %in% names(issues)) {
      stop(sprintf("The treated unit has missing data: %s.",
                   issues[[treated]]), call. = FALSE)
    }
    if (!drop_incomplete) {
      shown <- utils::head(names(issues), 8L)
      stop(sprintf(
        "%d donor(s) have missing data (use drop_incomplete = TRUE to drop them):\n%s%s",
        length(issues),
        paste0("  ", shown, ": ", unlist(issues[shown]), collapse = "\n"),
        if (length(issues) > 8L) "\n  ..." else ""), call. = FALSE)
    }
    dropped <- names(issues)
    message(sprintf("Dropped %d donor(s) with missing data: %s.",
                    length(dropped), .collapse(dropped)))
    units <- setdiff(units, dropped)
    donors <- setdiff(donors, dropped)
    Y <- Y[units, , drop = FALSE]
    Fm <- Fm[units, , drop = FALSE]
    X_raw <- X_raw[units, , drop = FALSE]
    Z_raw <- Z_raw[units, , drop = FALSE]
  }
  if (length(donors) < 2L) {
    stop("At least two donors with complete data are needed.", call. = FALSE)
  }

  structure(
    list(
      target = Y, reference = Fm,
      X = if (normalize) normalize_minmax(X_raw) else X_raw,
      Z = if (normalize) normalize_minmax(Z_raw) else Z_raw,
      X_raw = X_raw, Z_raw = Z_raw,
      treated = treated, donors = donors, dropped = dropped,
      labels = labels, axis_labels = axis, outcome = outcome,
      covariates = cov_names, normalize = normalize,
      layout = if (long) "domain column" else "domains as columns",
      call = cl
    ),
    class = "fusion_data"
  )
}

#' @export
print.fusion_data <- function(x, ...) {
  cat("Data for synthetic control data fusion\n")
  cat(sprintf("  layout: %s\n", x$layout))
  cat(sprintf("  treated unit: %s | donors: %d\n", x$treated,
              length(x$donors)))
  per <- function(M) {
    p <- colnames(M)
    if (length(p) == 1L) return(sprintf("1 period (%s)", p))
    sprintf("%d periods (%s to %s)", length(p), p[1], p[length(p)])
  }
  for (d in c("target", "reference")) {
    M <- if (d == "target") x$target else x$reference
    lab <- if (identical(tolower(x$labels[[d]]), d)) "" else
      sprintf(" (%s)", x$labels[[d]])
    cat(sprintf("  %s domain%s: outcome '%s', %s\n", d, lab,
                x$outcome[[d]], per(M)))
  }
  for (d in c("target", "reference")) {
    cn <- x$covariates[[d]]
    cat(sprintf("  %s covariates (%d): %s\n", d, length(cn),
                paste(cn, collapse = ", ")))
  }
  if (x$normalize) {
    cat(sprintf(paste0("  covariates min-max normalised across the %d units ",
                       "(original values kept for covariate_balance())\n"),
                nrow(x$target)))
  }
  if (length(x$dropped)) {
    cat(sprintf("  dropped for missing data: %s\n", .collapse(x$dropped)))
  }
  invisible(x)
}

# ---- helpers -------------------------------------------------------------------

.check_name <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf("`%s` must be a single column name.", name), call. = FALSE)
  }
}

.check_cols <- function(df, cols, where) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) {
    stop(sprintf("Column(s) not found in `%s`: %s.", where,
                 .collapse(miss)), call. = FALSE)
  }
}

.collapse <- function(x, max = 10L) {
  x <- as.character(x)
  if (length(x) > max) x <- c(x[seq_len(max)], "...")
  paste(x, collapse = ", ")
}

# Units x periods matrix of one domain's outcome (NA where absent).
.outcome_grid <- function(rows, unit, time, outcome, units, label) {
  y <- rows[[outcome]]
  if (!is.numeric(y)) {
    stop(sprintf("Outcome column '%s' must be numeric.", outcome),
         call. = FALSE)
  }
  u <- as.character(rows[[unit]])
  tt <- rows[[time]]
  key <- paste(u, as.character(tt), sep = "\r")
  if (anyDuplicated(key)) {
    d <- strsplit(key[duplicated(key)][1], "\r", fixed = TRUE)[[1]]
    stop(sprintf(paste0("More than one %s-domain row for unit '%s' at time ",
                        "'%s'."), label, d[1], d[2]), call. = FALSE)
  }
  periods <- as.character(sort(unique(tt)))
  m <- matrix(NA_real_, length(units), length(periods),
              dimnames = list(units, periods))
  if (length(u)) {
    m[cbind(match(u, units), match(as.character(tt), periods))] <-
      as.numeric(y)
  }
  m
}

.select_periods <- function(M, periods, label) {
  if (is.null(periods)) {
    obs <- colSums(!is.na(M)) > 0
    if (!any(obs)) {
      stop(sprintf("No observed %s-domain outcomes.", label), call. = FALSE)
    }
    return(M[, obs, drop = FALSE])
  }
  p <- as.character(periods)
  miss <- setdiff(p, colnames(M))
  if (length(miss)) {
    stop(sprintf("Requested %s periods not found: %s.", label,
                 .collapse(miss)), call. = FALSE)
  }
  M[, colnames(M) %in% p, drop = FALSE]
}

# Units x covariates matrix with one value per unit.
.covariate_block <- function(src, unit, time, units, cols, domain_col,
                             domain_value, covariate_time, summary_fun) {
  rows <- src
  if (!is.null(domain_col)) {
    rows <- rows[as.character(rows[[domain_col]]) == domain_value, ,
                 drop = FALSE]
  }
  if (!is.null(covariate_time)) {
    rows <- rows[as.character(rows[[time]]) %in% as.character(covariate_time),
                 , drop = FALSE]
  }
  u <- as.character(rows[[unit]])
  rows <- rows[u %in% units, , drop = FALSE]
  u <- factor(as.character(rows[[unit]]), levels = units)
  out <- matrix(NA_real_, length(units), length(cols),
                dimnames = list(units, cols))
  varying <- character()
  for (k in cols) {
    v <- rows[[k]]
    if (!is.numeric(v)) {
      stop(sprintf("Covariate column '%s' must be numeric.", k),
           call. = FALSE)
    }
    groups <- split(v, u)
    for (i in seq_along(units)) {
      g <- groups[[i]]
      g <- g[!is.na(g)]
      if (!length(g)) next
      if (is.null(summary_fun)) {
        ug <- unique(g)
        if (length(ug) > 1L) {
          varying <- c(varying, sprintf("%s (%s)", k, units[i]))
          next
        }
        out[i, k] <- ug
      } else {
        val <- summary_fun(g)
        if (!is.numeric(val) || length(val) != 1L) {
          stop("`covariate_summary` must return a single number.",
               call. = FALSE)
        }
        out[i, k] <- val
      }
    }
  }
  if (length(varying)) {
    stop(sprintf(paste0(
      "Covariates take more than one value for some units: %s. Choose rows ",
      "with `covariate_time` or combine them with `covariate_summary` ",
      "(e.g. mean)."), .collapse(varying, 5L)), call. = FALSE)
  }
  out
}

# Describe missing values per unit; returns a named list (unit -> text).
.missing_issues <- function(units, blocks) {
  out <- list()
  for (u in units) {
    parts <- character()
    for (nm in names(blocks)) {
      M <- blocks[[nm]]
      miss <- colnames(M)[is.na(M[u, ])]
      if (length(miss)) {
        parts <- c(parts, sprintf("%s (%s)", nm, .collapse(miss, 4L)))
      }
    }
    if (length(parts)) out[[u]] <- paste(parts, collapse = "; ")
  }
  out
}
