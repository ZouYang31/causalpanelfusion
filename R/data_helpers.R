#' Reshape a long panel into a units-by-periods outcome matrix
#'
#' @param data A data frame in long format (one row per unit and period).
#' @param unit,time,outcome Column names (strings) of the unit identifier,
#'   the period and the outcome.
#' @param units Optional character vector giving the row order (for example
#'   the treated unit first). Defaults to the order of first appearance.
#'
#' @return A numeric matrix with one row per unit (row names = units) and one
#'   column per period (sorted).
#' @examples
#' long <- data.frame(city = rep(c("A", "B", "C"), each = 3),
#'                    month = rep(1:3, 3), rate = runif(9))
#' panel_to_matrix(long, unit = "city", time = "month", outcome = "rate")
#' @export
panel_to_matrix <- function(data, unit, time, outcome, units = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  for (nm in c(unit, time, outcome)) {
    if (!is.character(nm) || length(nm) != 1L || !nm %in% names(data)) {
      stop(sprintf("Column '%s' not found in `data`.", nm), call. = FALSE)
    }
  }
  u <- as.character(data[[unit]])
  tt <- data[[time]]
  if (anyDuplicated(paste(u, tt, sep = "\r"))) {
    stop("`data` has more than one row for some unit-period pair.",
         call. = FALSE)
  }
  if (is.null(units)) units <- unique(u)
  if (!setequal(units, unique(u))) {
    stop("`units` must contain exactly the units in `data`.", call. = FALSE)
  }
  periods <- sort(unique(tt))
  m <- matrix(NA_real_, length(units), length(periods),
              dimnames = list(units, as.character(periods)))
  m[cbind(match(u, units), match(tt, periods))] <- as.numeric(data[[outcome]])
  if (anyNA(m)) {
    stop("The panel is unbalanced or has missing outcomes.", call. = FALSE)
  }
  m
}

#' Extract a units-by-covariates matrix
#'
#' @param data A data frame with one row per unit.
#' @param unit Column name of the unit identifier.
#' @param covariates Character vector of covariate column names.
#' @param units Optional row order.
#' @param normalize If `TRUE`, min-max normalise each covariate (see
#'   [normalize_minmax()]).
#' @return A numeric matrix with row names = units.
#' @examples
#' df <- data.frame(city = c("A", "B", "C"), income = c(50, 60, 55),
#'                  age = c(30, 41, 35))
#' covariate_matrix(df, "city", c("income", "age"), normalize = TRUE)
#' @export
covariate_matrix <- function(data, unit, covariates, units = NULL,
                             normalize = FALSE) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  miss <- setdiff(c(unit, covariates), names(data))
  if (length(miss)) {
    stop(sprintf("Column(s) not found in `data`: %s.",
                 paste(miss, collapse = ", ")), call. = FALSE)
  }
  u <- as.character(data[[unit]])
  if (anyDuplicated(u)) {
    stop("`data` must have one row per unit.", call. = FALSE)
  }
  if (is.null(units)) units <- u
  if (!setequal(units, u)) {
    stop("`units` must contain exactly the units in `data`.", call. = FALSE)
  }
  m <- as.matrix(data[match(units, u), covariates, drop = FALSE])
  if (!is.numeric(m)) stop("Covariates must be numeric.", call. = FALSE)
  rownames(m) <- units
  if (normalize) m <- normalize_minmax(m)
  m
}

#' Min-max normalisation of columns
#'
#' Rescales every column to \eqn{[0, 1]} via \eqn{(x - \min x) / (\max x -
#' \min x)}, as done for the covariates in the paper's application. Constant
#' columns are set to zero.
#'
#' @param x A numeric matrix or data frame.
#' @param columns Optional subset of columns (names or indices) to rescale.
#' @return `x` with the selected columns rescaled.
#' @examples
#' normalize_minmax(cbind(a = c(1, 2, 3), b = c(10, 20, 40)))
#' @export
normalize_minmax <- function(x, columns = NULL) {
  if (is.null(columns)) columns <- seq_len(ncol(x))
  for (j in columns) {
    v <- x[, j]
    if (!is.numeric(v)) {
      stop("Only numeric columns can be normalised.", call. = FALSE)
    }
    r <- range(v, na.rm = TRUE)
    x[, j] <- if (diff(r) > 0) (v - r[1]) / diff(r) else 0
  }
  x
}
