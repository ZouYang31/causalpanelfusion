#' Candidate budget grid for the data-fusion weights
#'
#' Builds the grid of budget vectors \eqn{b = (b_F, b_Z, b_X)} searched by
#' [scfusion()]. Every candidate lies on the probability simplex and each
#' component is at least `min_value`. The grid is generated with exactly the same floating-point arithmetic as
#' the research code (`generate_b_list()`), so results are reproducible. With
#' the defaults (`step = 0.01`, `min_value = 0.01`) this gives 4,796
#' candidates: floating-point rounding drops 55 of the 4,851 points with
#' \eqn{b_X = 0.01}, exactly as in the original implementation.
#'
#' @param step Grid spacing (positive).
#' @param min_value Smallest value allowed for each component (positive).
#'
#' @return A numeric matrix with columns `b_F`, `b_Z` and `b_X`, one row per
#'   candidate budget.
#' @examples
#' grid <- budget_grid(step = 0.05)
#' nrow(grid)
#' head(grid)
#' @export
budget_grid <- function(step = 0.01, min_value = 0.01) {
  .check_scalar(step, "step", lower = 0, strict = TRUE)
  .check_scalar(min_value, "min_value", lower = 0, strict = TRUE)
  if (3 * min_value > 1 + sqrt(.Machine$double.eps)) {
    stop("`min_value` must be at most 1/3.", call. = FALSE)
  }
  out <- list()
  k <- 1L
  for (b_F in seq(min_value, 1 - min_value, by = step)) {
    if (1 - b_F < min_value) next
    for (b_Z in seq(min_value, 1 - b_F, by = step)) {
      b_X <- 1 - b_F - b_Z
      if (b_X >= min_value) {
        out[[k]] <- c(b_F, b_Z, b_X)
        k <- k + 1L
      }
    }
  }
  if (length(out) == 0L) {
    stop("No budget candidates satisfy the requested `step` and `min_value`.",
         call. = FALSE)
  }
  grid <- do.call(rbind, out)
  colnames(grid) <- c("b_F", "b_Z", "b_X")
  grid
}
