#' Simulate data from the paper's linear factor model
#'
#' Generates reference- and target-domain panels for one treated unit and
#' `n_units - 1` donors from the interactive fixed-effects model used in the
#' simulation study of the paper:
#' \deqn{F_{it} = \rho_t + \phi_t' Z_i + \theta_t' \mu_i + \epsilon_{it},}
#' \deqn{Y_{is}(0) = \varrho_s + \varphi_s' X_i + \vartheta_s' \mu_i +
#'   \varepsilon_{is},}
#' with a shared latent confounder \eqn{\mu_i}. The treated unit (row 1)
#' receives `effect` in the target domain.
#'
#' With the default arguments and `seed = 45` the draws reproduce
#' `main_example.R` of the research code exactly.
#'
#' @param n_units Number of units (treated unit plus donors).
#' @param n_ref Number of reference-domain periods (\eqn{T}).
#' @param n_target Number of target-domain periods (\eqn{S}).
#' @param dr,dt,du Dimensions of \eqn{Z}, \eqn{X} and the latent \eqn{\mu}.
#' @param sd_ref,sd_target Standard deviations of the idiosyncratic shocks.
#' @param coef_max Upper limit of the uniform distribution of the factor
#'   loadings \eqn{\phi, \theta, \varphi, \vartheta}.
#' @param effect Per-period effects (length `n_target`, or a scalar). By
#'   default drawn as sorted \eqn{U(2, 5)} values.
#' @param seed Optional seed; the caller's random stream is left untouched.
#'
#' @return A list with `target` (observed target outcomes, effect added to
#'   unit 1), `target0` (untreated target potential outcomes), `reference`,
#'   `X`, `Z`, `mu`, `effect`, `psi` (the true average effect), `treated`
#'   (= 1) and `params` (the drawn coefficients).
#' @examples
#' sim <- simulate_fusion_data(seed = 45)
#' dim(sim$reference); dim(sim$target)
#' sim$psi
#' @export
simulate_fusion_data <- function(n_units = 31L, n_ref = 20L, n_target = 5L,
                                 dr = 3L, dt = 3L, du = 3L,
                                 sd_ref = sqrt(2), sd_target = sqrt(1 / 2),
                                 coef_max = 10, effect = NULL, seed = NULL) {
  .check_scalar(n_units, "n_units", lower = 3, integer = TRUE)
  .check_scalar(n_ref, "n_ref", lower = 1, integer = TRUE)
  .check_scalar(n_target, "n_target", lower = 1, integer = TRUE)
  .check_scalar(dr, "dr", lower = 1, integer = TRUE)
  .check_scalar(dt, "dt", lower = 1, integer = TRUE)
  .check_scalar(du, "du", lower = 1, integer = TRUE)
  .check_scalar(sd_ref, "sd_ref", lower = 0)
  .check_scalar(sd_target, "sd_target", lower = 0)
  .check_scalar(coef_max, "coef_max", lower = 0, strict = TRUE)
  if (!is.null(effect)) {
    if (length(effect) == 1L) effect <- rep(effect, n_target)
    if (!is.numeric(effect) || length(effect) != n_target ||
        any(!is.finite(effect))) {
      stop("`effect` must be a scalar or have length `n_target`.",
           call. = FALSE)
    }
  }
  N <- n_units; T <- n_ref; S <- n_target
  .with_seed(seed, {
    # draw order follows main_example.R
    mu <- matrix(stats::runif(N * du), N, du)
    Z <- matrix(stats::runif(N * dr), N, dr)
    X <- matrix(stats::runif(N * dt), N, dt)
    rho <- sort(stats::runif(T, 0, 20))
    phi <- matrix(stats::runif(dr * T, 0, coef_max), T, dr)
    theta <- matrix(stats::runif(du * T, 0, coef_max), T, du)
    eps <- matrix(stats::rnorm(N * T, 0, sd_ref), N, T)
    varrho <- sort(stats::runif(S, 0, 10))
    varphi <- matrix(stats::runif(dt * S, 0, coef_max), S, dt)
    vartheta <- matrix(stats::runif(du * S, 0, coef_max), S, du)
    vareps <- matrix(stats::rnorm(N * S, 0, sd_target), N, S)
    if (is.null(effect)) effect <- sort(stats::runif(S, 2, 5))
  })
  ones <- rep(1, N)
  Fm <- outer(ones, rho) + Z %*% t(phi) + mu %*% t(theta) + eps
  Y0 <- outer(ones, varrho) + X %*% t(varphi) + mu %*% t(vartheta) + vareps
  Y <- Y0
  Y[1L, ] <- Y0[1L, ] + effect
  units <- c("treated", sprintf("donor_%02d", seq_len(N - 1L)))
  rownames(Fm) <- rownames(Y) <- rownames(Y0) <- units
  rownames(Z) <- rownames(X) <- rownames(mu) <- units
  colnames(Fm) <- paste0("t", seq_len(T))
  colnames(Y) <- colnames(Y0) <- paste0("s", seq_len(S))
  colnames(Z) <- paste0("Z", seq_len(dr))
  colnames(X) <- paste0("X", seq_len(dt))
  list(
    target = Y, target0 = Y0, reference = Fm, X = X, Z = Z, mu = mu,
    effect = effect, psi = mean(effect), treated = 1L,
    params = list(rho = rho, phi = phi, theta = theta, varrho = varrho,
                  varphi = varphi, vartheta = vartheta)
  )
}
