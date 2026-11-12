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
#'   (= 1), `params` (the drawn coefficients), and the same observed data as
#'   two data frames for trying out [fusion_data()]:
#'   \describe{
#'     \item{`long`}{One row per unit, domain and period: columns `unit`,
#'       `domain` (`"target"`/`"reference"`), `time`, `outcome` and
#'       covariates `cov1`, `cov2`, ... (the `Z` values in reference rows,
#'       the `X` values in target rows).}
#'     \item{`wide`}{One row per unit and period: columns `unit`, `time`,
#'       `y_reference`, `y_target` and covariates `Z1`, ..., `X1`, ....}
#'   }
#'   The reference domain is observed in periods `1, ..., n_ref` and the
#'   target domain in the last `n_target` of them (`y_target` is `NA`
#'   before).
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
                  varphi = varphi, vartheta = vartheta),
    long = .sim_long(Y, Fm, X, Z),
    wide = .sim_wide(Y, Fm, X, Z)
  )
}

# Time indices: reference 1..T, target the last S of them.
.sim_times <- function(T, S) {
  list(ref = seq_len(T), tgt = seq_len(S) + max(T - S, 0L))
}

.sim_long <- function(Y, Fm, X, Z) {
  N <- nrow(Y)
  units <- rownames(Y)
  tm <- .sim_times(ncol(Fm), ncol(Y))
  K <- max(ncol(X), ncol(Z))
  pad <- function(M) {
    out <- matrix(NA_real_, N, K, dimnames = list(NULL, paste0("cov",
                                                                seq_len(K))))
    out[, seq_len(ncol(M))] <- M
    out
  }
  block <- function(M, C, dom, times) {
    nt <- length(times)
    data.frame(
      unit = rep(units, each = nt),
      domain = dom,
      time = rep(times, N),
      outcome = as.vector(t(M)),
      pad(C)[rep(seq_len(N), each = nt), , drop = FALSE],
      stringsAsFactors = FALSE
    )
  }
  out <- rbind(block(Fm, Z, "reference", tm$ref),
               block(Y, X, "target", tm$tgt))
  out <- out[order(match(out$unit, units), out$domain, out$time), ]
  rownames(out) <- NULL
  out
}

.sim_wide <- function(Y, Fm, X, Z) {
  N <- nrow(Y)
  units <- rownames(Y)
  tm <- .sim_times(ncol(Fm), ncol(Y))
  times <- sort(unique(c(tm$ref, tm$tgt)))
  nt <- length(times)
  iu <- rep(seq_len(N), each = nt)
  it <- rep(times, N)
  out <- data.frame(
    unit = units[iu], time = it,
    y_reference = Fm[cbind(iu, match(it, tm$ref))],
    y_target = Y[cbind(iu, match(it, tm$tgt))],
    stringsAsFactors = FALSE
  )
  out <- cbind(out, Z[iu, , drop = FALSE], X[iu, , drop = FALSE])
  rownames(out) <- NULL
  out
}
