test_that("prediction intervals have the expected structure", {
  fit <- small_fit()
  pi <- scfusion_pi(fit, B = 12, seed = 3, progress = FALSE)
  expect_s3_class(pi, "scfusion_pi")
  expect_equal(nrow(pi$intervals), 4L)
  expect_true(all(pi$intervals$lower <= pi$intervals$upper))
  expect_equal(pi$estimate, fit$estimate)
  expect_equal(pi$bootstrap$usable, 12L)
  expect_true(nrow(pi$placebo$gaps) >= 2L)
  expect_output(print(pi), "Prediction intervals")
})

test_that("intervals are reproducible and do not disturb the RNG", {
  fit <- small_fit()
  set.seed(5)
  a <- scfusion_pi(fit, B = 10, seed = 3, progress = FALSE)
  u1 <- runif(1)
  set.seed(5)
  b <- scfusion_pi(fit, B = 10, seed = 3, progress = FALSE)
  u2 <- runif(1)
  expect_identical(a$intervals, b$intervals)
  expect_identical(u1, u2)
})

test_that("interval formulas match the decomposition", {
  fit <- small_fit()
  pi <- scfusion_pi(fit, alpha1 = 0.1, alpha2 = 0.1, B = 15, seed = 2,
                    progress = FALSE)
  d <- pi$bootstrap$delta
  M1 <- unname(stats::quantile(d, c(0.05, 0.95)))
  avg <- rowMeans(pi$placebo$gaps)
  s <- stats::sd(avg - mean(avg))
  iv <- pi$intervals
  expect_equal(iv$lower[1], fit$estimate - M1[2])
  expect_equal(iv$upper[1], fit$estimate - M1[1])
  expect_equal(iv$lower[2], fit$estimate - M1[2] - s * qnorm(0.95))
  expect_equal(iv$upper[2], fit$estimate - M1[1] - s * qnorm(0.05))
  k <- sqrt(2 * log(2 / 0.1))
  expect_equal(iv$lower[4], fit$estimate - M1[2] - s * k)
  expect_equal(iv$upper[4], fit$estimate - M1[1] + s * k)
})

test_that("per-period bands and sensitivity analysis are coherent", {
  fit <- small_fit()
  pi <- scfusion_pi(fit, B = 12, seed = 3, progress = FALSE)
  for (m in c("gaussian", "empirical", "subgaussian")) {
    b <- pi_bands(pi, m)
    expect_equal(nrow(b), ncol(fit$data$target))
    expect_true(all(b$total_lo <= b$insample_lo + 1e-12))
    expect_true(all(b$total_hi >= b$insample_hi - 1e-12))
    expect_equal(b$effect, as.numeric(fit$effects))
  }
  sens <- sensitivity_pi(pi, rho = c(0, 1))
  expect_equal(sens$lower[2], pi$intervals$lower[2])
  expect_equal(sens$upper[2], pi$intervals$upper[2])
  expect_equal(sens$lower[1], pi$intervals$lower[1])
  expect_equal(sens$upper[1], pi$intervals$upper[1])
  sp <- sensitivity_pi(pi, rho = 1, method = "empirical", period = 2)
  b <- pi_bands(pi, "empirical")
  expect_equal(sp$lower, b$effect_lo[2])
  expect_equal(sp$upper, b$effect_hi[2])
  expect_s3_class(plot(sens), "ggplot")
})

test_that("bias allowance widens intervals symmetrically", {
  fit <- small_fit()
  a <- scfusion_pi(fit, B = 10, seed = 1, progress = FALSE)
  b <- scfusion_pi(fit, B = 10, seed = 1, bias_allowance = 0.5,
                   progress = FALSE)
  expect_equal(b$intervals$width[-1], a$intervals$width[-1] + 1)
})

test_that("fixed-budget refits and option checks", {
  fit <- small_fit()
  pi <- scfusion_pi(fit, B = 8, seed = 1, reselect_budget = FALSE,
                    progress = FALSE)
  expect_equal(pi$settings$n_budgets, 1L)
  expect_error(scfusion_pi(fit, alpha1 = 0.6, alpha2 = 0.5), "below 1")
  expect_error(scfusion_pi(list()), "scfusion")
})

test_that("parallel refits reproduce the serial result", {
  skip_on_cran()
  skip_if_not_installed("future.apply")
  fit <- small_fit()
  a <- scfusion_pi(fit, B = 6, seed = 4, progress = FALSE)
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)
  b <- scfusion_pi(fit, B = 6, seed = 4, parallel = TRUE, progress = FALSE)
  expect_identical(a$bootstrap$delta, b$bootstrap$delta)
  expect_identical(a$intervals, b$intervals)
})
