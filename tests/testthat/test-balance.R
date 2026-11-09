test_that("default table reproduces the fitted-scale balance", {
  fit <- small_fit()
  b <- covariate_balance(fit)
  expect_s3_class(b, "scfusion_balance")
  expect_equal(nrow(b), ncol(fit$data$Z) + ncol(fit$data$X))
  expect_equal(b$treated, fit$balance$treated)
  expect_equal(b$synthetic, fit$balance$synthetic)
  expect_equal(b$gap, b$gap_fit_scale)
  expect_equal(unname(attr(b, "nse")), unname(fit$nse[c("Z", "X")]))
  expect_false(attr(b, "original_units"))
})

test_that("covariates in original units are reported and aligned by name", {
  sim <- small_sim()
  X_raw <- sim$X * 1000 + 5
  Z_raw <- sim$Z * 50
  fit <- scfusion(sim$target, sim$reference, X = normalize_minmax(X_raw),
                  Z = normalize_minmax(Z_raw), budget_step = 0.1)
  b <- covariate_balance(fit, X = X_raw[rev(rownames(X_raw)), ], Z = Z_raw)
  w <- fit$weights
  tX <- b[b$domain == "target (X)", ]
  expect_equal(tX$treated, unname(X_raw["treated", ]))
  expect_equal(tX$synthetic, unname(colSums(X_raw[names(w), ] * w)))
  expect_equal(tX$donor_mean, unname(colMeans(X_raw[-1, ])))
  expect_equal(tX$gap_fit_scale,
               unname(fit$data$X["treated", ] -
                        colSums(fit$data$X[names(w), ] * w)))
  expect_equal(tX$treated_rank,
               unname(apply(X_raw, 2, function(v) rank(-v, ties.method = "min")[1])))
  expect_true(attr(b, "original_units"))
  expect_output(print(b), "Covariate balance")
  expect_output(print(b), "best attainable")
})

test_that("unnamed reporting matrices are matched by position", {
  sim <- small_sim()
  fit <- small_fit(sim)
  b <- covariate_balance(fit, X = unname(sim$X * 2), Z = unname(sim$Z))
  expect_equal(b$covariate, c(colnames(sim$Z), colnames(sim$X)))
  expect_equal(b$treated[b$domain == "target (X)"], unname(sim$X[1, ] * 2))
})

test_that("compare adds synthetic columns from other fits", {
  sim <- small_sim()
  fit <- small_fit(sim)
  fit2 <- small_fit(sim, eta_X = 0.01)
  b <- covariate_balance(fit, compare = list(tight = fit2))
  expect_true("synthetic_tight" %in% names(b))
  expect_equal(b$synthetic_tight, covariate_balance(fit2)$synthetic)
  expect_equal(unname(attr(b, "estimates")), c(fit$estimate, fit2$estimate))
  expect_output(print(b), "Estimates")
  b2 <- covariate_balance(fit, compare = fit2)
  expect_true("synthetic_fit2" %in% names(b2))
})

test_that("invalid inputs give clear errors", {
  sim <- small_sim()
  fit <- small_fit(sim)
  expect_error(covariate_balance(list()), "scfusion")
  expect_error(covariate_balance(fit, X = sim$X[, 1:2]), "same 3 columns")
  expect_error(covariate_balance(fit, X = sim$X[-1, ]), "Row names")
  expect_error(covariate_balance(fit, X = unname(sim$X)[-1, ]), "one row")
  other <- small_fit(small_sim(seed = 3, n_units = 10))
  expect_error(covariate_balance(fit, compare = list(other)), "same units")
  expect_error(covariate_balance(fit, compare = list(1)), "list of")
})

test_that("subsetting keeps a printable data frame", {
  b <- covariate_balance(small_fit())
  expect_output(print(b[1:2, ]), "treated")
})
