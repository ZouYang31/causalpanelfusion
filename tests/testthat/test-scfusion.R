test_that("weights live on the simplex and the estimate is consistent", {
  sim <- small_sim()
  fit <- small_fit(sim)
  w <- fit$weights
  expect_s3_class(fit, "scfusion")
  expect_true(all(w >= 0))
  expect_true(sum(w) <= 1 + 1e-6 && sum(w) > 0.95)
  expect_named(w, rownames(sim$target)[-1])
  expect_equal(fit$estimate, mean(fit$effects))
  expect_equal(
    fit$estimate,
    mean(sim$target[1, ]) - sum(w * rowMeans(sim$target[-1, ]))
  )
  expect_equal(sum(fit$budget), 1)
  expect_equal(unname(coef(fit)), unname(w))
})

test_that("covariate constraints hold at the solution", {
  fit <- small_fit()
  s <- fit$settings
  bZ <- (1 + s$eta_Z) * (1 + fit$baseline[["Z"]]) - 1
  bX <- (1 + s$eta_X) * (1 + fit$baseline[["X"]]) - 1
  # thresholding small weights can move NSE slightly; allow a little slack
  expect_lt(fit$nse[["Z"]], bZ + 0.01)
  expect_lt(fit$nse[["X"]], bX + 0.01)
})

test_that("results match the original research implementation", {
  # Reference values from the research code (find_best_B with Gurobi,
  # baselines from kernlab::ipop) on main_example.R data, grid step 0.05.
  FIXTURE_W <- c(0, 0, 0, 0, 0.22327289, 0.17987033, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0.21731518, 0.2351197, 0.09418022, 0, 0, 0, 0, 0, 0, 0,
                 0.05024169, 0, 0, 0, 0)
  FIXTURE_PSI <- 3.5087932954775
  FIXTURE_BUDGET <- c(b_F = 0.96, b_Z = 0.01, b_X = 0.03)
  FIXTURE_BASE <- c(0.00334825048636344, 8.04765111183881e-05)
  sim <- simulate_fusion_data(seed = 45)
  fit <- scfusion(sim$target, sim$reference, X = sim$X, Z = sim$Z,
                  budget_step = 0.05)
  expect_equal(unname(fit$weights), FIXTURE_W, tolerance = 1e-4)
  expect_equal(fit$estimate, FIXTURE_PSI, tolerance = 1e-4)
  expect_equal(fit$budget, FIXTURE_BUDGET, tolerance = 1e-8)
  # ipop (research) solves the baseline QP only to ~5 significant digits
  expect_lt(max(abs(unname(fit$baseline) - FIXTURE_BASE)), 1e-4)
})

test_that("row names align inputs and the treated unit can be named", {
  sim <- small_sim()
  fit1 <- small_fit(sim)
  perm <- c(3, 1, 2, 4:nrow(sim$X))
  fit2 <- scfusion(sim$target, sim$reference[rev(rownames(sim$reference)), ],
                   X = sim$X[perm, ], Z = sim$Z[rev(seq_len(nrow(sim$Z))), ],
                   treated = "treated", budget_step = 0.1)
  expect_equal(fit1$weights, fit2$weights)
  expect_equal(fit1$estimate, fit2$estimate)
})

test_that("the treated unit need not be the first row", {
  sim <- small_sim()
  fit1 <- small_fit(sim)
  ord <- c(2:5, 1, 6:nrow(sim$target))
  fit2 <- scfusion(sim$target[ord, ], sim$reference[ord, ], X = sim$X[ord, ],
                   Z = sim$Z[ord, ], treated = 5, budget_step = 0.1)
  expect_equal(fit2$data$unit_names[fit2$treated], "treated")
  expect_equal(fit1$estimate, fit2$estimate)
  expect_equal(fit1$weights[sort(names(fit1$weights))],
               fit2$weights[sort(names(fit2$weights))])
})

test_that("unnamed matrices and data frames are accepted", {
  sim <- small_sim()
  fit <- scfusion(unname(sim$target), as.data.frame(unname(sim$reference)),
                  X = unname(sim$X), Z = unname(sim$Z), budget_step = 0.1)
  expect_equal(fit$estimate, small_fit(sim)$estimate)
  expect_true(all(grepl("^unit_", names(fit$weights))))
})

test_that("custom budgets and keep_path work", {
  sim <- small_sim()
  b <- rbind(c(b_F = 0.9, b_Z = 0.05, b_X = 0.05),
             c(b_F = 0.5, b_Z = 0.25, b_X = 0.25))
  fit <- scfusion(sim$target, sim$reference, X = sim$X, Z = sim$Z,
                  budgets = b, keep_path = TRUE)
  expect_equal(nrow(fit$path), 2L)
  expect_equal(fit$nse[["F"]], min(fit$path[, "NSE_F"], na.rm = TRUE))
  expect_error(scfusion(sim$target, sim$reference, X = sim$X, Z = sim$Z,
                        budgets = cbind(1, 1, 1)), "sum to one")
})

test_that("input validation gives clear errors", {
  sim <- small_sim()
  expect_error(scfusion(sim$target[-1, ], sim$reference, sim$X, sim$Z),
               "row names")
  expect_error(scfusion(unname(sim$target)[-1, ], unname(sim$reference),
                        unname(sim$X), unname(sim$Z)), "same number")
  bad <- sim$target
  bad[2, 2] <- NA
  expect_error(scfusion(bad, sim$reference, sim$X, sim$Z), "non-finite")
  expect_error(scfusion(sim$target, sim$reference, sim$X, sim$Z,
                        treated = "nope"), "treated")
  expect_error(scfusion(sim$target, sim$reference, sim$X, sim$Z,
                        eta_Z = -1), "eta_Z")
})

test_that("print and summary methods run", {
  fit <- small_fit()
  expect_output(print(fit), "Synthetic control data fusion")
  s <- summary(fit)
  expect_s3_class(s, "summary.scfusion")
  expect_output(print(s), "Covariate balance")
  expect_equal(nrow(s$effects), ncol(fit$data$target))
})

test_that("ECOS calls cannot corrupt package constants", {
  # ECOSolveR rescales its inputs in place; the package passes copies.
  sim <- simulate_fusion_data(seed = 45)
  for (k in 1:3) {
    fit <- scfusion(sim$target, sim$reference, X = sim$X, Z = sim$Z,
                    budget_step = 0.05)
  }
  expect_equal(fit$estimate, 3.5087932954775, tolerance = 1e-4)
  m <- matrix(1, nrow = 1, ncol = 3)
  expect_equal(nrow(m), 1L)
})

test_that("solver = 'gurobi' fails informatively without Gurobi", {
  skip_if(requireNamespace("gurobi", quietly = TRUE))
  sim <- small_sim()
  expect_error(scfusion(sim$target, sim$reference, sim$X, sim$Z,
                        solver = "gurobi"), "Gurobi")
})
