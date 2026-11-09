fd_long <- function(sim, ...) {
  fusion_data(sim$long, unit = "unit", time = "time", treated = "treated",
              domain = "domain", target = "target", reference = "reference",
              outcome = "outcome", covariates = c("cov1", "cov2", "cov3"), ...)
}
fd_wide <- function(sim, ...) {
  fusion_data(sim$wide, unit = "unit", time = "time", treated = "treated",
              outcome = c(target = "y_target", reference = "y_reference"),
              covariates = list(target = c("X1", "X2", "X3"),
                                reference = c("Z1", "Z2", "Z3")), ...)
}

test_that("both layouts rebuild the simulated matrices", {
  sim <- small_sim()
  for (fd in list(fd_long(sim, normalize = FALSE),
                  fd_wide(sim, normalize = FALSE))) {
    expect_s3_class(fd, "fusion_data")
    expect_equal(unname(fd$target), unname(sim$target))
    expect_equal(unname(fd$reference), unname(sim$reference))
    expect_equal(unname(fd$X), unname(sim$X))
    expect_equal(unname(fd$Z), unname(sim$Z))
    expect_equal(rownames(fd$target)[1], "treated")
    expect_equal(fd$donors, rownames(sim$target)[-1])
  }
  # target observed only in the last periods: earlier all-NA periods dropped
  expect_equal(colnames(fd_wide(sim)$target),
               as.character(seq_len(4) + 6))
})

test_that("scfusion() accepts fusion_data and matches the matrix call", {
  sim <- small_sim()
  fd <- fd_long(sim, normalize = FALSE)
  f1 <- scfusion(fd, budget_step = 0.1)
  f0 <- small_fit(sim)
  expect_equal(f1$estimate, f0$estimate)
  expect_equal(unname(f1$weights), unname(f0$weights))
  expect_error(scfusion(fd, sim$reference), "do not supply")
  expect_error(scfusion(fd, treated = 2), "do not supply")
  expect_equal(equi_confounding(fd)$estimate,
               equi_confounding(sim$target, sim$reference)$estimate)
  expect_error(equi_confounding(fd, sim$reference), "do not supply")
})

test_that("normalisation keeps original values for covariate_balance()", {
  sim <- small_sim()
  X_raw <- sim$X * 1000
  wide <- sim$wide
  wide[c("X1", "X2", "X3")] <- wide[c("X1", "X2", "X3")] * 1000
  fd <- fusion_data(wide, unit = "unit", time = "time", treated = "treated",
                    outcome = c(target = "y_target", reference = "y_reference"),
                    covariates = list(target = c("X1", "X2", "X3"),
                                      reference = c("Z1", "Z2", "Z3")),
                    target = "Hispanic", reference = "Black")
  expect_true(all(fd$X >= 0 & fd$X <= 1))
  expect_equal(unname(fd$X_raw), unname(X_raw))
  fit <- scfusion(fd, budget_step = 0.1)
  b <- covariate_balance(fit)
  expect_true(attr(b, "original_units"))
  expect_equal(b$treated[b$domain == "target (Hispanic)"],
               unname(X_raw[1, ]))
  expect_equal(unique(b$domain), c("reference (Black)", "target (Hispanic)"))
  expect_match(plot(fit)$labels$y, "target domain")
  expect_output(print(fd), "Hispanic")
})

test_that("covariates from a separate table, shared or domain-specific", {
  sim <- small_sim()
  long <- sim$long[, c("unit", "domain", "time", "outcome")]
  covs <- data.frame(unit = rownames(sim$X), sim$X, sim$Z)
  fd <- fusion_data(long, unit = "unit", time = "time", treated = "treated",
                    domain = "domain", target = "target",
                    reference = "reference", outcome = "outcome",
                    covariates = list(target = c("X1", "X2", "X3"),
                                      reference = c("Z1", "Z2", "Z3")),
                    covariate_data = covs, normalize = FALSE)
  expect_equal(unname(fd$X), unname(sim$X))
  expect_equal(unname(fd$Z), unname(sim$Z))
  # a character vector without a domain column: shared by both domains
  fd2 <- fusion_data(long, unit = "unit", time = "time", treated = "treated",
                     domain = "domain", target = "target",
                     reference = "reference", outcome = "outcome",
                     covariates = c("X1", "X2"), covariate_data = covs,
                     normalize = FALSE)
  expect_equal(unname(fd2$X), unname(fd2$Z))
})

test_that("periods and donors can be selected", {
  sim <- small_sim()
  fd <- fd_long(sim, target_periods = 9:10, reference_periods = 1:5,
                donors = c("donor_03", "donor_01", "treated", "donor_07"))
  expect_equal(colnames(fd$target), c("9", "10"))
  expect_equal(colnames(fd$reference), as.character(1:5))
  expect_equal(rownames(fd$target), c("treated", "donor_03", "donor_01",
                                      "donor_07"))
  expect_error(fd_long(sim, donors = c("donor_01", "nope")), "nope")
  expect_error(fd_long(sim, target_periods = 99), "not found")
})

test_that("time-varying covariates need covariate_time or a summary", {
  sim <- small_sim()
  long <- sim$long
  long$cov1 <- long$cov1 + long$time / 100
  expect_error(fusion_data(long, unit = "unit", time = "time",
                           treated = "treated", domain = "domain",
                           target = "target", reference = "reference",
                           outcome = "outcome", covariates = "cov1"),
               "more than one value")
  fd <- fusion_data(long, unit = "unit", time = "time", treated = "treated",
                    domain = "domain", target = "target",
                    reference = "reference", outcome = "outcome",
                    covariates = "cov1", covariate_time = 10,
                    normalize = FALSE)
  expect_equal(unname(fd$X[, 1]), unname(sim$X[, 1] + 0.1))
  expect_equal(unname(fd$Z[, 1]), unname(sim$Z[, 1] + 0.1))
  fd2 <- fusion_data(long, unit = "unit", time = "time", treated = "treated",
                     domain = "domain", target = "target",
                     reference = "reference", outcome = "outcome",
                     covariates = "cov1", covariate_summary = mean,
                     normalize = FALSE)
  expect_equal(unname(fd2$Z[, 1]), unname(sim$Z[, 1] + mean(1:10) / 100))
})

test_that("missing data: informative error or dropped donors", {
  sim <- small_sim()
  long <- sim$long
  long$outcome[long$unit == "donor_02" & long$domain == "target" &
                 long$time == 8] <- NA
  long <- long[!(long$unit == "donor_05" & long$domain == "reference" &
                   long$time == 3), ]
  args <- list(long, unit = "unit", time = "time", treated = "treated",
               domain = "domain", target = "target", reference = "reference",
               outcome = "outcome", covariates = c("cov1", "cov2", "cov3"))
  expect_error(do.call(fusion_data, args), "donor_02: target outcome \\(8\\)")
  expect_message(fd <- do.call(fusion_data, c(args, drop_incomplete = TRUE)),
                 "Dropped 2 donor")
  expect_false(any(c("donor_02", "donor_05") %in% rownames(fd$target)))
  expect_equal(fd$dropped, c("donor_02", "donor_05"))
  expect_output(print(fd), "dropped for missing data")
  bad <- sim$long
  bad$outcome[bad$unit == "treated"][1] <- NA
  expect_error(fusion_data(bad, unit = "unit", time = "time",
                           treated = "treated", domain = "domain",
                           target = "target", reference = "reference",
                           outcome = "outcome", covariates = "cov1",
                           drop_incomplete = TRUE), "treated unit has missing")
})

test_that("input checks give clear errors", {
  sim <- small_sim()
  expect_error(fusion_data(sim$long, unit = "unit", time = "time",
                           treated = "nobody", domain = "domain",
                           target = "target", reference = "reference",
                           outcome = "outcome", covariates = "cov1"),
               "not found")
  expect_error(fusion_data(sim$long, unit = "unit", time = "time",
                           treated = "treated", domain = "domain",
                           target = "Hispanic", reference = "reference",
                           outcome = "outcome", covariates = "cov1"),
               "Hispanic")
  expect_error(fusion_data(sim$wide, unit = "unit", time = "time",
                           treated = "treated", outcome = "y_target",
                           covariates = "X1"), "one column per domain")
  expect_error(fusion_data(sim$long, unit = "unit", time = "time",
                           treated = "treated", domain = "domain",
                           target = "target", reference = "reference",
                           outcome = "outcome", covariates = "nope"),
               "nope")
  dup <- rbind(sim$long, sim$long[1, ])
  expect_error(fusion_data(dup, unit = "unit", time = "time",
                           treated = "treated", domain = "domain",
                           target = "target", reference = "reference",
                           outcome = "outcome", covariates = "cov1"),
               "More than one")
})

test_that("downstream functions work on fusion_data fits", {
  fit <- scfusion(fd_long(small_sim()), budget_step = 0.1)
  pi <- scfusion_pi(fit, B = 5, seed = 1, progress = FALSE)
  expect_equal(nrow(pi_bands(pi)), 4L)
  expect_s3_class(placebo_test(fit, type = "leave_one_out",
                               progress = FALSE), "scfusion_placebo")
  expect_output(print(summary(fit)), "Covariate balance")
})
