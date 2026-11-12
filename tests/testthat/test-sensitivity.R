test_that("eta sensitivity reproduces the fit at its own eta", {
  fit <- small_fit()
  es <- eta_sensitivity(fit, eta_Z = c(0.1, 0.3), eta_X = c(0.1, 0.3),
                        progress = FALSE)
  expect_s3_class(es, "scfusion_eta")
  expect_equal(nrow(es), 4L)
  row <- es[es$eta_Z == 0.1 & es$eta_X == 0.1, ]
  expect_equal(row$estimate, fit$estimate)
  expect_equal(row$NSE_F, fit$nse[["F"]])
  expect_true(all(c("absdiff_Z_Z1", "absdiff_X_X1") %in% names(es)))
  expect_s3_class(plot(es), "ggplot")
  expect_error(eta_sensitivity(fit, eta_Z = -1), "non-negative")
})
