test_that("equi-confounding estimators follow the paper's formulas", {
  sim <- small_sim()
  Y <- sim$target
  Fm <- sim$reference
  ec <- equi_confounding(Y, Fm)
  lin <- (mean(Y[1, ]) - mean(Y[-1, ])) - (mean(Fm[1, ]) - mean(Fm[-1, ]))
  lg <- mean(Y[1, ]) - mean(Fm[1, ]) / mean(Fm[-1, ]) * mean(Y[-1, ])
  expect_equal(ec$estimate[["linear"]], lin)
  expect_equal(ec$estimate[["log"]], lg)
  expect_named(equi_confounding(Y, Fm, type = "log")$estimate, "log")
  expect_equal(equi_confounding(Y, Fm, treated = "treated")$estimate,
               ec$estimate)
  expect_output(print(ec), "linear")
  expect_s3_class(plot(ec), "ggplot")
})
