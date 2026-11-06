test_that("panel_to_matrix reshapes long panels", {
  long <- data.frame(city = rep(c("B", "A"), each = 3), month = rep(3:1, 2),
                     y = 1:6)
  m <- panel_to_matrix(long, "city", "month", "y")
  expect_equal(rownames(m), c("B", "A"))
  expect_equal(colnames(m), c("1", "2", "3"))
  expect_equal(unname(m["B", ]), c(3, 2, 1))
  expect_error(panel_to_matrix(long[-1, ], "city", "month", "y"),
               "unbalanced")
  expect_error(panel_to_matrix(long, "town", "month", "y"), "town")
})

test_that("covariate_matrix and normalize_minmax", {
  df <- data.frame(id = c("a", "b", "c"), x = c(1, 3, 5), z = c(2, 2, 2))
  m <- covariate_matrix(df, "id", c("x", "z"), units = c("c", "a", "b"))
  expect_equal(rownames(m), c("c", "a", "b"))
  n <- normalize_minmax(m)
  expect_equal(unname(n[, "x"]), c(1, 0, 0.5))
  expect_equal(unname(n[, "z"]), c(0, 0, 0))
  expect_error(covariate_matrix(df, "id", "w"), "w")
})

test_that("plot methods return ggplot objects", {
  fit <- small_fit()
  expect_s3_class(plot(fit), "ggplot")
  expect_s3_class(plot(fit, type = "reference"), "ggplot")
  expect_s3_class(plot(fit, type = "gaps"), "ggplot")
  pi <- scfusion_pi(fit, B = 6, seed = 1, progress = FALSE)
  expect_s3_class(plot(pi), "ggplot")
  expect_s3_class(plot(pi, type = "effect", method = "empirical"), "ggplot")
  expect_s3_class(plot(fit, bands = pi, band_style = "ribbon"), "ggplot")
})
