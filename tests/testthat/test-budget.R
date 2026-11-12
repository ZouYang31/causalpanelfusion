test_that("default grid reproduces the research grid", {
  g <- budget_grid()
  expect_equal(ncol(g), 3L)
  expect_equal(colnames(g), c("b_F", "b_Z", "b_X"))
  # generate_b_list() in the research code yields 4,796 candidates
  expect_equal(nrow(g), 4796L)
  expect_true(all(abs(rowSums(g) - 1) < 1e-12))
  expect_true(all(g >= 0.01 - 1e-12))
})

test_that("coarser grids and invalid input", {
  expect_equal(nrow(budget_grid(0.05)), 210L)
  expect_error(budget_grid(step = 0), "step")
  expect_error(budget_grid(min_value = 0.5), "1/3")
})
