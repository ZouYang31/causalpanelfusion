test_that("simulate_fusion_data is reproducible and well formed", {
  a <- simulate_fusion_data(seed = 45)
  b <- simulate_fusion_data(seed = 45)
  expect_identical(a, b)
  expect_equal(dim(a$reference), c(31L, 20L))
  expect_equal(dim(a$target), c(31L, 5L))
  expect_equal(unname(a$target[1, ] - a$target0[1, ]), a$effect)
  expect_equal(a$target[-1, ], a$target0[-1, ])
  expect_equal(a$psi, mean(a$effect))
})

test_that("a seed leaves the caller's random stream untouched", {
  set.seed(99)
  x1 <- runif(1)
  set.seed(99)
  invisible(simulate_fusion_data(seed = 1))
  x2 <- runif(1)
  expect_identical(x1, x2)
})

test_that("scalar effects are recycled", {
  s <- simulate_fusion_data(n_units = 5, n_target = 3, effect = 2, seed = 1)
  expect_equal(s$effect, rep(2, 3))
  expect_error(simulate_fusion_data(n_target = 3, effect = 1:2), "effect")
})
