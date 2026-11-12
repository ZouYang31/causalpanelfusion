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

legend_labels <- function(p, aes = "colour") {
  b <- ggplot2::ggplot_build(p)
  b$plot$scales$get_scales(aes)$get_labels()
}

test_that("legend entries cover every coloured element", {
  sim <- small_sim()
  fit <- small_fit(sim)
  pi <- scfusion_pi(fit, B = 6, seed = 1, progress = FALSE)
  p <- plot(pi, counterfactual = sim$target0[1, ])
  labs <- legend_labels(p)
  expect_true(all(c("treated", "Synthetic treated", "Counterfactual treated",
                    "Donor pool units") %in% labs))
  expect_true(any(grepl("^Weight construction error", labs)))
  expect_true(any(grepl("^Prediction interval", labs)))
  pe <- plot(pi, type = "effect", counterfactual = sim$target0[1, ])
  expect_equal(legend_labels(pe)[1:2], c("Estimated effect", "True effect"))
  expect_true("Placebo units" %in% legend_labels(
    plot(placebo_test(fit, progress = FALSE))))
  loo <- placebo_test(fit, type = "leave_one_out", progress = FALSE)
  expect_equal(legend_labels(plot(loo)),
               c("Original synthetic treated",
                 "Leave-one-out synthetic treated"))
})

test_that("labels, axis titles and legend positions can be changed", {
  fit <- small_fit()
  p <- plot(fit, treated_label = "observed target unit",
            synthetic_label = "synthetic control target unit",
            donor_label = "donor pool units", xlab = "Period", ylab = "Y",
            ylim = c(0, 50))
  expect_equal(p$labels$x, "Period")
  expect_equal(p$labels$y, "Y")
  expect_equal(legend_labels(p), c("observed target unit",
                                   "synthetic control target unit",
                                   "donor pool units"))
  for (pos in list("bottomright", "bottomleft", "topright", "topleft",
                   "bottom", "right", "none", c(0.5, 0.5))) {
    expect_s3_class(plot(fit, legend_position = pos), "ggplot")
  }
  expect_error(plot(fit, legend_position = "middle"), "legend_position")
  expect_error(plot(fit, legend_position = c(2, 0)), "legend_position")
  expect_error(plot(fit, treated_label = "a", synthetic_label = "a"),
               "distinct")
  expect_error(plot(fit, counterfactual = 1:2), "finite values")
  expect_error(plot(fit, colours = c(nope = "red")), "colours")
  p2 <- plot(fit, colours = c(treated = "darkgreen"))
  expect_equal(unname(ggplot2::ggplot_build(p2)$data[[3]]$colour[1]),
               "darkgreen")
  pl <- plot(placebo_test(small_fit(), type = "leave_one_out",
                          progress = FALSE),
             xlab = "Month", ylab = "Causal effect")
  expect_equal(pl$labels$y, "Causal effect")
  expect_equal(pl$labels$x, "Month")
})

test_that("eta heat map colours and legend title can be changed", {
  es <- eta_sensitivity(small_fit(), eta_Z = c(0.1, 0.3), eta_X = c(0.1, 0.3),
                        progress = FALSE)
  tiles <- function(p) ggplot2::ggplot_build(p)$data[[1]]
  d0 <- tiles(plot(es))
  d1 <- tiles(plot(es, colours = c("white", "red")))
  if (length(unique(d0$fill)) > 1L) {
    expect_false(identical(d0$fill, d1$fill))
    expect_true(all(c("#FFFFFF", "#FF0000") %in% toupper(d1$fill)))
  }
  p <- plot(es, colours = c("white", "orange", "red"),
            legend_title = "Causal effect", text_colour = "grey20")
  expect_equal(p$scales$get_scales("fill")$name, "Causal effect")
  expect_equal(unique(ggplot2::ggplot_build(p)$data[[2]]$colour), "grey20")
  expect_error(plot(es, colours = "red"), "at least two")
  lim <- range(es$estimate) + c(-1, 1)
  d2 <- tiles(plot(es, limits = lim))
  expect_true(all(!is.na(d2$fill)))
  expect_equal(plot(es, limits = lim)$scales$get_scales("fill")$limits, lim)
  expect_error(plot(es, limits = c(2, 1)), "limits")
  labs3 <- ggplot2::ggplot_build(plot(es, digits = 3))$data[[2]]$label
  expect_equal(labs3, formatC(es$estimate, format = "f", digits = 3))
  expect_true(all(grepl("\\.[0-9]{3}$", labs3)))
  expect_error(plot(es, digits = -1), "digits")
})
