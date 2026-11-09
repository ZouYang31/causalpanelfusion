# causalpanelfusion 0.2.0

A rewrite of the package around a single model object, adding the inference
methods of the paper.

## New

* `scfusion()` fits synthetic control data fusion and returns an
  `"scfusion"` object with `print()`, `summary()`, `coef()` and `plot()`
  methods. It replaces the chain `optimize_w_ipop()` / `NSE_x()` /
  `generate_b_list()` / `find_best_B()`.
* `scfusion_pi()` computes prediction intervals for the average effect:
  bootstrap of the weight-construction error with full budget re-selection,
  donor-placebo calibration of the transfer error, Approach 2
  (location-scale; Gaussian and empirical) and Approach 1 (sub-Gaussian).
  `pi_bands()` gives pointwise per-period bands and `sensitivity_pi()` scales
  the transfer error.
* `placebo_test()` (in-space and leave-one-out) and `rmspe_test()`
  (RMSPE-ratio placebo p-value).
* `covariate_balance()` reports the covariate balance table (treated unit,
  synthetic control, donor pool summaries, gaps on the reported and fitting
  scales, rank of the treated unit) in the original units of the covariates,
  optionally side by side with other fits.
* `eta_sensitivity()` refits over a grid of balance slacks.
* `equi_confounding()` returns the linear and logarithmic estimators together
  with a `plot()` method for the equi-confounding diagram.
* `simulate_fusion_data()` draws from the paper's factor model
  (`seed = 45` reproduces `main_example.R`).
* `panel_to_matrix()`, `covariate_matrix()` and `normalize_minmax()` prepare
  data; inputs are aligned by row name and the treated unit can be named.

## Changes

* The optimisation now uses the free ECOS solver by default (package
  ECOSolveR), so Gurobi is no longer required. `solver = "gurobi"` reproduces
  the original implementation (kernlab::ipop baselines and Gurobi
  quadratically constrained programs) exactly.
* Bootstrap and placebo refits can run in parallel through future.apply
  (`parallel = TRUE`); results do not depend on the backend.
* Plot functions return ggplot objects instead of printing and saving files,
  and no longer depend on global variables.

## Correspondence with the research code

| research code                                   | package                              |
|-------------------------------------------------|--------------------------------------|
| `generate_b_list()`                             | `budget_grid()`                      |
| `optimize_w_ipop()` + `NSE_x()` + `find_best_B()` | `scfusion()`                       |
| `synthetic_ett()`                               | `scfusion()$estimate`                |
| `linear_/log_equi_confounding_ett()`            | `equi_confounding()`                 |
| `4_5_prediction_interval_approach2.R`           | `scfusion_pi()` (location-scale rows)|
| `4_5_prediction_interval_approach1.R`           | `scfusion_pi()` (sub-Gaussian row)   |
| `compute_prediction_bands()`                    | `pi_bands()`                         |
| `synthetic_target_plot*()` with `bands`         | `plot(<scfusion_pi>)`                |
| `placebo_test_data()`                           | `placebo_test(type = "in_space")`    |
| `placebo_test_LODO_data()`                      | `placebo_test(type = "leave_one_out")` |
| `6_rmspe_ratio_test.R`                          | `rmspe_test()`                       |
| eta heat-map loop in `main_real_data.R`         | `eta_sensitivity()`                  |
| `6_sensitivity_PI.R`                            | `sensitivity_pi()`                   |
