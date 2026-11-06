# Small problems keep the test suite fast.
small_sim <- function(seed = 11, n_units = 12, n_ref = 10, n_target = 4) {
  simulate_fusion_data(n_units = n_units, n_ref = n_ref, n_target = n_target,
                       seed = seed)
}

small_fit <- function(sim = small_sim(), step = 0.1, ...) {
  scfusion(sim$target, sim$reference, X = sim$X, Z = sim$Z,
           budget_step = step, ...)
}
