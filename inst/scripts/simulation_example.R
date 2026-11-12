##############################################################################
# causalpanelfusion: worked example on simulated data
#
# Follows main_example.R: one data set from the paper's factor model
# (31 units, T = 20 reference periods, S = 5 target periods, seed 45), taken
# through every step of the analysis and checked against the true effect,
# which is known in a simulation.
#
# Set `quick <- FALSE` for the paper's settings (500 bootstrap refits with the
# full budget grid for every refit). The quick settings run in a few minutes.
##############################################################################

library(causalpanelfusion)
library(ggplot2)

quick <- TRUE

settings <- if (quick) {
  list(B = 100, refit_step = 0.05)     # coarser budget grid for the refits
} else {
  list(B = 500, refit_step = NULL)     # NULL = full 0.01 grid
}

# Figures are shown on screen and saved here
out_dir <- file.path(getwd(), "simulation_example")
dir.create(out_dir, showWarnings = FALSE)
save_fig <- function(p, name, width = 8, height = 6) {
  print(p)
  ggsave(file.path(out_dir, name), p, width = width, height = height)
}

# Refits (bootstrap, placebos) run on all cores; results are identical to a
# serial run.
future::plan(future::multisession)

# ---- 1. Simulate -----------------------------------------------------------
# Same draws as main_example.R with set.seed(45). Unit 1 ("treated") receives
# the effect in the target domain; target0 holds its untreated outcomes.
sim <- simulate_fusion_data(seed = 45)
dim(sim$reference)          # 31 x 20  (F)
dim(sim$target)             # 31 x 5   (Y)
sim$effect                  # true per-period effects
sim$psi                     # true average effect

# The same data are also returned as data frames:
head(sim$long)              # one row per unit, domain and period
head(sim$wide)              # one row per unit and period

# ---- 2. Synthetic control data fusion --------------------------------------
# Matrices can be passed directly ...
fit <- scfusion(sim$target, sim$reference, X = sim$X, Z = sim$Z,
                eta_Z = 0.1, eta_X = 0.1)

# ... or built from a data frame with fusion_data(). The simulated
# covariates are already on [0, 1], so normalize = FALSE gives the same fit.
fd <- fusion_data(sim$long, unit = "unit", time = "time", treated = "treated",
                  domain = "domain", target = "target", reference = "reference",
                  outcome = "outcome", covariates = c("cov1", "cov2", "cov3"),
                  normalize = FALSE)
fd
stopifnot(isTRUE(all.equal(scfusion(fd)$estimate, fit$estimate)))

fit
summary(fit)
coef(fit)[coef(fit) > 0]              # donors with positive weight
covariate_balance(fit)                # balance in both domains

# Target domain with the true counterfactual (cf. Y_plot_multilines.pdf)
save_fig(plot(fit, counterfactual = sim$target0[1, ],
              treated_label = "observed target unit",
              synthetic_label = "synthetic control target unit",
              counterfactual_label = "counterfactual target unit",
              donor_label = "donor pool units",
              ylab = "Outcome in the target domain", ylim = c(5, 50)),
         "sim_target.pdf")

# Reference domain, where the weights are fitted (cf. F_plot_multilines.pdf)
save_fig(plot(fit, type = "reference",
              treated_label = "observed target unit",
              synthetic_label = "synthetic control target unit",
              donor_label = "donor pool units",
              ylab = "Outcome in the reference domain", ylim = c(5, 50)),
         "sim_reference.pdf")

# ---- 3. Equi-confounding data fusion ---------------------------------------
ec <- equi_confounding(sim$target, sim$reference)
ec
save_fig(plot(ec), "sim_equi_confounding.pdf", 6, 5)

# Compare the three estimators with the truth
estimates <- data.frame(
  estimator = c("synthetic control fusion", "linear equi-confounding",
                "logarithmic equi-confounding"),
  estimate = c(fit$estimate, ec$estimate[["linear"]], ec$estimate[["log"]])
)
estimates$bias <- estimates$estimate - sim$psi
estimates

# ---- 4. Prediction intervals -----------------------------------------------
# alpha1 = alpha2 = 0.1 as in main_real_data.R (nominal level 80%).
pi <- scfusion_pi(fit, alpha1 = 0.1, alpha2 = 0.1, B = settings$B,
                  budget_step = settings$refit_step, seed = 1234,
                  parallel = TRUE)
pi

# Does each interval cover the true average effect?
cbind(pi$intervals[, c("method", "level", "lower", "upper")],
      covers_truth = pi$intervals$lower <= sim$psi &
        sim$psi <= pi$intervals$upper)

# Per-period bands, and whether they cover the true per-period effects
for (m in c("gaussian", "empirical", "subgaussian")) {
  b <- pi_bands(pi, method = m)
  cat(sprintf("%-11s bands cover the true effect in %d of %d periods\n", m,
              sum(b$effect_lo <= sim$effect & sim$effect <= b$effect_hi),
              nrow(b)))
}
pi_bands(pi, method = "gaussian")

save_fig(plot(pi, method = "gaussian", counterfactual = sim$target0[1, ],
              treated_label = "observed target unit",
              synthetic_label = "synthetic control target unit",
              counterfactual_label = "counterfactual target unit",
              donor_label = "donor pool units",
              ylab = "Outcome in the target domain",
              legend_position = "topleft"),
         "sim_PI_target.pdf")
save_fig(plot(pi, type = "effect", method = "gaussian",
              counterfactual = sim$target0[1, ], legend_position = "topleft"),
         "sim_PI_effect.pdf")

# Sensitivity to the transfer-error scale
sens <- sensitivity_pi(pi)
sens
save_fig(plot(sens), "sim_sensitivity_PI.pdf", 7, 5)

# ---- 5. Placebo tests ------------------------------------------------------
# In-space placebos (treated unit kept in the donor pools, as in Abadie et
# al. 2010; cf. placebo_F.pdf / placebo_Y.pdf)
pl <- placebo_test(fit, parallel = TRUE)
save_fig(plot(pl, domain = "reference", treated_label = "target unit",
              placebo_label = "control units", ylab = "Causal Effect"),
         "sim_placebo_reference.pdf")
save_fig(plot(pl, domain = "target", treated_label = "target unit",
              placebo_label = "control units", ylab = "Causal Effect"),
         "sim_placebo_target.pdf")

# RMSPE-ratio test (treated unit excluded from the placebo donor pools)
rt <- rmspe_test(fit, parallel = TRUE)
rt
save_fig(plot(rt), "sim_rmspe_ratio.pdf", 9, 5)

# Leave-one-donor-out (cf. Placebo_test_LODO)
loo <- placebo_test(fit, type = "leave_one_out", parallel = TRUE)
loo
save_fig(plot(loo, domain = "target",
              treated_label = "Original target unit",
              placebo_label = "Leave-one-out target unit",
              ylab = "Causal Effect"),
         "sim_LODO_target.pdf")

# ---- 6. Sensitivity to the balance slack -----------------------------------
es <- eta_sensitivity(fit, eta_Z = c(0.05, 0.1, 0.2), eta_X = c(0.05, 0.1, 0.2),
                      parallel = TRUE)
es[, c("eta_Z", "eta_X", "estimate", "NSE_F", "NSE_Z", "NSE_X", "n_donors")]
save_fig(plot(es), "sim_eta_heatmap.pdf", 7, 5)

# ---- 7. Save ---------------------------------------------------------------
saveRDS(list(sim = sim, fit = fit, equi = ec, estimates = estimates, pi = pi,
             placebo = pl, rmspe = rt, loo = loo, eta = es, sens = sens,
             settings = settings),
        file.path(out_dir, "simulation_example_results.rds"))

future::plan(future::sequential)
