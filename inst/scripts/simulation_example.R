##############################################################################
# causalpanelfusion: worked example on simulated data
#
# Part A follows main_example.R: one data set from the paper's factor model
#         (31 units, T = 20 reference periods, S = 5 target periods,
#         seed 45), taken through every step of the analysis and checked
#         against the true effect, which is known in a simulation.
# Part B follows main_simulation.R: a Monte Carlo study of the bias of the
#         synthetic-control and equi-confounding estimators as the number of
#         reference periods T grows (design held fixed, noise redrawn).
#
# Set `quick <- FALSE` for the paper's settings (500 bootstrap refits, full
# budget grid for every refit, 300 Monte Carlo replications for
# T = 10, 20, ..., 100). The quick settings run in a few minutes.
##############################################################################

library(causalpanelfusion)
library(ggplot2)

quick <- TRUE

settings <- if (quick) {
  list(B = 100, refit_step = 0.05,                 # prediction intervals
       mc_T = c(10, 20, 50), mc_reps = 30, mc_step = 0.05)
} else {
  list(B = 500, refit_step = NULL,                 # NULL = full 0.01 grid
       mc_T = seq(10, 100, by = 10), mc_reps = 300, mc_step = 0.01)
}

# Legend inside the panel for the hand-made Monte Carlo plots (ggplot2 3.5
# changed how inside legends are specified)
inside_legend <- function(pos = c(0.98, 0.98)) {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    theme(legend.position = "inside", legend.position.inside = pos,
          legend.justification = pos)
  } else {
    theme(legend.position = pos, legend.justification = pos)
  }
}

# Figures are shown on screen and saved here
out_dir <- file.path(getwd(), "simulation_example")
dir.create(out_dir, showWarnings = FALSE)
save_fig <- function(p, name, width = 8, height = 6) {
  print(p)
  ggsave(file.path(out_dir, name), p, width = width, height = height)
}

# Refits (bootstrap, placebos, Monte Carlo) run on all cores; results are
# identical to a serial run.
future::plan(future::multisession)


# ============================================================================
# Part A. One simulated data set
# ============================================================================

# ---- A1. Simulate ------------------------------------------------------------
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

# ---- A2. Synthetic control data fusion ----------------------------------------
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

# ---- A3. Equi-confounding data fusion ------------------------------------------
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

# ---- A4. Prediction intervals --------------------------------------------------
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

# ---- A5. Placebo tests ------------------------------------------------------------
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

# ---- A6. Sensitivity to the balance slack -------------------------------------------
es <- eta_sensitivity(fit, eta_Z = c(0.05, 0.1, 0.2), eta_X = c(0.05, 0.1, 0.2),
                      parallel = TRUE)
es[, c("eta_Z", "eta_X", "estimate", "NSE_F", "NSE_Z", "NSE_X", "n_donors")]
save_fig(plot(es), "sim_eta_heatmap.pdf", 7, 5)


# ============================================================================
# Part B. Monte Carlo study (cf. main_simulation.R)
# ============================================================================
# The design (latent confounder mu, covariates Z and X, time effects and
# loadings, true effects) is drawn once and held fixed; each replication
# redraws the idiosyncratic noise in both domains and refits. The first T
# of max(T) reference periods are used. Bias = estimate - truth for all three
# estimators.

design <- simulate_fusion_data(n_ref = max(settings$mc_T), seed = 45)

one_replication <- function(T, design, budget_step) {
  N <- nrow(design$Z)
  S <- ncol(design$target)
  p <- design$params
  ones <- rep(1, N)
  Fm <- outer(ones, p$rho[seq_len(T)]) +
    design$Z %*% t(p$phi[seq_len(T), , drop = FALSE]) +
    design$mu %*% t(p$theta[seq_len(T), , drop = FALSE]) +
    matrix(rnorm(N * T, 0, sqrt(2)), N, T)
  Y <- outer(ones, p$varrho) + design$X %*% t(p$varphi) +
    design$mu %*% t(p$vartheta) + matrix(rnorm(N * S, 0, sqrt(1 / 2)), N, S)
  Y[1, ] <- Y[1, ] + design$effect          # treated unit receives the effect

  fit <- scfusion(Y, Fm, X = design$X, Z = design$Z, budget_step = budget_step)
  ec <- equi_confounding(Y, Fm)
  w <- fit$weights
  mu_gap <- design$mu[1, ] - colSums(design$mu[-1, , drop = FALSE] * w)
  c(T = T,
    bias_sc = fit$estimate - design$psi,
    bias_linear = ec$estimate[["linear"]] - design$psi,
    bias_log = ec$estimate[["log"]] - design$psi,
    nse_mu = sum(mu_gap^2) / length(mu_gap),
    nse_Z = fit$nse[["Z"]], nse_X = fit$nse[["X"]], nse_F = fit$nse[["F"]])
}

set.seed(2026)
jobs <- rep(settings$mc_T, each = settings$mc_reps)
mc <- future.apply::future_lapply(jobs, one_replication, design = design,
                                  budget_step = settings$mc_step,
                                  future.seed = TRUE)
mc <- as.data.frame(do.call(rbind, mc))

# Mean bias, its Monte Carlo standard error, and RMSE by T
summarise_mc <- function(v) {
  c(mean = mean(v), se = stats::sd(v) / sqrt(length(v)),
    rmse = sqrt(mean(v^2)))
}
mc_summary <- do.call(rbind, lapply(split(mc, mc$T), function(d) {
  data.frame(
    T = d$T[1],
    estimator = c("synthetic control", "linear equi-conf.", "log equi-conf."),
    t(sapply(d[c("bias_sc", "bias_linear", "bias_log")], summarise_mc)),
    row.names = NULL)
}))
mc_summary

# Box plots of the bias by T (cf. bias_plot())
bias_long <- data.frame(
  T = factor(rep(mc$T, 3), levels = settings$mc_T),
  estimator = factor(rep(c("synthetic control", "linear equi-confounding",
                           "log equi-confounding"), each = nrow(mc)),
                     levels = c("synthetic control", "linear equi-confounding",
                                "log equi-confounding")),
  bias = c(mc$bias_sc, mc$bias_linear, mc$bias_log)
)
p_bias <- ggplot(bias_long, aes(x = T, y = bias, fill = estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("grey70", "#3B75AF", "#C0504D"), name = NULL) +
  labs(x = "Number of reference periods T", y = "Bias (estimate - truth)") +
  theme_minimal(base_size = 14) +
  inside_legend() +
  theme(legend.background = element_rect(fill = alpha("grey", 0.5),
                                         colour = NA))
save_fig(p_bias, "sim_mc_bias.pdf")

# Box plots of the normalised squared errors of mu, Z and X (cf. NSE_plot_ZX)
nse_long <- data.frame(
  T = factor(rep(mc$T, 3), levels = settings$mc_T),
  quantity = factor(rep(c("NSE(mu)", "NSE(Z)", "NSE(X)"), each = nrow(mc)),
                    levels = c("NSE(mu)", "NSE(Z)", "NSE(X)")),
  value = c(mc$nse_mu, mc$nse_Z, mc$nse_X)
)
p_nse <- ggplot(nse_long, aes(x = T, y = value, fill = quantity)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#59A14F", "#9C6ADE", "#F28E2B"), name = NULL,
                    labels = c(expression(NSE(mu[1])), expression(NSE(Z[1])),
                               expression(NSE(X[1])))) +
  labs(x = "Number of reference periods T", y = "NSE") +
  theme_minimal(base_size = 14) +
  inside_legend() +
  theme(legend.background = element_rect(fill = alpha("grey", 0.5),
                                         colour = NA))
save_fig(p_nse, "sim_mc_nse.pdf")

# Keep everything for later
saveRDS(list(sim = sim, fit = fit, equi = ec, estimates = estimates, pi = pi,
             placebo = pl, rmspe = rt, loo = loo, eta = es, sens = sens,
             mc = mc, mc_summary = mc_summary, settings = settings),
        file.path(out_dir, "simulation_example_results.rds"))

future::plan(future::sequential)
