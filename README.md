# causalpanelfusion

**causalpanelfusion** estimates the effect of an intervention on a single
treated unit when the outcome of interest (the *target domain*) has **no
pre-intervention period**. It fuses panel data from a related *reference
domain*, observed for the same units and unaffected by the intervention, so
that counterfactual outcomes can still be recovered.

It implements the methods of

> Yang, Z., Lee, S. H., Köhler, J. R. and Ghassami, A. *Causal Data Fusion for
> Panel Data without a Pre-Intervention Period.*
> [arXiv:2410.16391](https://arxiv.org/abs/2410.16391)

| | |
|---|---|
| **Synthetic control data fusion** | `scfusion()`: donor weights from a data-driven constrained optimisation in the reference domain, transported to the target domain |
| **Prediction intervals** | `scfusion_pi()`, `pi_bands()`: bootstrap of the weight-construction error + donor-placebo calibration of the transfer error (location-scale and sub-Gaussian constructions), average and per-period |
| **Diagnostics** | `placebo_test()` (in-space, leave-one-out), `rmspe_test()`, `eta_sensitivity()`, `sensitivity_pi()` |
| **Equi-confounding data fusion** | `equi_confounding()`: linear and logarithmic estimators |
| **Utilities** | `simulate_fusion_data()`, `panel_to_matrix()`, `covariate_matrix()`, `normalize_minmax()` |

## Installation

```r
# install.packages("remotes")
remotes::install_github("ZouYang31/causalpanelfusion", build_vignettes = TRUE)
```

No commercial solver is needed: the optimisation uses the free
[ECOS](https://github.com/embotech/ecos) solver. If you have a Gurobi
licence, `solver = "gurobi"` reproduces the original research implementation
exactly.

## Example

Inputs are four matrices with one row per unit: target outcomes (`target`,
the paper's *Y*), reference outcomes (`reference`, *F*), target covariates
`X` and reference covariates `Z`.

```r
library(causalpanelfusion)

sim <- simulate_fusion_data(seed = 45)   # 31 units, 20 reference / 5 target periods
fit <- scfusion(sim$target, sim$reference, X = sim$X, Z = sim$Z)
fit
#> Synthetic control data fusion
#>   treated unit: treated | donors: 30 | target periods: 5 | reference periods: 20
#>   estimated average effect (psi_hat): 3.5099
#>   donors with positive weight: 6
#> donor_17 donor_05 donor_16 donor_06 donor_18 donor_26
#>   0.2352   0.2234   0.2171   0.1797   0.0942   0.0504
sim$psi   # true effect
#> [1] 3.719489
```

`summary(fit)` adds per-period effects, the selected budget, fit criteria and
covariate balance; `plot(fit)` draws the trajectories.

### Prediction intervals

```r
pi <- scfusion_pi(fit, B = 200, budget_step = 0.05, seed = 1)
pi
#>                      method   approach level estimate   lower   upper   width
#>              in-sample only diagnostic   95%   3.5099  2.3122  4.9424  2.6302
#>   location-scale (Gaussian) Approach 2   90%   3.5099 -3.7708 11.0253 14.7961
#>  location-scale (empirical) Approach 2   90%   3.5099 -1.2463 10.3690 11.6152
#>          sub-Gaussian bound Approach 1   90%   3.5099 -6.1179 13.3724 19.4903

pi_bands(pi, method = "gaussian")   # pointwise per-period bands
plot(pi)
```

<img src="man/figures/README-target.png" width="640" />

Each bootstrap replicate re-runs the whole algorithm, including the budget
search, so this is the slow step (about 5 s per refit with the default grid
and 30 donors). For reported results use the default `B = 500` and run the
refits in parallel; the results are identical to a serial run:

```r
future::plan(future::multisession)
pi <- scfusion_pi(fit, parallel = TRUE, seed = 1)
```

`budget_step = 0.05` (refits only) or `reselect_budget = FALSE` give faster
approximations.

### Diagnostics and sensitivity

```r
rmspe_test(fit)               # RMSPE-ratio placebo test
#> Treated unit (treated) rank: 4 of 31 -> exact placebo p-value = 0.129

plot(placebo_test(fit))                          # in-space placebo gaps
placebo_test(fit, type = "leave_one_out")        # drop each donor in turn
eta_sensitivity(fit)                             # balance-slack grid
sensitivity_pi(pi)                               # scale the transfer error

equi_confounding(sim$target, sim$reference)
#> Equi-confounding data fusion
#>   linear:     psi_hat = 1.6277
#>   logarithmic: psi_hat = 2.2327
```

See `vignette("causalpanelfusion")` for a full walk-through, including how
to prepare your own long-format data.

## Interpreting the intervals

The intervals are feasible plug-in versions of the paper's oracle intervals:
the transfer scale is estimated from a small number of donor placebos and its
estimation error is not propagated, and the identifying assumptions are
maintained rather than tested. With `bias_allowance = 0` (default) they
target the effect plus the conditional transfer bias; supply a valid bound
through `bias_allowance` to target the effect itself.

## Citation

```r
citation("causalpanelfusion")
```
