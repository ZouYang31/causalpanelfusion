#' causalpanelfusion: Causal Data Fusion for Panel Data without a
#' Pre-Intervention Period
#'
#' Estimates the effect of an intervention on a single treated unit when no
#' pre-intervention data exist in the target domain, by fusing panel data
#' from a related reference domain observed for the same units.
#'
#' Main functions:
#' \itemize{
#'   \item [fusion_data()]: build the inputs from long or wide data frames.
#'   \item [scfusion()]: synthetic control data fusion (data-driven
#'     constrained weights); [covariate_balance()] for the balance table.
#'   \item [scfusion_pi()], [pi_bands()]: prediction intervals for the
#'     average effect and pointwise per-period bands.
#'   \item [placebo_test()], [rmspe_test()]: in-space placebo,
#'     leave-one-out and RMSPE-ratio diagnostics.
#'   \item [eta_sensitivity()], [sensitivity_pi()]: sensitivity analyses.
#'   \item [equi_confounding()]: linear and logarithmic equi-confounding
#'     estimators.
#'   \item [simulate_fusion_data()]: data from the paper's factor model.
#' }
#'
#' @references Yang, Z., Lee, S. H., Kohler, J. R. and Ghassami, A. Causal
#'   data fusion for panel data without a pre-intervention period.
#'   arXiv:2410.16391.
#' @keywords internal
#' @importFrom ggplot2 .data
"_PACKAGE"
