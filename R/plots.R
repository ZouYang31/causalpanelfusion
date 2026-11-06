#' Plot a synthetic control data-fusion fit
#'
#' Draws the treated unit, its synthetic control and the donors in the target
#' or reference domain, or the per-period gaps. Prediction bands from
#' [scfusion_pi()] can be overlaid on the synthetic target trajectory.
#'
#' @param x An `"scfusion"` object.
#' @param type `"target"` (default), `"reference"` or `"gaps"` (treated minus
#'   synthetic in both domains).
#' @param bands Optional bands to overlay on the target plot: an
#'   `"scfusion_pi"` object (Gaussian bands are used) or a data frame from
#'   [pi_bands()].
#' @param band_style `"errorbar"` (default) or `"ribbon"`.
#' @param show_donors Draw the donor trajectories in grey.
#' @param treated_label,donor_label Legend labels.
#' @param xlab,ylab Axis labels.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @examples
#' sim <- simulate_fusion_data(n_units = 16, n_ref = 12, n_target = 4,
#'                             seed = 1)
#' fit <- scfusion(sim$target, sim$reference, X = sim$X, Z = sim$Z,
#'                 budget_step = 0.1)
#' plot(fit)
#' plot(fit, type = "reference")
#' @export
plot.scfusion <- function(x, type = c("target", "reference", "gaps"),
                          bands = NULL, band_style = c("errorbar", "ribbon"),
                          show_donors = TRUE, treated_label = NULL,
                          donor_label = "Donors", xlab = "Time",
                          ylab = NULL, ...) {
  type <- match.arg(type)
  band_style <- match.arg(band_style)
  tr_name <- x$data$unit_names[x$treated]
  if (is.null(treated_label)) treated_label <- tr_name
  syn_label <- paste("Synthetic", treated_label)
  if (type == "gaps") return(.plot_gaps(x, xlab, ylab))

  M <- if (type == "target") x$data$target else x$data$reference
  syn <- if (type == "target") x$synthetic_target else x$synthetic_reference
  if (is.null(ylab)) {
    ylab <- if (type == "target") "Outcome (target domain)" else
      "Outcome (reference domain)"
  }
  n_t <- ncol(M)
  donors <- data.frame(
    time = rep(seq_len(n_t), each = length(x$donors)),
    value = as.vector(M[x$donors, , drop = FALSE]),
    unit = rep(x$data$unit_names[x$donors], n_t),
    series = donor_label, stringsAsFactors = FALSE
  )
  main <- rbind(
    data.frame(time = seq_len(n_t), value = as.numeric(M[x$treated, ]),
               unit = tr_name, series = treated_label,
               stringsAsFactors = FALSE),
    data.frame(time = seq_len(n_t), value = as.numeric(syn),
               unit = "synthetic", series = syn_label,
               stringsAsFactors = FALSE)
  )
  lev <- c(treated_label, syn_label, donor_label)
  colours <- stats::setNames(c("black", "black", "grey60"), lev)
  ltys <- stats::setNames(c("solid", "dashed", "solid"), lev)

  p <- ggplot2::ggplot()
  if (!is.null(bands) && type == "target") {
    bd <- .as_bands(bands)
    if (band_style == "ribbon") p <- p + .band_layers(bd, "ribbon")
  }
  if (show_donors) {
    p <- p + ggplot2::geom_line(
      data = donors,
      ggplot2::aes(x = .data$time, y = .data$value, group = .data$unit,
                   colour = .data$series, linetype = .data$series),
      linewidth = 0.6, alpha = 0.5)
  }
  p <- p +
    ggplot2::geom_line(
      data = main,
      ggplot2::aes(x = .data$time, y = .data$value, group = .data$series,
                   colour = .data$series, linetype = .data$series),
      linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = colours, breaks = lev, name = NULL) +
    ggplot2::scale_linetype_manual(values = ltys, breaks = lev, name = NULL) +
    ggplot2::scale_x_continuous(breaks = .int_breaks(n_t)) +
    ggplot2::labs(x = xlab, y = ylab) +
    .cpf_theme()
  if (!is.null(bands) && type == "target") {
    if (band_style == "errorbar") p <- p + .band_layers(bd, "errorbar")
    p <- p + ggplot2::labs(caption = sprintf(
      "Bands on the synthetic control: blue = weight-construction error (%s%%), red = total (%s%%, %s)",
      format(100 * (1 - attr(bd, "alpha")[["alpha1"]])),
      format(100 * attr(bd, "level")), attr(bd, "method")))
  }
  p
}

.as_bands <- function(bands) {
  if (inherits(bands, "scfusion_pi")) bands <- pi_bands(bands, "gaussian")
  need <- c("time", "insample_lo", "insample_hi", "total_lo", "total_hi")
  if (!is.data.frame(bands) || !all(need %in% names(bands))) {
    stop("`bands` must be an scfusion_pi object or a pi_bands() data frame.",
         call. = FALSE)
  }
  bands
}

.band_layers <- function(bd, style, col_in = "#1f4ee0", col_out = "#E06666") {
  if (style == "ribbon") {
    list(
      ggplot2::geom_ribbon(
        data = bd, ggplot2::aes(x = .data$time, ymin = .data$total_lo,
                                ymax = .data$total_hi),
        inherit.aes = FALSE, fill = col_out, alpha = 0.18),
      ggplot2::geom_ribbon(
        data = bd, ggplot2::aes(x = .data$time, ymin = .data$insample_lo,
                                ymax = .data$insample_hi),
        inherit.aes = FALSE, fill = col_in, alpha = 0.22)
    )
  } else {
    list(
      ggplot2::geom_errorbar(
        data = bd, ggplot2::aes(x = .data$time, ymin = .data$total_lo,
                                ymax = .data$total_hi),
        inherit.aes = FALSE, colour = col_out, width = 0.39,
        linewidth = 0.7),
      ggplot2::geom_errorbar(
        data = bd, ggplot2::aes(x = .data$time, ymin = .data$insample_lo,
                                ymax = .data$insample_hi),
        inherit.aes = FALSE, colour = col_in, width = 0.3, linewidth = 0.8)
    )
  }
}

.plot_gaps <- function(x, xlab, ylab) {
  df <- rbind(
    data.frame(domain = "Reference domain",
               time = seq_along(x$synthetic_reference),
               gap = as.numeric(x$data$reference[x$treated, ] -
                                  x$synthetic_reference)),
    data.frame(domain = "Target domain", time = seq_along(x$effects),
               gap = as.numeric(x$effects))
  )
  ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$gap)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "steelblue") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::facet_wrap(~domain, scales = "free") +
    ggplot2::labs(x = xlab, y = if (is.null(ylab)) "Treated - synthetic" else
      ylab) +
    .cpf_theme()
}

#' Plot prediction intervals
#'
#' @param x An `"scfusion_pi"` object.
#' @param method Band construction passed to [pi_bands()].
#' @param type `"counterfactual"` (default): the target-domain trajectories
#'   with bands on the synthetic control; `"effect"`: per-period effects with
#'   their bands.
#' @param ... Passed to [plot.scfusion()] for `type = "counterfactual"`.
#' @return A `ggplot` object.
#' @export
plot.scfusion_pi <- function(x, method = c("gaussian", "empirical",
                                           "subgaussian"),
                             type = c("counterfactual", "effect"), ...) {
  method <- match.arg(method)
  type <- match.arg(type)
  bd <- pi_bands(x, method = method)
  if (type == "counterfactual") return(plot.scfusion(x$fit, bands = bd, ...))
  ggplot2::ggplot(bd, ggplot2::aes(x = .data$time, y = .data$effect)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey40") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$effect_lo,
                                        ymax = .data$effect_hi),
                           colour = "#E06666", width = 0.3, linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = .int_breaks(nrow(bd))) +
    ggplot2::labs(x = "Time", y = "Effect on the treated unit",
                  caption = sprintf("Pointwise %s%% prediction bands (%s)",
                                    format(100 * attr(bd, "level")),
                                    method)) +
    .cpf_theme()
}

#' Plot placebo-test gaps
#'
#' @param x An `"scfusion_placebo"` object.
#' @param domain `"target"` (default) or `"reference"`.
#' @param ylim Optional y-axis limits.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
plot.scfusion_placebo <- function(x, domain = c("target", "reference"),
                                  ylim = NULL, ...) {
  domain <- match.arg(domain)
  df <- x$gaps[x$gaps$domain == domain, , drop = FALSE]
  if (x$type == "in_space") {
    lev <- c(x$treated, "Placebo units")
  } else {
    lev <- c("All donors", "One donor left out")
  }
  df$series <- ifelse(df$role == "treated", lev[1], lev[2])
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "steelblue") +
    ggplot2::geom_line(
      data = df[df$role != "treated", ],
      ggplot2::aes(x = .data$time, y = .data$gap, group = .data$unit,
                   colour = .data$series), linewidth = 0.6, alpha = 0.6) +
    ggplot2::geom_line(
      data = df[df$role == "treated", ],
      ggplot2::aes(x = .data$time, y = .data$gap, group = .data$unit,
                   colour = .data$series), linewidth = 1) +
    ggplot2::scale_colour_manual(values = stats::setNames(c("black",
                                                            "grey60"), lev),
                                 breaks = lev, name = NULL) +
    ggplot2::scale_x_continuous(breaks = .int_breaks(max(df$time))) +
    ggplot2::labs(x = "Time", y = sprintf("Gap (%s domain)", domain)) +
    .cpf_theme()
  if (!is.null(ylim)) p <- p + ggplot2::coord_cartesian(ylim = ylim)
  p
}

#' Plot the RMSPE ratios
#' @param x An `"scfusion_rmspe"` object.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
plot.scfusion_rmspe <- function(x, ...) {
  tab <- x$table
  tab$unit <- factor(tab$unit, levels = tab$unit)
  ggplot2::ggplot(tab, ggplot2::aes(x = .data$unit, y = .data$ratio,
                                    fill = .data$role)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c(treated = "black",
                                          placebo = "grey70"), name = NULL) +
    ggplot2::labs(x = NULL, y = "RMSPE ratio (target / reference)",
                  caption = sprintf("Exact placebo p-value = %d/%d = %s",
                                    x$rank, x$n_units, .fmt(x$p_value, 3))) +
    .cpf_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60,
                                                       hjust = 1))
}

#' Heat map of the eta sensitivity analysis
#' @param x An `"scfusion_eta"` object.
#' @param value Column to display (default `"estimate"`).
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
plot.scfusion_eta <- function(x, value = "estimate", ...) {
  if (!value %in% names(x)) {
    stop(sprintf("Column '%s' not found.", value), call. = FALSE)
  }
  df <- data.frame(eta_X = factor(x$eta_X), eta_Z = factor(x$eta_Z),
                   value = x[[value]])
  ggplot2::ggplot(df, ggplot2::aes(x = .data$eta_X, y = .data$eta_Z,
                                   fill = .data$value)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = signif(.data$value, 3)),
                       size = 4) +
    ggplot2::scale_fill_gradient(low = "#EAF2FB", high = "#6BAED6",
                                 name = value) +
    ggplot2::labs(x = expression(eta[X]), y = expression(eta[Z])) +
    .cpf_theme()
}

#' Plot the transfer-scale sensitivity analysis
#' @param x An `"scfusion_sensitivity"` object.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
plot.scfusion_sensitivity <- function(x, ...) {
  per <- attr(x, "period")
  ggplot2::ggplot(x, ggplot2::aes(x = .data$rho)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower,
                                        ymax = .data$upper),
                           width = 0.1, colour = "#7B2D26",
                           linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(y = .data$estimate), size = 1.8) +
    ggplot2::scale_x_continuous(breaks = x$rho) +
    ggplot2::labs(
      x = expression(rho ~ "(multiplier of the transfer scale)"),
      y = if (is.null(per)) "Average effect" else
        sprintf("Effect in target period %d", per),
      caption = sprintf(
        "%s%% intervals (%s). rho = 0: weight-construction error only; rho = 1: estimated scale.",
        format(100 * attr(x, "level")), attr(x, "method"))) +
    .cpf_theme()
}

#' Plot the equi-confounding diagram
#' @param x An `"equi_confounding"` object.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
plot.equi_confounding <- function(x, ...) {
  m <- x$means
  df <- data.frame(
    domain = c(0, 1, 0, 1),
    value = c(m[["F0"]], m[["Y0"]], m[["F1"]], m[["Y1"]]),
    group = c("Donors", "Donors", "Treated", "Treated"),
    label = c("bar(F)[0]", "bar(Y)[0]", "bar(F)[1]", "bar(Y)[1]")
  )
  ggplot2::ggplot(df, ggplot2::aes(x = .data$domain, y = .data$value,
                                   group = .data$group,
                                   linetype = .data$group)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), parse = TRUE,
                       vjust = -0.6) +
    ggplot2::scale_x_continuous(breaks = c(0, 1),
                                labels = c("Reference", "Target"),
                                expand = ggplot2::expansion(add = 0.2)) +
    ggplot2::labs(x = NULL, y = "Average outcome", linetype = NULL) +
    .cpf_theme()
}

.cpf_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.minor = ggplot2::element_blank())
}

.int_breaks <- function(n) {
  if (n <= 12) seq_len(n) else unique(round(pretty(seq_len(n))))
}
