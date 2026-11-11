# Plot methods. All return ggplot objects; every coloured element has a
# legend entry, the legend sits inside the panel by default
# (`legend_position`), and axis labels can be changed with `xlab`/`ylab`.

.col <- list(
  treated = "black", synthetic = "black", donor = "grey60",
  counterfactual = "red", insample = "#1f4ee0", total = "#E06666",
  zero = "blue", sensitivity = "#7B2D26"
)

#' Plot a synthetic control data-fusion fit
#'
#' Draws the treated unit, its synthetic control and the donor pool in the
#' target or reference domain, or the per-period gaps. Prediction bands from
#' [scfusion_pi()] and, in simulations, the true counterfactual trajectory can
#' be added; each gets its own legend entry.
#'
#' @param x An `"scfusion"` object.
#' @param type `"target"` (default), `"reference"` or `"gaps"` (treated minus
#'   synthetic in both domains).
#' @param bands Optional bands to overlay on the target plot: an
#'   `"scfusion_pi"` object (Gaussian bands are used) or a data frame from
#'   [pi_bands()].
#' @param band_style `"errorbar"` (default) or `"ribbon"`.
#' @param counterfactual Optional numeric vector with the treated unit's true
#'   untreated target-domain outcomes (known in simulations, e.g.
#'   `sim$target0[1, ]`), drawn as a red line.
#' @param show_donors Draw the donor trajectories in grey.
#' @param treated_label,synthetic_label,donor_label,counterfactual_label
#'   Legend labels. Defaults: the treated unit's name, "Synthetic <name>",
#'   "Donor pool units" and "Counterfactual <name>".
#' @param band_labels Optional labels for the two bands,
#'   `c(insample = , total = )`. The blue `insample` band covers the
#'   weight-construction error only (level \eqn{1 - \alpha_1}); the red
#'   `total` band is the prediction interval, which adds the transfer error
#'   on both sides (level \eqn{1 - \alpha_1 - \alpha_2}). Defaults:
#'   "Weight construction error (95%)" and "Prediction interval (90%):
#'   weight construction + transfer error" (levels follow the fit).
#' @param colours Optional named vector changing any of the colours
#'   `treated` ("black"), `synthetic` ("black"), `counterfactual` ("red"),
#'   `donor` ("grey60"), `insample` ("#1f4ee0") and `total` ("#E06666"),
#'   e.g. `colours = c(total = "orange")`.
#' @param xlab,ylab Axis labels.
#' @param ylim Optional y-axis limits (zooms without dropping data).
#' @param legend_position Where to put the legend: inside the panel at
#'   `"bottomright"` (default), `"bottomleft"`, `"topright"` or `"topleft"`;
#'   a numeric `c(x, y)` in `[0, 1]` for any point inside the panel; or
#'   outside with `"bottom"`, `"top"`, `"left"`, `"right"`; or `"none"`.
#' @param ... Unused.
#' @return A `ggplot` object (add layers or themes with `+`).
#' @examples
#' sim <- simulate_fusion_data(n_units = 16, n_ref = 12, n_target = 4,
#'                             seed = 1)
#' fit <- scfusion(sim$target, sim$reference, X = sim$X, Z = sim$Z,
#'                 budget_step = 0.1)
#' plot(fit)
#' plot(fit, counterfactual = sim$target0[1, ], treated_label = "target unit",
#'      xlab = "Period", ylab = "Outcome", legend_position = "topleft")
#' plot(fit, type = "reference")
#' @export
plot.scfusion <- function(x, type = c("target", "reference", "gaps"),
                          bands = NULL, band_style = c("errorbar", "ribbon"),
                          counterfactual = NULL, show_donors = TRUE,
                          treated_label = NULL, synthetic_label = NULL,
                          donor_label = "Donor pool units",
                          counterfactual_label = NULL, band_labels = NULL,
                          colours = NULL,
                          xlab = "Time", ylab = NULL, ylim = NULL,
                          legend_position = "bottomright", ...) {
  type <- match.arg(type)
  band_style <- match.arg(band_style)
  col <- .resolve_colours(colours)
  tr_name <- x$data$unit_names[x$treated]
  if (type == "gaps") return(.plot_gaps(x, xlab, ylab))

  if (is.null(treated_label)) treated_label <- tr_name
  if (is.null(synthetic_label)) {
    synthetic_label <- paste("Synthetic", treated_label)
  }
  if (is.null(counterfactual_label)) {
    counterfactual_label <- paste("Counterfactual", treated_label)
  }
  M <- if (type == "target") x$data$target else x$data$reference
  syn <- if (type == "target") x$synthetic_target else x$synthetic_reference
  n_t <- ncol(M)
  if (is.null(ylab)) {
    ylab <- if (is.null(x$data$axis_labels)) {
      sprintf("Outcome (%s domain)", type)
    } else {
      sprintf("%s (%s domain)", x$data$axis_labels[[type]], type)
    }
  }

  lines <- list(
    .series_df(seq_len(n_t), M[x$treated, ], treated_label),
    .series_df(seq_len(n_t), syn, synthetic_label)
  )
  styles <- list(c(treated_label, col$treated, "solid"),
                 c(synthetic_label, col$synthetic, "dashed"))
  if (!is.null(counterfactual)) {
    if (type != "target") {
      warning("`counterfactual` is only drawn for type = \"target\".",
              call. = FALSE)
    } else {
      counterfactual <- as.numeric(counterfactual)
      if (length(counterfactual) != n_t || any(!is.finite(counterfactual))) {
        stop(sprintf("`counterfactual` must have %d finite values.", n_t),
             call. = FALSE)
      }
      lines[[3]] <- .series_df(seq_len(n_t), counterfactual,
                               counterfactual_label)
      styles[[3]] <- c(counterfactual_label, col$counterfactual, "solid")
    }
  }
  donors <- NULL
  if (show_donors) {
    donors <- data.frame(
      time = rep(seq_len(n_t), each = length(x$donors)),
      value = as.vector(M[x$donors, , drop = FALSE]),
      unit = rep(x$data$unit_names[x$donors], n_t),
      series = donor_label, stringsAsFactors = FALSE
    )
    styles[[length(styles) + 1L]] <- c(donor_label, col$donor, "solid")
  }

  bd <- NULL
  if (!is.null(bands) && type == "target") {
    bd <- .as_bands(bands)
    if (is.null(band_labels)) {
      band_labels <- c(
        insample = sprintf("Weight construction error (%s%%)",
                           format(100 * (1 - attr(bd, "alpha")[["alpha1"]]))),
        total = sprintf(
          "Prediction interval (%s%%):\nweight construction + transfer error",
          format(100 * attr(bd, "level"))))
    }
    if (length(band_labels) != 2L) {
      stop("`band_labels` must have two elements (insample, total).",
           call. = FALSE)
    }
    if (is.null(names(band_labels))) {
      names(band_labels) <- c("insample", "total")
    }
    if (band_style == "errorbar") {
      styles[[length(styles) + 1L]] <- c(band_labels[["insample"]],
                                         col$insample, "solid")
      styles[[length(styles) + 1L]] <- c(band_labels[["total"]],
                                         col$total, "solid")
    }
  }
  sty <- .style_table(styles)

  p <- ggplot2::ggplot()
  if (!is.null(bd) && band_style == "ribbon") {
    p <- p + .band_ribbons(bd, band_labels, col)
  }
  if (!is.null(donors)) {
    p <- p + ggplot2::geom_line(
      data = donors,
      ggplot2::aes(x = .data$time, y = .data$value, group = .data$unit,
                   colour = .data$series, linetype = .data$series),
      linewidth = 0.8, alpha = 0.5)
  }
  # draw synthetic, counterfactual, then the treated unit on top
  for (k in c(2L, if (length(lines) > 2L) 3L, 1L)) {
    p <- p + ggplot2::geom_line(
      data = lines[[k]],
      ggplot2::aes(x = .data$time, y = .data$value, colour = .data$series,
                   linetype = .data$series),
      linewidth = 0.8)
  }
  if (!is.null(bd) && band_style == "errorbar") {
    p <- p + .band_errorbars(bd, band_labels)
  }
  p <- p +
    ggplot2::scale_colour_manual(values = sty$colour, breaks = sty$label,
                                 name = NULL) +
    ggplot2::scale_linetype_manual(values = sty$linetype, breaks = sty$label,
                                   name = NULL) +
    ggplot2::scale_x_continuous(breaks = .int_breaks(n_t)) +
    ggplot2::labs(x = xlab, y = ylab) +
    .cpf_theme(legend_position)
  if (!is.null(ylim)) p <- p + ggplot2::coord_cartesian(ylim = ylim)
  p
}

.series_df <- function(time, value, label) {
  data.frame(time = time, value = as.numeric(value), series = label,
             stringsAsFactors = FALSE)
}

# Turn a list of c(label, colour, linetype) into named value vectors.
.style_table <- function(styles) {
  labels <- vapply(styles, `[[`, "", 1L)
  if (anyDuplicated(labels)) {
    stop(sprintf("Legend labels must be distinct; duplicated: %s.",
                 paste(unique(labels[duplicated(labels)]), collapse = ", ")),
         call. = FALSE)
  }
  list(label = labels,
       colour = stats::setNames(vapply(styles, `[[`, "", 2L), labels),
       linetype = stats::setNames(vapply(styles, `[[`, "", 3L), labels))
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

.band_errorbars <- function(bd, labels) {
  tot <- data.frame(time = bd$time, lo = bd$total_lo, hi = bd$total_hi,
                    series = labels[["total"]], stringsAsFactors = FALSE)
  ins <- data.frame(time = bd$time, lo = bd$insample_lo, hi = bd$insample_hi,
                    series = labels[["insample"]], stringsAsFactors = FALSE)
  list(
    ggplot2::geom_errorbar(
      data = tot, ggplot2::aes(x = .data$time, ymin = .data$lo,
                               ymax = .data$hi, colour = .data$series,
                               linetype = .data$series),
      inherit.aes = FALSE, width = 0.39, linewidth = 0.7),
    ggplot2::geom_errorbar(
      data = ins, ggplot2::aes(x = .data$time, ymin = .data$lo,
                               ymax = .data$hi, colour = .data$series,
                               linetype = .data$series),
      inherit.aes = FALSE, width = 0.3, linewidth = 0.8)
  )
}

.band_ribbons <- function(bd, labels, col) {
  df <- rbind(
    data.frame(time = bd$time, lo = bd$total_lo, hi = bd$total_hi,
               series = labels[["total"]], stringsAsFactors = FALSE),
    data.frame(time = bd$time, lo = bd$insample_lo, hi = bd$insample_hi,
               series = labels[["insample"]], stringsAsFactors = FALSE)
  )
  fills <- stats::setNames(c(col$total, col$insample),
                           c(labels[["total"]], labels[["insample"]]))
  list(
    ggplot2::geom_ribbon(
      data = df, ggplot2::aes(x = .data$time, ymin = .data$lo,
                              ymax = .data$hi, fill = .data$series),
      inherit.aes = FALSE, alpha = 0.25),
    ggplot2::scale_fill_manual(values = fills,
                               breaks = c(labels[["insample"]],
                                          labels[["total"]]),
                               name = NULL)
  )
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
                        colour = .col$zero, linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::facet_wrap(~domain, scales = "free") +
    ggplot2::labs(x = xlab, y = if (is.null(ylab)) "Treated - synthetic" else
      ylab) +
    .cpf_theme("none")
}

#' Plot prediction intervals
#'
#' @param x An `"scfusion_pi"` object.
#' @param method Band construction passed to [pi_bands()].
#' @param type `"counterfactual"` (default): the target-domain trajectories
#'   with the weight-construction interval and the prediction interval
#'   (weight construction + transfer error) on the synthetic control;
#'   `"effect"`: per-period effects with their prediction intervals.
#' @param counterfactual Optional true untreated outcomes of the treated
#'   unit (simulations). With `type = "counterfactual"` it is drawn as a red
#'   line; with `type = "effect"` the true effects are shown as red
#'   triangles.
#' @param effect_label,true_label,band_label Legend labels for
#'   `type = "effect"`.
#' @param colours Optional named vector of colours, as in [plot.scfusion()]
#'   (`type = "effect"` uses `treated` for the estimates, `counterfactual`
#'   for the true effects and `total` for the intervals).
#' @param xlab,ylab,ylim,legend_position As in [plot.scfusion()].
#' @param ... Further arguments passed to [plot.scfusion()] for
#'   `type = "counterfactual"` (e.g. `treated_label`, `donor_label`,
#'   `band_labels`, `band_style`).
#' @return A `ggplot` object.
#' @export
plot.scfusion_pi <- function(x, method = c("gaussian", "empirical",
                                           "subgaussian"),
                             type = c("counterfactual", "effect"),
                             counterfactual = NULL,
                             effect_label = "Estimated effect",
                             true_label = "True effect", band_label = NULL,
                             colours = NULL,
                             xlab = "Time", ylab = NULL, ylim = NULL,
                             legend_position = "bottomright", ...) {
  method <- match.arg(method)
  type <- match.arg(type)
  bd <- pi_bands(x, method = method)
  if (type == "counterfactual") {
    return(plot.scfusion(x$fit, bands = bd, counterfactual = counterfactual,
                         colours = colours, xlab = xlab, ylab = ylab,
                         ylim = ylim, legend_position = legend_position, ...))
  }
  if (is.null(band_label)) {
    band_label <- sprintf("Prediction interval (%s%%)",
                          format(100 * attr(bd, "level")))
  }
  col <- .resolve_colours(colours)
  styles <- list(c(effect_label, col$treated, "blank"))
  shapes <- 16
  est <- .series_df(bd$time, bd$effect, effect_label)
  truth <- NULL
  if (!is.null(counterfactual)) {
    counterfactual <- as.numeric(counterfactual)
    if (length(counterfactual) != nrow(bd)) {
      stop(sprintf("`counterfactual` must have %d values.", nrow(bd)),
           call. = FALSE)
    }
    truth <- .series_df(bd$time, bd$treated - counterfactual, true_label)
    styles[[2]] <- c(true_label, col$counterfactual, "blank")
    shapes <- c(shapes, 17)
  }
  styles[[length(styles) + 1L]] <- c(band_label, col$total, "solid")
  shapes <- c(shapes, NA)
  sty <- .style_table(styles)
  bars <- data.frame(time = bd$time, lo = bd$effect_lo, hi = bd$effect_hi,
                     series = band_label, stringsAsFactors = FALSE)

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey40") +
    ggplot2::geom_errorbar(
      data = bars, ggplot2::aes(x = .data$time, ymin = .data$lo,
                                ymax = .data$hi, colour = .data$series),
      width = 0.3, linewidth = 0.7) +
    ggplot2::geom_point(
      data = est, ggplot2::aes(x = .data$time, y = .data$value,
                               colour = .data$series), size = 2.2)
  if (!is.null(truth)) {
    p <- p + ggplot2::geom_point(
      data = truth, ggplot2::aes(x = .data$time, y = .data$value,
                                 colour = .data$series),
      shape = 17, size = 2.4)
  }
  p <- p +
    ggplot2::scale_colour_manual(values = sty$colour, breaks = sty$label,
                                 name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(shape = shapes, linetype = unname(sty$linetype)))) +
    ggplot2::scale_x_continuous(breaks = .int_breaks(nrow(bd))) +
    ggplot2::labs(x = xlab,
                  y = if (is.null(ylab)) "Effect on the treated unit" else
                    ylab) +
    .cpf_theme(legend_position)
  if (!is.null(ylim)) p <- p + ggplot2::coord_cartesian(ylim = ylim)
  p
}

#' Plot placebo-test gaps
#'
#' For in-space placebos, the treated unit's gap (black) is drawn against
#' the placebo units' gaps (grey). For the leave-one-out test, the gap of the
#' original synthetic control (thick black) is drawn against the gaps
#' obtained when each positive-weight donor is left out in turn (grey).
#'
#' @param x An `"scfusion_placebo"` object.
#' @param domain `"target"` (default) or `"reference"`.
#' @param treated_label,placebo_label Legend labels. Defaults: the treated
#'   unit's name and "Placebo units" (in-space), or "Original synthetic
#'   <name>" and "Leave-one-out synthetic <name>" (leave-one-out).
#' @param xlab,ylab,ylim,legend_position As in [plot.scfusion()].
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
plot.scfusion_placebo <- function(x, domain = c("target", "reference"),
                                  treated_label = NULL, placebo_label = NULL,
                                  xlab = "Time", ylab = NULL, ylim = NULL,
                                  legend_position = "bottomright", ...) {
  domain <- match.arg(domain)
  loo <- x$type == "leave_one_out"
  if (is.null(treated_label)) {
    treated_label <- if (loo) paste("Original synthetic", x$treated) else
      x$treated
  }
  if (is.null(placebo_label)) {
    placebo_label <- if (loo) paste("Leave-one-out synthetic", x$treated) else
      "Placebo units"
  }
  if (is.null(ylab)) {
    ylab <- if (domain == "target") "Effect (target domain)" else
      "Gap (reference domain)"
  }
  df <- x$gaps[x$gaps$domain == domain, , drop = FALSE]
  df$series <- ifelse(df$role == "treated", treated_label, placebo_label)
  sty <- .style_table(list(c(treated_label, .col$treated, "solid"),
                           c(placebo_label, .col$donor, "solid")))
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = .col$zero, linewidth = 0.5) +
    ggplot2::geom_line(
      data = df[df$role != "treated", ],
      ggplot2::aes(x = .data$time, y = .data$gap, group = .data$unit,
                   colour = .data$series, linetype = .data$series),
      linewidth = 0.8, alpha = 0.6) +
    ggplot2::geom_line(
      data = df[df$role == "treated", ],
      ggplot2::aes(x = .data$time, y = .data$gap, group = .data$unit,
                   colour = .data$series, linetype = .data$series),
      linewidth = if (loo) 1.2 else 0.8) +
    ggplot2::scale_colour_manual(values = sty$colour, breaks = sty$label,
                                 name = NULL) +
    ggplot2::scale_linetype_manual(values = sty$linetype, breaks = sty$label,
                                   name = NULL) +
    ggplot2::scale_x_continuous(breaks = .int_breaks(max(df$time))) +
    ggplot2::labs(x = xlab, y = ylab) +
    .cpf_theme(legend_position)
  if (!is.null(ylim)) p <- p + ggplot2::coord_cartesian(ylim = ylim)
  p
}

#' Plot the RMSPE ratios
#' @param x An `"scfusion_rmspe"` object.
#' @param treated_label,placebo_label Legend labels.
#' @param xlab,ylab,legend_position As in [plot.scfusion()]; the legend
#'   defaults to the top right, where the bars (sorted from largest) leave
#'   room.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
plot.scfusion_rmspe <- function(x, treated_label = NULL,
                                placebo_label = "Placebo units", xlab = NULL,
                                ylab = "RMSPE ratio (target / reference)",
                                legend_position = "topright", ...) {
  if (is.null(treated_label)) treated_label <- x$treated
  tab <- x$table
  tab$unit <- factor(tab$unit, levels = tab$unit)
  tab$series <- ifelse(tab$role == "treated", treated_label, placebo_label)
  fills <- stats::setNames(c("black", "grey70"),
                           c(treated_label, placebo_label))
  ggplot2::ggplot(tab, ggplot2::aes(x = .data$unit, y = .data$ratio,
                                    fill = .data$series)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = fills, breaks = names(fills),
                               name = NULL) +
    ggplot2::labs(x = xlab, y = ylab,
                  caption = sprintf("Exact placebo p-value = %d/%d = %s",
                                    x$rank, x$n_units, .fmt(x$p_value, 3))) +
    .cpf_theme(legend_position) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60,
                                                       hjust = 1, size = 11))
}

#' Heat map of the eta sensitivity analysis
#' @param x An `"scfusion_eta"` object.
#' @param value Column to display (default `"estimate"`).
#' @param colours Colours of the scale, from the smallest to the largest
#'   value: two colours for a simple gradient or more for a multi-step one,
#'   e.g. `colours = c("white", "orange", "red")`.
#' @param limits Optional range `c(low, high)` of the colour scale. By
#'   default the scale stretches over the values shown, so even negligible
#'   differences span the full range of colours; set e.g.
#'   `limits = c(0, 5)` to colour the tiles on a fixed range instead (values
#'   outside it are drawn grey).
#' @param legend_title Title of the colour bar; defaults to the name of the
#'   displayed column.
#' @param text_colour Colour of the numbers printed in the tiles (use a
#'   light colour on a dark scale).
#' @param xlab,ylab Axis labels.
#' @param legend_position As in [plot.scfusion()]; defaults to `"right"`
#'   (outside), since an inside legend would cover tiles.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
plot.scfusion_eta <- function(x, value = "estimate",
                              colours = c("#EAF2FB", "#6BAED6"),
                              limits = NULL,
                              legend_title = value, text_colour = "black",
                              xlab = expression(eta[X]),
                              ylab = expression(eta[Z]),
                              legend_position = "right", ...) {
  if (!value %in% names(x)) {
    stop(sprintf("Column '%s' not found.", value), call. = FALSE)
  }
  if (!is.character(colours) || length(colours) < 2L || anyNA(colours)) {
    stop("`colours` must be a vector of at least two colours ",
         "(smallest to largest value).", call. = FALSE)
  }
  if (!is.null(limits) &&
      (!is.numeric(limits) || length(limits) != 2L || anyNA(limits) ||
       limits[1L] >= limits[2L])) {
    stop("`limits` must be an increasing numeric c(low, high).",
         call. = FALSE)
  }
  df <- data.frame(eta_X = factor(x$eta_X), eta_Z = factor(x$eta_Z),
                   value = x[[value]])
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$eta_X, y = .data$eta_Z,
                                        fill = .data$value)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = signif(.data$value, 3)),
                       size = 4.5, colour = text_colour) +
    ggplot2::scale_fill_gradientn(colours = colours, limits = limits,
                                  name = legend_title) +
    ggplot2::labs(x = xlab, y = ylab) +
    .cpf_theme(legend_position)
  # the colour bar keeps its title (the displayed quantity)
  p + ggplot2::theme(legend.title = ggplot2::element_text(size = 13))
}

#' Plot the transfer-scale sensitivity analysis
#' @param x An `"scfusion_sensitivity"` object.
#' @param estimate_label,interval_label Legend labels.
#' @param colour Colour of the intervals.
#' @param xlab,ylab,legend_position As in [plot.scfusion()].
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
plot.scfusion_sensitivity <- function(x, estimate_label = "Estimate",
                                      interval_label = NULL,
                                      colour = "#7B2D26",
                                      xlab = expression(rho ~
                                        "(multiplier of the transfer scale)"),
                                      ylab = NULL,
                                      legend_position = "topleft", ...) {
  per <- attr(x, "period")
  if (is.null(interval_label)) {
    interval_label <- sprintf("Prediction interval (%s%%, %s)",
                              format(100 * attr(x, "level")),
                              attr(x, "method"))
  }
  if (is.null(ylab)) {
    ylab <- if (is.null(per)) "Average effect" else
      sprintf("Effect in target period %d", per)
  }
  sty <- .style_table(list(c(estimate_label, "black", "blank"),
                           c(interval_label, colour, "solid")))
  df <- as.data.frame(unclass(x))
  est <- .series_df(df$rho, df$estimate, estimate_label)
  bars <- data.frame(rho = df$rho, lo = df$lower, hi = df$upper,
                     series = interval_label, stringsAsFactors = FALSE)
  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::geom_errorbar(
      data = bars, ggplot2::aes(x = .data$rho, ymin = .data$lo,
                                ymax = .data$hi, colour = .data$series),
      width = 0.1, linewidth = 0.7) +
    ggplot2::geom_point(
      data = est, ggplot2::aes(x = .data$time, y = .data$value,
                               colour = .data$series), size = 2) +
    ggplot2::scale_colour_manual(values = sty$colour, breaks = sty$label,
                                 name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(shape = c(16, NA),
                          linetype = c("blank", "solid")))) +
    ggplot2::scale_x_continuous(breaks = df$rho) +
    ggplot2::labs(x = xlab, y = ylab) +
    .cpf_theme(legend_position)
}

#' Plot the equi-confounding diagram
#' @param x An `"equi_confounding"` object.
#' @param treated_label,donor_label Legend labels.
#' @param xlab,ylab,legend_position As in [plot.scfusion()].
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
plot.equi_confounding <- function(x, treated_label = "Treated",
                                  donor_label = "Donors", xlab = NULL,
                                  ylab = "Average outcome",
                                  legend_position = "bottomright", ...) {
  m <- x$means
  df <- data.frame(
    domain = c(0, 1, 0, 1),
    value = c(m[["F0"]], m[["Y0"]], m[["F1"]], m[["Y1"]]),
    series = c(donor_label, donor_label, treated_label, treated_label),
    label = c("bar(F)[0]", "bar(Y)[0]", "bar(F)[1]", "bar(Y)[1]"),
    stringsAsFactors = FALSE
  )
  sty <- .style_table(list(c(treated_label, "black", "dashed"),
                           c(donor_label, "black", "solid")))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$domain, y = .data$value,
                                   group = .data$series,
                                   linetype = .data$series)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), parse = TRUE,
                       vjust = -0.6, size = 5) +
    ggplot2::scale_linetype_manual(values = sty$linetype, breaks = sty$label,
                                   name = NULL) +
    ggplot2::scale_x_continuous(breaks = c(0, 1),
                                labels = c("Reference", "Target"),
                                expand = ggplot2::expansion(add = 0.2)) +
    ggplot2::labs(x = xlab, y = ylab) +
    .cpf_theme(legend_position)
}

# Default colours, optionally overridden by a named vector/list.
.resolve_colours <- function(colours) {
  col <- .col
  if (is.null(colours)) return(col)
  colours <- as.list(colours)
  bad <- setdiff(names(colours), c("treated", "synthetic", "counterfactual",
                                   "donor", "insample", "total"))
  if (is.null(names(colours)) || any(names(colours) == "") || length(bad)) {
    stop("`colours` must be named with any of: treated, synthetic, ",
         "counterfactual, donor, insample, total.", call. = FALSE)
  }
  col[names(colours)] <- lapply(colours, as.character)
  col
}

# ---- theme -----------------------------------------------------------------

.cpf_theme <- function(legend_position = "bottomright") {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans"),
      axis.title = ggplot2::element_text(size = 18),
      axis.text = ggplot2::element_text(size = 14),
      strip.text = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 13),
      legend.title = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(1, "lines"),
      legend.key.width = ggplot2::unit(2, "lines"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    .legend_theme(legend_position)
}

# Legend inside the panel (corner keyword or c(x, y)) or outside.
.legend_theme <- function(position) {
  corners <- list(bottomright = c(0.98, 0.02), bottomleft = c(0.02, 0.02),
                  topright = c(0.98, 0.98), topleft = c(0.02, 0.98))
  outside <- c("none", "bottom", "top", "left", "right")
  if (is.character(position) && length(position) == 1L) {
    if (position %in% outside) {
      return(ggplot2::theme(legend.position = position))
    }
    if (!position %in% names(corners)) {
      stop("`legend_position` must be one of \"bottomright\", ",
           "\"bottomleft\", \"topright\", \"topleft\", \"bottom\", \"top\", ",
           "\"left\", \"right\", \"none\", or a numeric c(x, y) in [0, 1].",
           call. = FALSE)
    }
    pos <- corners[[position]]
  } else if (is.numeric(position) && length(position) == 2L &&
             all(is.finite(position)) && all(position >= 0 & position <= 1)) {
    pos <- as.numeric(position)
  } else {
    stop("`legend_position` must be a position keyword or a numeric c(x, y) ",
         "in [0, 1].", call. = FALSE)
  }
  inside <- ggplot2::theme(
    legend.background = ggplot2::element_rect(
      fill = ggplot2::alpha("grey", 0.5), colour = NA),
    legend.margin = ggplot2::margin(3, 6, 3, 6),
    legend.justification = pos
  )
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    inside + ggplot2::theme(legend.position = "inside",
                            legend.position.inside = pos)
  } else {
    inside + ggplot2::theme(legend.position = pos)
  }
}

.int_breaks <- function(n) {
  if (n <= 12) seq_len(n) else unique(round(pretty(seq_len(n))))
}
