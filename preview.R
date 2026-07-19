# ══════════════════════════════════════════════════════════════════════════════
# preview.R — interactive plot preview harness for rapid visual iteration
# ══════════════════════════════════════════════════════════════════════════════
#
# Not part of the deployed app. Source it AFTER app.R from an RStudio session:
#
#     source("app.R")
#     source("preview.R")
#
# Then render a plot straight to the Plots pane, flipping the same switches the
# app exposes plus a `trace_alpha` knob for the faint background traces:
#
#     preview_psup("glycolytic=PGK1,FBA1,TDH3; ribosomal=RPL4A,RPL19A,RPS9B",
#                  traces = TRUE, errorbars = TRUE, error_type = "SEM")
#
#     preview_psup("glycolytic=PGK1,FBA1,TDH3; ribosomal=RPL4A,RPL19A,RPS9B",
#                  dataset = "dia", plot = "temperature",
#                  traces = TRUE, trace_alpha = 0.2)
#
# Sweep several alpha values to compare (each opens in the plot history):
#
#     preview_alpha_sweep("glycolytic=PGK1,FBA1,TDH3; ribosomal=RPL4A,RPL19A,RPS9B",
#                         alphas = c(0.05, 0.1, 0.2, 0.3))
#
# `trace_alpha` sets options(psup.trace_alpha=) for the duration of the render,
# which is the exact knob app.R reads — so what you see here is what the app
# draws. Set it globally (options(psup.trace_alpha = 0.15)) to preview it live
# in the running app too. Every call returns the ggplot invisibly, so you can
# `p <- preview_psup(...); ggsave("fig.pdf", p)`.

# ── Main entry point ──────────────────────────────────────────────────────────
preview_psup <- function(ids,
                         dataset     = c("wallace", "dia"),
                         plot        = c("time", "temperature"),
                         traces      = TRUE,
                         errorbars   = FALSE,
                         error_type  = "SEM",          # "SEM" / "SD" (or se/sd)
                         trace_alpha = NULL,           # NULL = leave option as-is
                         idType      = c("gene", "orf"),  # wallace only
                         bioreps     = FALSE,          # dia only
                         species     = NULL,           # dia: NULL = all species
                         base_temps  = NULL,           # dia: NULL = all start temps
                         end_temp    = NULL,           # dia time plots: heat-shock temp
                         temp_range  = c(23, 50),      # dia temperature x-range
                         time_range  = c(0, 20),       # dia time x-range
                         print       = TRUE) {
  dataset <- match.arg(dataset)
  plot    <- match.arg(plot)
  idType  <- match.arg(idType)
  et      <- if (toupper(as.character(error_type)[1]) == "SD") "sd" else "se"

  # Apply the trace-alpha override for this render only, then restore it. The
  # value is baked into the geom when the plot is built, so restoring is safe.
  if (!is.null(trace_alpha)) {
    old <- getOption("psup.trace_alpha")
    options(psup.trace_alpha = trace_alpha)
    on.exit(options(psup.trace_alpha = old), add = TRUE)
  }

  p <- if (dataset == "wallace") {
    plots <- wallace_plot_from_input(ids, show_traces = traces,
               errorbars = errorbars, error_type = et, idType = idType)
    if (plot == "time") plots$plot_time else plots$plot_temp
  } else {
    .preview_dia(ids, plot = plot, traces = traces, bioreps = bioreps,
                 errorbars = errorbars, error_type = et,
                 species = species, base_temps = base_temps, end_temp = end_temp,
                 temp_range = temp_range, time_range = time_range)
  }

  if (isTRUE(print)) print(p)
  invisible(p)
}

# ── Alpha sweep: render the same view at several trace alphas ─────────────────
preview_alpha_sweep <- function(ids, alphas = c(0.05, 0.1, 0.2, 0.3), ...) {
  plots <- lapply(alphas, function(a) {
    p <- preview_psup(ids, trace_alpha = a, print = FALSE, ...) +
      ggplot2::ggtitle(sprintf("trace_alpha = %.2f", a))
    print(p)          # each lands in the RStudio plot history to scrub through
    p
  })
  names(plots) <- sprintf("alpha_%.2f", alphas)
  invisible(plots)
}

# ── DIA preview ───────────────────────────────────────────────────────────────
# The DIA plot is assembled inside app.R's Shiny server, so it can't be called
# directly. This reproduces that pipeline using the app's own shared helpers
# (dia_match_genes / dia_add_categories / scale_colour_dark2 / …) and the same
# getOption("psup.trace_alpha") knob. Keep .dia_build() in sync with
# dia_summarize_plot_data() + dia_build_plot() in app.R if those change.
.preview_dia <- function(ids, plot, traces, bioreps, errorbars, error_type,
                         species, base_temps, end_temp, temp_range, time_range) {
  gene_list   <- splitstring(ids)
  genes       <- unique(unlist(gene_list))
  is_category <- length(gene_list) > 1 || any(vapply(gene_list, length, integer(1)) > 1)

  # Species / base-temperature selection (default: everything available).
  sel <- dia_species_base_temps
  if (!is.null(species))    sel <- dplyr::filter(sel, species %in% .env$species)
  if (!is.null(base_temps)) sel <- dplyr::filter(sel, start_temp %in% .env$base_temps)
  if (nrow(sel) == 0) stop("No species/base-temperature combinations selected.")

  dat <- dia_dt %>%
    dplyr::inner_join(sel, by = c("species", "start_temp")) %>%
    dia_match_genes(genes) %>%
    dia_get_display_label()
  if (nrow(dat) == 0) stop("No matching data for those genes / species / base temps.")

  if (plot == "temperature") {
    x_var   <- "end_temp"
    dat     <- dplyr::filter(dat, timepoint == 8 | (start_temp == end_temp & timepoint == 0))
    x_label <- "Temperature (°C) after 8 min."
    x_min   <- temp_range[1]; x_max <- temp_range[2]
    x_breaks <- dia_temp_axis_breaks(x_min, x_max)
  } else {
    x_var <- "timepoint"
    et_shock <- if (!is.null(end_temp)) as.numeric(end_temp) else {
      av <- dat %>% dplyr::filter(end_temp != start_temp) %>% dplyr::pull(end_temp)
      if (length(av)) max(av) else 37
    }
    dat     <- dplyr::filter(dat, end_temp == et_shock | (end_temp == start_temp & timepoint == 0))
    x_label <- paste0("Minutes at ", et_shock, "°C")
    x_min   <- time_range[1]; x_max <- time_range[2]
    x_breaks <- sort(unique(c(seq(x_min, x_max, by = 4),
                  unique(dat$timepoint[dat$timepoint >= x_min & dat$timepoint <= x_max]))))
  }
  if (is_category) dat <- dia_add_categories(dat, gene_list)
  if (nrow(dat) == 0) stop("No data left after filtering for that plot type / range.")

  .dia_build(dat, x_var, x_label, x_min, x_max, x_breaks,
             is_category, gene_list, traces, bioreps, errorbars, error_type)
}

.dia_build <- function(dat, x_var, x_label, x_min, x_max, x_breaks,
                       is_category, gene_list, show_traces, show_bioreps,
                       show_errorbars, error_type) {
  label_col <- if (is_category && "category" %in% names(dat)) "category" else "display_label"

  mean_dat <- dat %>%
    dplyr::group_by(.data[[label_col]], species, .data[[x_var]]) %>%
    dplyr::summarise(sd = sd(pSup, na.rm = TRUE),
                     se = sd / sqrt(sum(!is.na(pSup))),
                     pSup = mean(pSup, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(display_label = .data[[label_col]])

  trace_dat <- NULL
  if (show_traces && is_category) {
    multi <- names(gene_list)[vapply(gene_list, length, integer(1)) > 1]
    if (length(multi) > 0) {
      trace_dat <- dat %>%
        dplyr::filter(category %in% multi) %>%
        dplyr::group_by(category, species, orf, .data[[x_var]]) %>%
        dplyr::summarise(pSup = mean(pSup, na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(display_label = category)
    }
  }

  p <- ggplot(mean_dat, aes(
    x = .data[[x_var]], y = pSup,
    colour = display_label, shape = species,
    group = interaction(display_label, species)
  ))

  if (!is.null(trace_dat) && nrow(trace_dat) > 0) {
    p <- p + geom_line(
      data = trace_dat,
      aes(x = .data[[x_var]], y = pSup, colour = display_label,
          group = interaction(display_label, species, orf)),
      alpha = getOption("psup.trace_alpha", 0.3), linewidth = 0.8,
      inherit.aes = FALSE)
  }

  p <- p + geom_line(linewidth = 0.8) + geom_point(size = 3)

  if (show_errorbars) {
    err_dat <- mean_dat %>% dplyr::mutate(
      .err = .data[[error_type]],
      pSup_lo = pmax(0, pSup - .err), pSup_hi = pmin(1, pSup + .err))
    p <- p + geom_errorbar(
      data = err_dat, aes(ymin = pSup_lo, ymax = pSup_hi),
      width = (x_max - x_min) * 0.015, linewidth = 0.4)
  }

  if (show_bioreps) {
    p <- p + geom_point(
      data = dat %>% dplyr::mutate(display_label = .data[[label_col]]),
      aes(x = .data[[x_var]], y = pSup, colour = display_label, shape = species),
      size = 2, alpha = 0.4, inherit.aes = FALSE)
  }

  p +
    scale_colour_dark2() +
    scale_shape_manual(values = dia_species_shapes,
                       labels = dia_species_labels_md, name = "Species") +
    labs(colour = NULL) +
    xlab(x_label) +
    scale_y_pSup() +
    scale_x_continuous(breaks = x_breaks, expand = expansion(mult = 0.02)) +
    coord_cartesian(xlim = c(x_min, x_max)) +
    theme(legend.position = "right", legend.text = element_markdown())
}
