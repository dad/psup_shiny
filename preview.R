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
# This file holds no plotting logic of its own. It only translates convenient
# arguments into calls to app.R's own functions — wallace_plot_from_input() and
# the shared dia_plot() pipeline — which are the same entry points server()
# uses. So a preview is the app's output, not a lookalike; verified by building
# both and comparing the rendered layer data. Add plotting changes in app.R.
#
# `trace_alpha` sets options(psup.trace_alpha=) for the duration of the render,
# which is the exact knob app.R reads. Set it globally
# (options(psup.trace_alpha = 0.15)) to preview it live in the running app too.
# Every call returns the ggplot invisibly, so you can
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
# No plotting logic here: this only turns preview_psup()'s convenience arguments
# into the arguments dia_plot() expects. The rendering itself is app.R's shared
# dia_plot() pipeline, the exact code the Shiny server calls.
.preview_dia <- function(ids, plot, traces, bioreps, errorbars, error_type,
                         species, base_temps, end_temp, temp_range, time_range) {
  gene_list <- splitstring(ids)

  # Species / base-temperature selection (default: everything available).
  sel <- dia_species_base_temps
  if (!is.null(species))    sel <- dplyr::filter(sel, species %in% .env$species)
  if (!is.null(base_temps)) sel <- dplyr::filter(sel, start_temp %in% .env$base_temps)
  if (nrow(sel) == 0) stop("No species/base-temperature combinations selected.")

  # Time plots need a heat-shock temperature; default to the hottest available
  # for the selected genes, mirroring the app's end-temperature selector.
  if (plot == "time" && is.null(end_temp)) {
    avail <- dia_filter_data(sel, unique(unlist(gene_list)))
    avail <- avail$end_temp[avail$end_temp != avail$start_temp]
    end_temp <- if (length(avail)) max(avail) else 37
  }

  rng <- if (plot == "temperature") temp_range else time_range
  p <- dia_plot(gene_list = gene_list, sel = sel, plot_type = plot,
                shock_temp = end_temp, x_min = rng[1], x_max = rng[2],
                show_traces = traces, show_bioreps = bioreps,
                show_errorbars = errorbars, error_type = error_type)
  if (is.null(p))
    stop("No data for those genes / species / base temps at that plot type.")
  p
}
