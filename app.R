# Merged pSup Data Browser
# Supports two datasets via URL parameter ?dataset=wallace (default) or ?dataset=dia
#
# Wallace et al., Cell 162(6), 2015 — S. cerevisiae heat stress
# Keyport Kik et al., Nature Communications 15:3127, 2024 — three-species DIA

library(shiny)
library(ggplot2)
#library(ggrepel)
library(tidyverse)
library(ggtext)

# ══════════════════════════════════════════════════════════════════════════════
# Data loading
# ══════════════════════════════════════════════════════════════════════════════

wallace_dt <- read_tsv("scer_aggregation_psup_long.txt", comment = "#",
                       show_col_types = FALSE)

dia_dt <- read_tsv("dia/data/dia_psup_data.tsv", show_col_types = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# Shared utilities
# ══════════════════════════════════════════════════════════════════════════════

theme_set(theme_minimal(base_size = 14))

splitstring <- function(s) {
  classes <- str_split(s, ";")[[1]]
  sapply(classes, function(ss) {
    res <- str_split(ss, "=")[[1]]
    if (length(res) == 2) {
      # Keep spaces inside an explicit label ("glycolytic enzymes"), trimming
      # only the edges. Gene names below still have all whitespace removed.
      label <- trimws(res[1])
      genenames <- res[2]
    } else {
      trimname <- gsub("[[:space:]]", "", res[1])
      label <- toupper(trimname)
      genenames <- trimname
    }
    genes <- toupper(gsub("[[:space:]]", "", str_split(genenames, ",")[[1]]))
    g <- list(g = genes)
    names(g)[1] <- label
    g
  }, USE.NAMES = FALSE)
}

scale_y_pSup <- function() {
  list(
    scale_y_continuous(
      name = "Proportion in the\nsupernatant (pSup)",
      expand = c(0.01, 0), limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1)
    ),
    theme(axis.title.y = element_text(angle = 90, vjust = 0.5))
  )
}

# Colorblind-friendly ColorBrewer Dark2 palette. Unlike scale_colour_brewer(),
# this recycles the 8 Dark2 colors when there are more than 8 series so that
# every series is still drawn (brewer drops series past the 8th to NA).
scale_colour_dark2 <- function(...) {
  dark2 <- scales::brewer_pal(palette = "Dark2")(8)
  discrete_scale(aesthetics = "colour",
                 palette = function(n) rep(dark2, length.out = n), ...)
}

# ── Shared plot styling ───────────────────────────────────────────────────────
# The single place these two knobs are defined. Both apps read them at draw
# time, so editing a default here — or setting the option for a session, e.g.
# options(psup.dodge_frac = 0.04) — propagates everywhere, including preview.R.

# Opacity of the faint individual-gene traces drawn behind each group mean.
psup_trace_alpha <- function() getOption("psup.trace_alpha", 0.15)

# Sideways offset applied to error bars so that overlapping bars stay
# distinguishable, as a fraction of the x extent. A dodge rather than a random
# jitter, so figures reproduce between renders.
psup_dodge_frac <- function() getOption("psup.dodge_frac", 0.02)

# Minimum vertical gap between end-of-line labels, in pSup units. Approximate:
# the space a label really needs also depends on the plot's pixel height.
psup_label_gap <- function() getOption("psup.label_gap", 0.045)

# ══════════════════════════════════════════════════════════════════════════════
# Wallace-specific code
# ══════════════════════════════════════════════════════════════════════════════

scale_time <- function(name = expression("Minutes at " * 46 * degree ~ C * ""),
                       text = TRUE, min_x = 0, max_x = NULL) {
  # limits here filter data (unlike coord_cartesian), so min_x must leave room
  # for any x dodge applied to the error bars or those points get dropped.
  if (is.null(max_x)) max_x <- if (text) 9.9 else 8
  scale_x_continuous(name, expand = c(0.01, 0.1), limits = c(min_x, max_x),
                     breaks = c(0, 2, 4, 8))
}

# Room to reserve to the right of the last timepoint for the end-of-line labels,
# in x units. Long category names ("glycolytic enzymes") would otherwise be
# clipped at the panel edge. 9.9 is kept as a floor so short labels are
# unaffected. Approximate: label width also depends on the plot's pixel width.
wallace_label_room <- function(labels, times) {
  widest <- suppressWarnings(max(nchar(as.character(labels)), 0))
  max(9.9, max(times) + 0.03 * diff(range(times)) + 0.15 * widest)
}

# Nudge end-of-line labels apart along y so none overlap, moving each as little
# as possible; labels already clear of their neighbours do not move at all.
# This is a pure function of the input positions, so the layout is identical on
# every re-render — unlike ggrepel, which runs a force simulation at draw time
# (it drifts labels off their own line even when nothing collides, re-lays them
# out whenever the plot is resized, and silently drops labels when crowded).
spread_label_y <- function(y, min_gap = psup_label_gap(), lo = 0, hi = 1) {
  n <- length(y)
  if (n < 2 || min_gap <= 0) return(y)
  ord <- order(y)
  ys <- y[ord]
  # Sweep up, pushing each label clear of the one below it.
  for (i in 2:n) ys[i] <- max(ys[i], ys[i - 1] + min_gap)
  # If the stack overflowed the top, shift it down and sweep back the other way.
  over <- ys[n] - hi
  if (over > 0) {
    ys <- ys - over
    for (i in rev(seq_len(n - 1))) ys[i] <- min(ys[i], ys[i + 1] - min_gap)
    ys <- pmax(ys, lo)
  }
  out <- numeric(n)
  out[ord] <- ys
  out
}

plotgenes_categories <- function(gene_cat_list, data = wallace_dt,
                                 tempexps = c("30C.rep2", "37C.8min", "42C.8min", "46C.8min"),
                                 temps = c(30, 37, 42, 46), tempbreaks = temps,
                                 timeexps = c("30C.rep1", "46C.2min", "46C.4min", "46C.8min"),
                                 times = c(0, 2, 4, 8),
                                 errorbars = FALSE, error_type = "se",
                                 show_traces = FALSE, linewidth = 0.8,
                                 idType = c("gene", "orf")) {
  names(temps) <- tempexps
  names(times) <- timeexps

  if (!is_list(gene_cat_list)) {
    gene_cat_list <- list("<unnamed>" = gene_cat_list)
  }

  ps_dt_temp_all <- bind_rows(lapply(names(gene_cat_list), function(cat) {
    mygenes <- gene_cat_list[[cat]]
    wallace_dt |> filter(.data[[idType]] %in% mygenes & experiment %in% tempexps) |>
      mutate(category = cat)
  }))

  ps_dt_time_all <- bind_rows(lapply(names(gene_cat_list), function(cat) {
    mygenes <- gene_cat_list[[cat]]
    wallace_dt |> filter(.data[[idType]] %in% mygenes & experiment %in% timeexps) |>
      mutate(category = cat)
  }))

  # Per-gene x positions, used to draw individual traces within each category.
  ps_dt_temp_all$temp <- temps[ps_dt_temp_all$experiment]
  ps_dt_time_all$time <- times[ps_dt_time_all$experiment]

  # Colour is a discrete scale, so ggplot assigns palette slots in level order.
  # Ordering categories by their position in the input (rather than the default
  # alphabetical order) keeps existing colours stable when a category is added.
  # The summaries below carry this factor into the gene/orf colour columns.
  ps_dt_temp_all$category <- factor(ps_dt_temp_all$category,
                                    levels = names(gene_cat_list))
  ps_dt_time_all$category <- factor(ps_dt_time_all$category,
                                    levels = names(gene_cat_list))

  # SEM uses the number of genes actually detected (non-NA) in the category.
  ps_dt_temp <- ps_dt_temp_all |>
    group_by(experiment, category) |>
    summarise(sd = sd(psup, na.rm = TRUE), se = sd / sqrt(sum(!is.na(psup))),
              psup = mean(psup, na.rm = TRUE),
              gene = first(category), orf = first(category), .groups = "drop")

  ps_dt_time <- ps_dt_time_all |>
    group_by(experiment, category) |>
    summarise(sd = sd(psup, na.rm = TRUE), se = sd / sqrt(sum(!is.na(psup))),
              psup = mean(psup, na.rm = TRUE),
              gene = first(category), orf = first(category), .groups = "drop")

  ps_dt_temp$temp <- temps[ps_dt_temp$experiment]
  ps_dt_time$time <- times[ps_dt_time$experiment]

  plotgenes_wallace(ps_dt_time, ps_dt_temp, ps_dt_time_all, ps_dt_temp_all,
                    tempexps, temps, tempbreaks, timeexps, times,
                    errorbars, error_type, show_traces, linewidth, idType)
}

plotgenes_wallace <- function(ps_dt_time, ps_dt_temp,
                              ps_dt_time_all, ps_dt_temp_all,
                              tempexps = c("30C.rep2", "37C.8min", "42C.8min", "46C.8min"),
                              temps = c(30, 37, 42, 46), tempbreaks = temps,
                              timeexps = c("30C.rep1", "46C.2min", "46C.4min", "46C.8min"),
                              times = c(0, 2, 4, 8),
                              errorbars = FALSE, error_type = "se",
                              show_traces = FALSE, linewidth = 0.8,
                              idType = c("gene", "orf")) {
  # Selected error statistic (SEM or SD), used for the error-bar half-widths.
  ps_dt_temp$err <- ps_dt_temp[[error_type]]
  ps_dt_time$err <- ps_dt_time[[error_type]]

  # Offset each category's error bars slightly along x so that overlapping bars
  # stay readable. A dodge (not a random jitter) keeps figures reproducible.
  # Width is a fraction of the x span, so it scales with the axis.
  dodge_frac <- psup_dodge_frac()
  dodge_temp <- position_dodge(width = diff(range(temps)) * dodge_frac)
  dodge_time <- position_dodge(width = diff(range(times)) * dodge_frac)

  # End-of-line labels sit at each series' final value, nudged apart in y only
  # where they would otherwise collide.
  lab_temp <- ps_dt_temp |> filter(temp == max(temps))
  lab_time <- ps_dt_time |> filter(time == max(times))
  lab_temp$psup_lab <- spread_label_y(lab_temp$psup)
  lab_time$psup_lab <- spread_label_y(lab_time$psup)

  # Individual traces are only meaningful for categories with >1 gene; a
  # single-gene category's trace would just overprint its own mean line.
  traces_temp <- ps_dt_temp_all |> group_by(category) |>
    filter(n_distinct(orf) > 1) |> ungroup()
  traces_time <- ps_dt_time_all |> group_by(category) |>
    filter(n_distinct(orf) > 1) |> ungroup()

  # \u2500\u2500 Temperature plot \u2500\u2500
  plot_temp <- ggplot(data = ps_dt_temp,
    aes(x = .data[["temp"]], y = .data[["psup"]],
        colour = .data[[idType]], label = .data[[idType]]))

  if (show_traces && nrow(traces_temp) > 0) {
    plot_temp <- plot_temp + geom_line(
      data = traces_temp,
      aes(x = temp, y = psup, colour = category,
          group = interaction(category, orf)),
      alpha = psup_trace_alpha(), linewidth = linewidth,
      inherit.aes = FALSE)
  }

  plot_temp <- plot_temp + geom_line(linewidth = linewidth)

  if (errorbars) {
    plot_temp <- plot_temp +
      geom_pointrange(aes(ymin = psup - err, ymax = psup + err),
                      position = dodge_temp)
  }

  plot_temp <- plot_temp +
    # Labels in a fixed column at each series' final value, nudged vertically
    # only where they would collide (see spread_label_y).
    geom_text(size = 4, hjust = 0, data = lab_temp,
      aes(x = max(temps) + 0.03 * diff(range(temps)), y = psup_lab)) +
    # Start just below the lowest temperature: the x scale has no expansion, so
    # data sitting exactly on the panel edge would be clipped in half.
    coord_cartesian(xlim = c(min(temps) - 1, 52)) +
    scale_x_continuous("Temperature (\u00b0C) of 8 min. treatment",
                       breaks = tempbreaks, labels = tempbreaks, expand = c(0, 0)) +
    scale_y_pSup() +
    scale_colour_dark2() +
    theme(legend.position = "none")

  # \u2500\u2500 Time plot \u2500\u2500
  plot_time <- ggplot(data = ps_dt_time,
    aes(x = .data[["time"]], y = .data[["psup"]],
        colour = .data[[idType]], label = .data[[idType]]))

  if (show_traces && nrow(traces_time) > 0) {
    plot_time <- plot_time + geom_line(
      data = traces_time,
      aes(x = time, y = psup, colour = category,
          group = interaction(category, orf)),
      alpha = psup_trace_alpha(), linewidth = linewidth,
      inherit.aes = FALSE)
  }

  plot_time <- plot_time + geom_line(linewidth = linewidth)

  if (errorbars) {
    plot_time <- plot_time +
      geom_pointrange(aes(ymin = psup - err, ymax = psup + err),
                      position = dodge_time)
  }

  time_max_x <- wallace_label_room(ps_dt_time[[idType]], times)

  plot_time <- plot_time +
    geom_text(size = 4, hjust = 0, data = lab_time,
      aes(x = max(times) + 0.03 * diff(range(times)), y = psup_lab)) +
    scale_y_pSup() +
    # Only widen the lower limit when dodged error bars need the room, so the
    # default view still starts flush at time 0.
    scale_time(min_x = if (errorbars)
                 min(times) - diff(range(times)) * dodge_frac else min(times),
               max_x = time_max_x) +
    scale_colour_dark2() +
    theme(legend.position = "none")

  list(plot_time = plot_time, plot_temp = plot_temp)
}

wallace_plot_from_input <- function(s, show_traces = FALSE, errorbars = FALSE,
                                    error_type = "se", idType) {
  ids <- splitstring(s)
  plotgenes_categories(ids, errorbars = errorbars, error_type = error_type,
                       show_traces = show_traces, idType = idType)
}

# ══════════════════════════════════════════════════════════════════════════════
# DIA-specific code
# ══════════════════════════════════════════════════════════════════════════════

dia_all_species <- c("S. cerevisiae", "S. kudriavzevii", "K. marxianus")

dia_species_shapes <- c(
  "S. cerevisiae"   = 16,
  "S. kudriavzevii" = 17,
  "K. marxianus"    = 15
)

dia_species_labels_md <- c(
  "S. cerevisiae"   = "*S. cerevisiae*",
  "S. kudriavzevii" = "*S. kudriavzevii*",
  "K. marxianus"    = "*K. marxianus*"
)

dia_data_temps <- c(23, 30, 37, 42, 46, 50)

dia_species_base_temps <- dia_dt %>%
  distinct(species, start_temp) %>%
  arrange(species, start_temp)

dia_condition_map <- dia_dt %>%
  distinct(species, start_temp, end_temp, timepoint) %>%
  arrange(species, start_temp, end_temp, timepoint)

dia_match_genes <- function(data, query_genes) {
  data %>% filter(
    toupper(gene) %in% query_genes |
      toupper(scer_ortholog) %in% query_genes |
      toupper(orf) %in% query_genes
  )
}

dia_get_display_label <- function(data) {
  data %>% mutate(display_label = case_when(
    !is.na(gene) & gene != "" ~ gene,
    !is.na(scer_ortholog) & scer_ortholog != "" ~ scer_ortholog,
    TRUE ~ orf
  ))
}

dia_temp_axis_breaks <- function(t_min, t_max) {
  regular <- seq(floor(t_min / 5) * 5, ceiling(t_max / 5) * 5, by = 5)
  sort(unique(c(regular, dia_data_temps[dia_data_temps >= t_min & dia_data_temps <= t_max])))
}

dia_add_categories <- function(dat, cats) {
  gene_to_cat <- bind_rows(lapply(names(cats), function(cat_name) {
    tibble(query_gene = cats[[cat_name]], category = cat_name)
  }))
  dat %>%
    mutate(query_match = case_when(
      toupper(gene) %in% gene_to_cat$query_gene ~ toupper(gene),
      toupper(scer_ortholog) %in% gene_to_cat$query_gene ~ toupper(scer_ortholog),
      toupper(orf) %in% gene_to_cat$query_gene ~ toupper(orf),
      TRUE ~ NA_character_
    )) %>%
    left_join(gene_to_cat, by = c("query_match" = "query_gene")) %>%
    filter(!is.na(category)) %>%
    # Order categories by their position in the input so that adding one leaves
    # the existing categories' colours unchanged (see plotgenes_categories).
    mutate(category = factor(category, levels = names(cats))) %>%
    select(-query_match)
}

# ── Shared DIA plot pipeline ──────────────────────────────────────────────────
# Plain functions, deliberately free of Shiny reactives, so the app's server and
# command-line callers (preview.R) drive exactly the same code. All DIA plotting
# logic belongs here rather than inside server().

dia_is_category_list <- function(gene_list) {
  length(gene_list) > 1 || any(vapply(gene_list, length, integer(1)) > 1)
}

# Rows for the selected species / base temperatures and the query genes.
dia_filter_data <- function(sel, genes) {
  dia_dt %>%
    inner_join(sel, by = c("species", "start_temp")) %>%
    dia_match_genes(genes) %>%
    dia_get_display_label()
}

# Narrow to the rows a temperature or time plot draws, then attach categories.
dia_plot_subset <- function(dat, gene_list, plot_type, shock_temp = NULL) {
  dat <- if (plot_type == "temperature") {
    dat %>% filter(timepoint == 8 | (start_temp == end_temp & timepoint == 0))
  } else {
    st <- as.numeric(shock_temp)
    dat %>% filter(end_temp == st | (end_temp == start_temp & timepoint == 0))
  }
  if (nrow(dat) == 0) return(dat)
  if (dia_is_category_list(gene_list)) dat <- dia_add_categories(dat, gene_list)
  dat
}

# Axis variable, label, limits and breaks for each plot type.
dia_plot_axis <- function(dat, plot_type, shock_temp = NULL,
                          x_min = NULL, x_max = NULL) {
  if (plot_type == "temperature") {
    if (is.null(x_min)) x_min <- 23
    if (is.null(x_max)) x_max <- 50
    list(x_var = "end_temp", x_min = x_min, x_max = x_max,
         x_label = "Temperature (°C) after 8 min.",
         x_breaks = dia_temp_axis_breaks(x_min, x_max))
  } else {
    if (is.null(x_min)) x_min <- 0
    if (is.null(x_max)) x_max <- 20
    list(x_var = "timepoint", x_min = x_min, x_max = x_max,
         x_label = paste0("Minutes at ", as.numeric(shock_temp), "°C"),
         x_breaks = sort(unique(c(
           seq(x_min, x_max, by = 4),
           unique(dat$timepoint[dat$timepoint >= x_min & dat$timepoint <= x_max])))))
  }
}

# Mean per series, plus optional raw biorep points and per-gene traces.
dia_summarize_plot_data <- function(dat, x_var, gene_list,
                                    show_bioreps = FALSE, show_traces = FALSE) {
  is_category <- dia_is_category_list(gene_list)
  label_col <- if (is_category && "category" %in% names(dat)) "category" else "display_label"

  mean_dat <- dat %>%
    group_by(.data[[label_col]], species, .data[[x_var]]) %>%
    summarise(sd = sd(pSup, na.rm = TRUE),
              se = sd / sqrt(sum(!is.na(pSup))),
              pSup = mean(pSup, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(display_label = .data[[label_col]])

  biorep_dat <- if (isTRUE(show_bioreps)) {
    dat %>% mutate(display_label = .data[[label_col]])
  } else NULL

  # Traces only for categories holding more than one gene; a single-gene
  # category's trace would just overprint its own mean line.
  trace_dat <- if (isTRUE(show_traces) && is_category) {
    multi <- names(gene_list)[vapply(gene_list, length, integer(1)) > 1]
    if (length(multi) > 0) {
      dat %>%
        filter(category %in% multi) %>%
        group_by(category, species, orf, .data[[x_var]]) %>%
        summarise(pSup = mean(pSup, na.rm = TRUE), .groups = "drop") %>%
        mutate(display_label = category)
    } else NULL
  } else NULL

  list(mean = mean_dat, bioreps = biorep_dat, traces = trace_dat)
}

# Assemble the ggplot from an already-subset data frame.
dia_assemble_plot <- function(dat, axis, gene_list, show_traces = FALSE,
                              show_bioreps = FALSE, show_errorbars = FALSE,
                              error_type = "se") {
  x_var <- axis$x_var
  plot_data <- dia_summarize_plot_data(dat, x_var, gene_list,
                                       show_bioreps = show_bioreps,
                                       show_traces = show_traces)
  mean_dat   <- plot_data$mean
  biorep_dat <- plot_data$bioreps
  trace_dat  <- plot_data$traces

  p <- ggplot(mean_dat, aes(
    x = .data[[x_var]], y = pSup,
    colour = display_label, shape = species,
    group = interaction(display_label, species)
  ))

  # Faint individual gene traces, drawn under the group mean line.
  if (!is.null(trace_dat) && nrow(trace_dat) > 0) {
    p <- p + geom_line(
      data = trace_dat,
      aes(x = .data[[x_var]], y = pSup, colour = display_label,
          group = interaction(display_label, species, orf)),
      alpha = psup_trace_alpha(), linewidth = 0.8,
      inherit.aes = FALSE)
  }

  # With error bars shown, nudge each series' mean point and bar sideways
  # together so overlapping bars stay readable; the line stays at its true x.
  dodge <- position_dodge(width = (axis$x_max - axis$x_min) * psup_dodge_frac())
  point_pos <- if (isTRUE(show_errorbars)) dodge else "identity"

  p <- p + geom_line(linewidth = 0.8) +
    geom_point(size = 3, position = point_pos)

  if (isTRUE(show_errorbars)) {
    err_dat <- mean_dat %>% mutate(
      .err = .data[[error_type]],
      pSup_lo = pmax(0, pSup - .err), pSup_hi = pmin(1, pSup + .err))
    p <- p + geom_errorbar(
      data = err_dat,
      aes(ymin = pSup_lo, ymax = pSup_hi),
      width = (axis$x_max - axis$x_min) * 0.015, linewidth = 0.4,
      position = dodge)
  }

  if (!is.null(biorep_dat)) {
    p <- p + geom_point(
      data = biorep_dat,
      aes(x = .data[[x_var]], y = pSup,
          colour = display_label, shape = species),
      size = 2, alpha = 0.4, inherit.aes = FALSE)
  }

  p +
    scale_colour_dark2() +
    scale_shape_manual(values = dia_species_shapes,
                       labels = dia_species_labels_md,
                       name = "Species") +
    labs(colour = NULL) +
    xlab(axis$x_label) +
    scale_y_pSup() +
    scale_x_continuous(breaks = axis$x_breaks, expand = expansion(mult = 0.02)) +
    coord_cartesian(xlim = c(axis$x_min, axis$x_max)) +
    theme(legend.position = "right", legend.text = element_markdown())
}

# Gene list + species/base-temp selection -> ggplot. The single entry point used
# by both the Shiny server and preview.R. Returns NULL when nothing matches.
dia_plot <- function(gene_list, sel, plot_type = c("temperature", "time"),
                     shock_temp = NULL, x_min = NULL, x_max = NULL,
                     show_traces = FALSE, show_bioreps = FALSE,
                     show_errorbars = FALSE, error_type = "se") {
  plot_type <- match.arg(plot_type)
  if (nrow(sel) == 0) return(NULL)
  dat <- dia_filter_data(sel, unique(unlist(gene_list)))
  if (nrow(dat) == 0) return(NULL)
  dat <- dia_plot_subset(dat, gene_list, plot_type, shock_temp)
  if (nrow(dat) == 0) return(NULL)
  axis <- dia_plot_axis(dat, plot_type, shock_temp, x_min, x_max)
  dia_assemble_plot(dat, axis, gene_list, show_traces = show_traces,
                    show_bioreps = show_bioreps,
                    show_errorbars = show_errorbars, error_type = error_type)
}

# ══════════════════════════════════════════════════════════════════════════════
# UI builder functions
# ══════════════════════════════════════════════════════════════════════════════

wallace_ui <- function() {
  sidebarLayout(
    sidebarPanel(
      textAreaInput("ids", "Gene identifiers separated by semicolons:",
                    value = "PGK1;PAB1;PMA1", rows = 4),
      helpText("Examples (click to update):"),
      actionLink("flat_examples", "Individual proteins"),
      br(),
      actionLink("category_examples", "Protein categories"),
      helpText(" "),
      selectInput("idType", "Identify by:", c("gene", "orf"), selected = "gene"),
      checkboxInput("show_traces", "Show individual traces", FALSE),
      checkboxInput("show_errorbars", "Show error bars", FALSE),
      conditionalPanel(
        condition = "input.show_errorbars",
        radioButtons("error_type", NULL,
                     choices = c("SEM" = "se", "SD" = "sd"),
                     selected = "se", inline = TRUE)
      ),
      selectInput("plotType", "Plot against:", c("time", "temperature"), selected = "time"),
      hr(),
      h5("Plot settings"),
      fluidRow(
        column(6, numericInput("plot_width", "Width (px):", value = 700, min = 300, max = 1600, step = 50)),
        column(6, numericInput("plot_height", "Height (px):", value = 500, min = 200, max = 1200, step = 50))
      )
    ),
    mainPanel(
      uiOutput("plot_ui"),
      fluidRow(
        column(3, selectInput("dl_format", NULL,
          choices = c("PNG" = "png", "SVG" = "svg", "PDF" = "pdf"),
          selected = "png", width = "100%")),
        column(3, downloadButton("dl_plot", "Download figure"))
      )
    )
  )
}

dia_ui <- function() {
  sidebarLayout(
    sidebarPanel(
      textAreaInput("ids", "Gene identifiers separated by semicolons:",
                    value = "glycolytic=PGK1,FBA1,TDH3;superaggregator=DED1,NUG1,OLA1;ribosomal=RPL4A,RPL19A,RPS9B", rows = 4),
      helpText("Examples (click to update):"),
      actionLink("flat_examples", "Individual proteins"),
      actionLink("category_examples", "Protein categories"),
      helpText(" "),
      selectInput("plotType", "Plot against:",
                  c("temperature", "time"), selected = "temperature"),

      h5("Species & base temperatures:"),
      uiOutput("species_base_temp_ui"),

      conditionalPanel(
        condition = "input.plotType == 'time'",
        uiOutput("end_temp_ui")
      ),

      checkboxInput("show_traces", "Show individual traces", FALSE),
      checkboxInput("show_bioreps", "Show individual bioreps", FALSE),
      checkboxInput("show_errorbars", "Show error bars", FALSE),
      conditionalPanel(
        condition = "input.show_errorbars",
        radioButtons("error_type", NULL,
                     choices = c("SEM" = "se", "SD" = "sd"),
                     selected = "se", inline = TRUE)
      ),

      hr(),
      h5("Plot settings"),
      fluidRow(
        column(6, numericInput("plot_width", "Width (px):", value = 700, min = 300, max = 1600, step = 50)),
        column(6, numericInput("plot_height", "Height (px):", value = 500, min = 200, max = 1200, step = 50))
      ),
      conditionalPanel(
        condition = "input.plotType == 'temperature'",
        fluidRow(
          column(6, numericInput("temp_min", "Temp min (\u00b0C):", value = 23, min = 0, max = 60, step = 1)),
          column(6, numericInput("temp_max", "Temp max (\u00b0C):", value = 50, min = 0, max = 60, step = 1))
        )
      ),
      conditionalPanel(
        condition = "input.plotType == 'time'",
        fluidRow(
          column(6, numericInput("time_min", "Time min (min):", value = 0, min = 0, max = 60, step = 1)),
          column(6, numericInput("time_max", "Time max (min):", value = 20, min = 0, max = 60, step = 1))
        )
      ),

      hr(),
      helpText(
        "Enter S. cerevisiae gene names, systematic ORFs (YxxNNNx),",
        "or native K. marxianus / S. kudriavzevii identifiers.",
        "Use category=gene1,gene2 syntax for grouped averages."
      ),
      width = 3
    ),
    mainPanel(
      uiOutput("plot_ui"),
      fluidRow(
        column(3, selectInput("dl_format", NULL,
          choices = c("PNG" = "png", "SVG" = "svg", "PDF" = "pdf"),
          selected = "png", width = "100%")),
        column(3, downloadButton("dl_plot", "Download figure"))
      ),
      htmlOutput("data_info"),
      width = 9
    )
  )
}

wallace_citation <- function() {
  p("Data from: Wallace EWJ, Kear-Scott JL, Pilipenko EV, Schwartz MH, Laskowski PR, Rojek AE,",
    "Katanski CD, Riback JA, Dion MF, Franks AM, Airoldi EM, Pan T, Budnik BA, Drummond DA.",
    "\u201cReversible, specific, active aggregates of endogenous proteins assemble upon heat stress.\u201d",
    tags$em("Cell"),
    "162(6), 2015. ",
    a("doi:10.1016/j.cell.2015.08.041", href = "https://doi.org/10.1016/j.cell.2015.08.041"),
    " | ",
    a("Paper page", href = "https://drummondlab.org/papers/paper/endogenous-aggregates")
  )
}

dia_citation <- function() {
  p("Data from: Keyport Kik S, Christopher D, Glauninger H, Wong Hickernell C, Bard JAM, Lin KM, Squires AH, Ford M, Sosnick TR, Drummond DA.",
    br(),
    "\u201cAn adaptive biomolecular condensation response is conserved across environmentally divergent species.\u201d",
    br(),
    tags$em("Nature Communications"),
    " 15:3127 (2024). ",
    a("doi:10.1038/s41467-024-47355-9", href = "https://doi.org/10.1038/s41467-024-47355-9"),
    " | ",
    a("Paper page", href = "https://drummondlab.org/papers/paper/conserved-condensation")
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════

ui <- function(request) {
  fluidPage(
    uiOutput("app_title"),
    uiOutput("app_body"),
    mainPanel(uiOutput("citation_ui"))
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# Server
# ══════════════════════════════════════════════════════════════════════════════

server <- function(input, output, session) {

  # ── Dataset routing ──────────────────────────────────────────────────────
  dataset <- reactive({
    query <- parseQueryString(session$clientData$url_search)
    d <- query$dataset
    if (is.null(d) || !(d %in% c("wallace", "dia"))) "wallace" else d
  })

  # ── Dynamic UI ───────────────────────────────────────────────────────────
  output$app_title <- renderUI({
    if (dataset() == "dia") {
      titlePanel("Proportion in supernatant (pSup) across yeast species (DIA)")
    } else {
      titlePanel("Proportion in supernatant (pSup) for yeast genes")
    }
  })

  output$app_body <- renderUI({
    if (dataset() == "dia") dia_ui() else wallace_ui()
  })

  output$citation_ui <- renderUI({
    if (dataset() == "dia") dia_citation() else wallace_citation()
  })

  # ── Shared: dynamic plot sizing ──────────────────────────────────────────
  output$plot_ui <- renderUI({
    req(input$plot_width, input$plot_height)
    w <- paste0(input$plot_width, "px")
    h <- paste0(input$plot_height, "px")
    plotOutput("plot", width = w, height = h)
  })

  # ── Shared: current plot reactive ────────────────────────────────────────
  current_plot <- reactive({
    if (dataset() == "dia") {
      dia_current_plot()
    } else {
      plots <- wallace_plot_from_input(input$ids,
        show_traces = isTRUE(input$show_traces),
        errorbars = isTRUE(input$show_errorbars),
        error_type = if (is.null(input$error_type)) "se" else input$error_type,
        idType = input$idType)
      if (input$plotType == "time") plots$plot_time
      else plots$plot_temp
    }
  })

  output$plot <- renderPlot({
    current_plot()
  })

  # ── Shared: figure download ─────────────────────────────────────────────
  output$dl_plot <- downloadHandler(
    filename = function() {
      prefix <- if (dataset() == "dia") "psup_dia_" else "psup_"
      paste0(prefix, Sys.Date(), ".", input$dl_format)
    },
    content = function(file) {
      w_in <- input$plot_width / 150
      h_in <- input$plot_height / 150
      p <- current_plot()
      ggsave(file, plot = p, device = input$dl_format,
             width = w_in, height = h_in,
             dpi = if (input$dl_format == "png") 300 else 72)
    }
  )

  # ── Shared: example links ───────────────────────────────────────────────
  observeEvent(input$flat_examples, {
    updateTextAreaInput(session, "ids", value = "PGK1;PMA1;PAB1")
  })

  observeEvent(input$category_examples, {
    if (dataset() == "dia") {
      updateTextAreaInput(session, "ids",
        value = "glycolytic=PGK1,FBA1,TDH3;superaggregator=DED1,NUG1,OLA1;ribosomal=RPL4A,RPL19A,RPS9B")
    } else {
      updateTextAreaInput(session, "ids",
        value = "glycolytic=PGK1,FBA1,TDH3; aggregator=pab1,xrn1,dcp2; superaggregator=ded1,nug1,ola1; ribosomal=RPL32,RPL4A,RPL19A,RPS9B,RPL21A; membrane=pma1,por1,srp61")
    }
  })

  # ══════════════════════════════════════════════════════════════════════════
  # DIA server logic
  # ══════════════════════════════════════════════════════════════════════════

  # ── Per-species base temperature checkboxes ─────────────────────────────
  output$species_base_temp_ui <- renderUI({
    req(dataset() == "dia")
    ui_elements <- lapply(dia_all_species, function(sp) {
      avail_temps <- dia_species_base_temps %>%
        filter(species == sp) %>% pull(start_temp) %>% sort()
      if (length(avail_temps) == 0) return(NULL)
      choices <- setNames(paste0(sp, "|", avail_temps),
                          paste0(avail_temps, "\u00b0C"))
      checkboxGroupInput(
        inputId = paste0("base_temp_", gsub("[^a-zA-Z]", "", sp)),
        label = tags$i(sp),
        choices = choices, selected = choices, inline = TRUE
      )
    })
    do.call(tagList, ui_elements)
  })

  dia_selected_species_temps <- reactive({
    req(dataset() == "dia")
    combos <- c()
    for (sp in dia_all_species) {
      input_id <- paste0("base_temp_", gsub("[^a-zA-Z]", "", sp))
      vals <- input[[input_id]]
      if (!is.null(vals)) combos <- c(combos, vals)
    }
    if (length(combos) == 0)
      return(tibble(species = character(), start_temp = numeric()))
    tibble(combo = combos) %>%
      separate(combo, into = c("species", "start_temp"), sep = "\\|") %>%
      mutate(start_temp = as.numeric(start_temp))
  })

  # ── End temp selector (time plots) ──────────────────────────────────────
  output$end_temp_ui <- renderUI({
    req(dataset() == "dia")
    sel <- dia_selected_species_temps()
    req(nrow(sel) > 0)
    avail <- dia_condition_map %>%
      inner_join(sel, by = c("species", "start_temp")) %>%
      filter(end_temp != start_temp) %>%
      distinct(end_temp) %>% pull(end_temp) %>% sort()
    if (length(avail) == 0) avail <- c(37)
    selectInput("end_temp", "Heat shock temperature (\u00b0C):",
                choices = avail, selected = max(avail))
  })

  # ── Gene parsing ────────────────────────────────────────────────────────
  dia_gene_list <- reactive({
    req(dataset() == "dia", input$ids)
    splitstring(input$ids)
  })

  dia_all_genes <- reactive({
    unique(unlist(dia_gene_list()))
  })

  # ── DIA plot ────────────────────────────────────────────────────────────
  # All plotting logic lives in the shared dia_plot() pipeline near the top of
  # this file, so the app and preview.R render through identical code.
  dia_current_plot <- function() {
    sel <- dia_selected_species_temps()
    req(nrow(sel) > 0)
    ptype <- input$plotType
    if (ptype == "time") req(input$end_temp)
    is_temp <- ptype == "temperature"
    p <- dia_plot(
      gene_list      = dia_gene_list(),
      sel            = sel,
      plot_type      = ptype,
      shock_temp     = input$end_temp,
      x_min          = if (is_temp) input$temp_min else input$time_min,
      x_max          = if (is_temp) input$temp_max else input$time_max,
      show_traces    = isTRUE(input$show_traces),
      show_bioreps   = isTRUE(input$show_bioreps),
      show_errorbars = isTRUE(input$show_errorbars),
      error_type     = if (is.null(input$error_type)) "se" else input$error_type
    )
    req(!is.null(p))
    p
  }

  # ── DIA data info panel ─────────────────────────────────────────────────
  output$data_info <- renderUI({
    req(dataset() == "dia")
    genes <- dia_all_genes()
    found <- dia_dt %>% dia_match_genes(genes) %>%
      distinct(gene, scer_ortholog, orf)
    found_upper <- c(toupper(found$gene), toupper(found$scer_ortholog),
                     toupper(found$orf))
    missing <- setdiff(genes, found_upper)
    if (length(missing) > 0) {
      HTML(paste0("<b>Not found:</b> ", paste(missing, collapse = ", ")))
    }
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# Run
# ══════════════════════════════════════════════════════════════════════════════

shinyApp(ui, server)
