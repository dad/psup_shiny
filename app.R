# Merged pSup Data Browser
# Supports two datasets via URL parameter ?dataset=wallace (default) or ?dataset=dia
#
# Wallace et al., Cell 162(6), 2015 — S. cerevisiae heat stress
# Keyport Kik et al., Nature Communications 15:3127, 2024 — three-species DIA

library(shiny)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(ggtext)
library(shinyURL)

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
      label <- gsub("[[:space:]]", "", res[1])
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

# ══════════════════════════════════════════════════════════════════════════════
# Wallace-specific code
# ══════════════════════════════════════════════════════════════════════════════

scale_time <- function(name = expression("Minutes at " * 46 * degree ~ C * ""),
                       text = TRUE) {
  mylim <- if (text) c(0, 9.9) else c(0, 8)
  scale_x_continuous(name, expand = c(0.01, 0.1), limits = mylim,
                     breaks = c(0, 2, 4, 8))
}

plotgenes_categories <- function(gene_cat_list, data = wallace_dt,
                                 tempexps = c("30C.rep2", "37C.8min", "42C.8min", "46C.8min"),
                                 temps = c(30, 37, 42, 46), tempbreaks = temps,
                                 timeexps = c("30C.rep1", "46C.2min", "46C.4min", "46C.8min"),
                                 times = c(0, 2, 4, 8),
                                 errorbars = FALSE, linewidth = 0.8,
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

  ps_dt_temp <- ps_dt_temp_all |>
    group_by(experiment, category) |>
    summarise(sd = sd(psup, na.rm = TRUE), se = sd / sqrt(n()),
              psup = mean(psup), psup.lo = psup - sd, psup.hi = psup + sd,
              gene = first(category), orf = first(category), .groups = "drop")

  ps_dt_time <- ps_dt_time_all |>
    group_by(experiment, category) |>
    summarise(sd = sd(psup, na.rm = TRUE), se = sd / sqrt(n()),
              psup = mean(psup), psup.lo = psup - sd, psup.hi = psup + sd,
              gene = first(category), orf = first(category), .groups = "drop")

  ps_dt_temp$temp <- temps[ps_dt_temp$experiment]
  ps_dt_time$time <- times[ps_dt_time$experiment]

  plotgenes_wallace(ps_dt_time, ps_dt_temp, tempexps, temps, tempbreaks,
                    timeexps, times, errorbars, linewidth, idType)
}

plotgenes_wallace <- function(ps_dt_time, ps_dt_temp,
                              tempexps = c("30C.rep2", "37C.8min", "42C.8min", "46C.8min"),
                              temps = c(30, 37, 42, 46), tempbreaks = temps,
                              timeexps = c("30C.rep1", "46C.2min", "46C.4min", "46C.8min"),
                              times = c(0, 2, 4, 8),
                              errorbars = FALSE, linewidth = 0.8,
                              idType = c("gene", "orf")) {
  plot_temp <- ggplot(data = ps_dt_temp,
    aes(x = .data[["temp"]], y = .data[["psup"]],
        ymin = .data[["psup.lo"]], ymax = .data[["psup.hi"]],
        colour = .data[[idType]], label = .data[[idType]])) +
    geom_line(linewidth = linewidth) +
    geom_text_repel(size = 4,
      data = ps_dt_temp |> filter(temp == max(temps)),
      aes(x = max(temps) + 0.5, y = psup), xlim = c(46, 52)) +
    coord_cartesian(xlim = c(30, 52)) +
    scale_x_continuous("Temperature (\u00b0C) of 8 min. treatment",
                       breaks = tempbreaks, labels = tempbreaks, expand = c(0, 0)) +
    scale_y_pSup() +
    scale_colour_brewer(palette = "Dark2") +
    theme(legend.position = "none")

  plot_time <- ggplot(data = ps_dt_time,
    aes(x = .data[["time"]], y = .data[["psup"]],
        ymin = .data[["psup.lo"]], ymax = .data[["psup.hi"]],
        colour = .data[[idType]], label = .data[[idType]])) +
    geom_line(linewidth = linewidth) +
    geom_text_repel(size = 4,
      data = subset(ps_dt_time, time == max(times)),
      aes(x = max(times) + 0.1, y = psup), xlim = c(8, 12)) +
    scale_y_pSup() + scale_time() +
    scale_colour_brewer(palette = "Dark2") +
    theme(legend.position = "none")

  if (errorbars) {
    plot_temp <- plot_temp + geom_pointrange()
    plot_time <- plot_time + geom_pointrange()
  }

  list(plot_time = plot_time, plot_temp = plot_temp)
}

wallace_plot_from_input <- function(s, errorbars, idType) {
  ids <- splitstring(s)
  plotgenes_categories(ids, errorbars = errorbars, idType = idType)
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
    select(-query_match)
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
      actionLink("category_examples", "Protein categories"),
      helpText(" "),
      selectInput("idType", "Identify by:", c("gene", "orf"), selected = "gene"),
      checkboxInput("interval", "Show 95% intervals", FALSE),
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
                    value = "PGK1;PAB1;DED1", rows = 4),
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

      checkboxInput("show_bioreps", "Show individual bioreps", FALSE),
      conditionalPanel(
        condition = "!input.show_bioreps",
        checkboxInput("show_errorbars", "Show \u00b1SD error bars", FALSE)
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
    mainPanel(uiOutput("citation_ui")),
    shinyURL.ui()
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
      if (input$plotType == "temperature") dia_build_temp_plot()
      else dia_build_time_plot()
    } else {
      plots <- wallace_plot_from_input(input$ids,
        errorbars = input$interval, idType = input$idType)
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

  dia_is_category_mode <- reactive({
    cats <- dia_gene_list()
    length(cats) > 1 || any(sapply(cats, length) > 1)
  })

  # ── Filtered data ───────────────────────────────────────────────────────
  dia_filtered_data <- reactive({
    sel <- dia_selected_species_temps()
    req(nrow(sel) > 0)
    genes <- dia_all_genes()
    dia_dt %>%
      inner_join(sel, by = c("species", "start_temp")) %>%
      dia_match_genes(genes) %>%
      dia_get_display_label()
  })

  # ── Temperature plot data ───────────────────────────────────────────────
  dia_temp_plot_data <- reactive({
    dat <- dia_filtered_data()
    req(nrow(dat) > 0)
    dat <- dat %>%
      filter(timepoint == 8 | (start_temp == end_temp & timepoint == 0))
    if (nrow(dat) == 0) return(dat)
    if (dia_is_category_mode()) dat <- dia_add_categories(dat, dia_gene_list())
    dat
  })

  # ── Time plot data ──────────────────────────────────────────────────────
  dia_time_plot_data <- reactive({
    req(input$end_temp)
    dat <- dia_filtered_data()
    req(nrow(dat) > 0)
    et <- as.numeric(input$end_temp)
    dat <- dat %>%
      filter((end_temp == et) | (end_temp == start_temp & timepoint == 0))
    if (nrow(dat) == 0) return(dat)
    if (dia_is_category_mode()) dat <- dia_add_categories(dat, dia_gene_list())
    dat
  })

  # ── Summarize (mean + optional biorep points) ──────────────────────────
  dia_summarize_plot_data <- function(dat, x_var) {
    label_col <- if (dia_is_category_mode() && "category" %in% names(dat)) "category" else "display_label"
    mean_dat <- dat %>%
      group_by(.data[[label_col]], species, .data[[x_var]]) %>%
      summarise(sd = sd(pSup, na.rm = TRUE), pSup = mean(pSup, na.rm = TRUE),
                .groups = "drop") %>%
      mutate(pSup_lo = pmax(0, pSup - sd), pSup_hi = pmin(1, pSup + sd),
             display_label = .data[[label_col]])
    if (isTRUE(input$show_bioreps)) {
      biorep_dat <- dat %>% mutate(display_label = .data[[label_col]])
      list(mean = mean_dat, bioreps = biorep_dat)
    } else {
      list(mean = mean_dat, bioreps = NULL)
    }
  }

  # ── DIA plot builder ────────────────────────────────────────────────────
  dia_build_plot <- function(dat, x_var, x_label, x_min, x_max, x_breaks) {
    req(nrow(dat) > 0)
    plot_data <- dia_summarize_plot_data(dat, x_var)
    mean_dat <- plot_data$mean
    biorep_dat <- plot_data$bioreps

    p <- ggplot(mean_dat, aes(
      x = .data[[x_var]], y = pSup,
      colour = display_label, shape = species,
      group = interaction(display_label, species)
    )) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 3)

    if (!isTRUE(input$show_bioreps) && isTRUE(input$show_errorbars)) {
      p <- p + geom_errorbar(
        aes(ymin = pSup_lo, ymax = pSup_hi),
        width = (x_max - x_min) * 0.015, linewidth = 0.4)
    }

    if (!is.null(biorep_dat)) {
      p <- p + geom_point(
        data = biorep_dat,
        aes(x = .data[[x_var]], y = pSup,
            colour = display_label, shape = species),
        size = 2, alpha = 0.4, inherit.aes = FALSE)
    }

    p +
      scale_colour_brewer(palette = "Dark2") +
      scale_shape_manual(values = dia_species_shapes,
                         labels = dia_species_labels_md,
                         name = "Species") +
      labs(colour = NULL) +
      xlab(x_label) +
      scale_y_pSup() +
      scale_x_continuous(breaks = x_breaks, expand = expansion(mult = 0.02)) +
      coord_cartesian(xlim = c(x_min, x_max)) +
      theme(legend.position = "right", legend.text = element_markdown())
  }

  dia_build_temp_plot <- function() {
    dat <- dia_temp_plot_data()
    req(nrow(dat) > 0)
    t_min <- if (!is.null(input$temp_min)) input$temp_min else 23
    t_max <- if (!is.null(input$temp_max)) input$temp_max else 50
    dia_build_plot(dat, x_var = "end_temp",
      x_label = "Temperature (\u00b0C) after 8 min.",
      x_min = t_min, x_max = t_max,
      x_breaks = dia_temp_axis_breaks(t_min, t_max))
  }

  dia_build_time_plot <- function() {
    dat <- dia_time_plot_data()
    req(nrow(dat) > 0)
    et <- as.numeric(input$end_temp)
    ti_min <- if (!is.null(input$time_min)) input$time_min else 0
    ti_max <- if (!is.null(input$time_max)) input$time_max else 20
    time_breaks <- sort(unique(c(
      seq(ti_min, ti_max, by = 4),
      unique(dat$timepoint[dat$timepoint >= ti_min & dat$timepoint <= ti_max])
    )))
    dia_build_plot(dat, x_var = "timepoint",
      x_label = paste0("Minutes at ", et, "\u00b0C"),
      x_min = ti_min, x_max = ti_max, x_breaks = time_breaks)
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

  # ── shinyURL ─────────────────────────────────────────────────────────────
  shinyURL.server(session)
}

# ══════════════════════════════════════════════════════════════════════════════
# Run
# ══════════════════════════════════════════════════════════════════════════════

shinyApp(ui, server)
