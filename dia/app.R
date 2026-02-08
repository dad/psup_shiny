# DIA Three-Species pSup Data Browser
# Shiny app for browsing protein supernatant proportion (pSup) data
# from DIA mass spectrometry across S. cerevisiae, S. kudriavzevii, and K. marxianus.

library(shiny)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(viridis)
library(ggtext)

# ── Data loading ─────────────────────────────────────────────────────────────
ps_dt <- read_tsv("data/dia_psup_data.tsv", show_col_types = FALSE)

# ── Constants ────────────────────────────────────────────────────────────────

all_species <- c("S. cerevisiae", "S. kudriavzevii", "K. marxianus")

species_shapes <- c(
  "S. cerevisiae"   = 16,  # filled circle
  "S. kudriavzevii" = 17,  # filled triangle
  "K. marxianus"    = 15   # filled square
)

# Italic markdown labels for legend (ggtext)
species_labels_md <- c(
  "S. cerevisiae"   = "*S. cerevisiae*",
  "S. kudriavzevii" = "*S. kudriavzevii*",
  "K. marxianus"    = "*K. marxianus*"
)

# Canonical temperatures where data were collected
data_temps <- c(23, 30, 37, 42, 46, 50)

# Pre-compute which base temperatures exist for each species
species_base_temps <- ps_dt %>%
  distinct(species, start_temp) %>%
  arrange(species, start_temp)

# Pre-compute full condition map
condition_map <- ps_dt %>%
  distinct(species, start_temp, end_temp, timepoint) %>%
  arrange(species, start_temp, end_temp, timepoint)

# ── Theme ────────────────────────────────────────────────────────────────────
theme_set(
  theme_minimal(base_size = 14) %+replace%
    theme(legend.position = "right")
)

# ── Utility functions ────────────────────────────────────────────────────────

# Reused from existing psup_shiny app (server.R)
splitstring <- function(s) {
  classes <- str_split(s, ";")[[1]]
  x <- sapply(classes, function(ss) {
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
  x
}

scale_y_pSup <- function() {
  list(
    scale_y_continuous(
      name = "Proportion in the\nsupernatant (pSup)",
      expand = c(0.01, 0),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1)
    ),
    theme(axis.title.y = element_text(angle = 90, vjust = 0.5))
  )
}

# Match gene query across gene, scer_ortholog, and orf columns
match_genes <- function(data, query_genes) {
  data %>%
    filter(
      toupper(gene) %in% query_genes |
        toupper(scer_ortholog) %in% query_genes |
        toupper(orf) %in% query_genes
    )
}

# Get a display label for matched rows (prefer gene name)
get_display_label <- function(data) {
  data %>%
    mutate(display_label = case_when(
      !is.na(gene) & gene != "" ~ gene,
      !is.na(scer_ortholog) & scer_ortholog != "" ~ scer_ortholog,
      TRUE ~ orf
    ))
}

# Build temperature axis breaks: always include data temps that fall in range
temp_axis_breaks <- function(t_min, t_max) {
  regular <- seq(
    floor(t_min / 5) * 5,
    ceiling(t_max / 5) * 5,
    by = 5
  )
  sort(unique(c(regular, data_temps[data_temps >= t_min & data_temps <= t_max])))
}

# ── UI ───────────────────────────────────────────────────────────────────────
ui <- fluidPage(

  titlePanel("Proportion in supernatant (pSup) across yeast species (DIA)"),

  sidebarLayout(
    sidebarPanel(
      textAreaInput(
        "ids",
        "Gene identifiers separated by semicolons:",
        value = "PGK1;PAB1;DED1",
        rows = 4
      ),
      helpText("Examples (click to update):"),
      actionLink("flat_examples", "Individual proteins"),
      actionLink("category_examples", "Protein categories"),
      helpText(" "),

      selectInput(
        "plotType", "Plot against:",
        c("temperature", "time"),
        selected = "temperature"
      ),

      # Per-species base temperature selection
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
  ),

  mainPanel(
    p(
      "Data from: Keyport Kik S, Christopher D, Glauninger H, Wong Hickernell C, Bard JAM, Lin KM, Squires AH, Ford M, Sosnick TR, Drummond DA.",
      br(),
      "\u201cAn adaptive biomolecular condensation response is conserved across environmentally divergent species.\u201d",
      br(),
      tags$em("Nature Communications"),
      " 15:3127 (2024). ",
      a("doi:10.1038/s41467-024-47355-9", href = "https://doi.org/10.1038/s41467-024-47355-9"),
      " | ",
      a("Paper page", href = "https://drummondlab.org/papers/paper/conserved-condensation")
    )
  )
)

# ── Server ───────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Dynamic UI: per-species base temperature checkboxes ──────────────────
  output$species_base_temp_ui <- renderUI({
    # Build one checkboxGroup per species showing available base temps
    ui_elements <- lapply(all_species, function(sp) {
      avail_temps <- species_base_temps %>%
        filter(species == sp) %>%
        pull(start_temp) %>%
        sort()

      if (length(avail_temps) == 0) return(NULL)

      choices <- setNames(
        paste0(sp, "|", avail_temps),
        paste0(avail_temps, "\u00b0C")
      )

      # Default: select all available base temps
      checkboxGroupInput(
        inputId = paste0("base_temp_", gsub("[^a-zA-Z]", "", sp)),
        label = tags$i(sp),
        choices = choices,
        selected = choices,
        inline = TRUE
      )
    })

    do.call(tagList, ui_elements)
  })

  # Reactive: collect all selected species|base_temp combos
  selected_species_temps <- reactive({
    combos <- c()
    for (sp in all_species) {
      input_id <- paste0("base_temp_", gsub("[^a-zA-Z]", "", sp))
      vals <- input[[input_id]]
      if (!is.null(vals)) combos <- c(combos, vals)
    }
    if (length(combos) == 0) return(tibble(species = character(), start_temp = numeric()))

    tibble(combo = combos) %>%
      separate(combo, into = c("species", "start_temp"), sep = "\\|") %>%
      mutate(start_temp = as.numeric(start_temp))
  })

  # ── Dynamic UI: end_temp (for time plots) ────────────────────────────────
  output$end_temp_ui <- renderUI({
    sel <- selected_species_temps()
    req(nrow(sel) > 0)
    avail <- condition_map %>%
      inner_join(sel, by = c("species", "start_temp")) %>%
      filter(end_temp != start_temp) %>%
      distinct(end_temp) %>%
      pull(end_temp) %>%
      sort()
    if (length(avail) == 0) avail <- c(37)
    selectInput(
      "end_temp",
      "Heat shock temperature (\u00b0C):",
      choices = avail,
      selected = max(avail)
    )
  })

  # ── Parse gene input ─────────────────────────────────────────────────────
  gene_list <- reactive({
    req(input$ids)
    splitstring(input$ids)
  })

  all_genes <- reactive({
    unique(unlist(gene_list()))
  })

  is_category_mode <- reactive({
    cats <- gene_list()
    length(cats) > 1 || any(sapply(cats, length) > 1)
  })

  # ── Data preparation: filter to selected species/base_temp combos ────────
  filtered_data <- reactive({
    sel <- selected_species_temps()
    req(nrow(sel) > 0)
    genes <- all_genes()

    dat <- ps_dt %>%
      inner_join(sel, by = c("species", "start_temp")) %>%
      match_genes(genes) %>%
      get_display_label()

    dat
  })

  # ── Data preparation: temperature plot ───────────────────────────────────
  temp_plot_data <- reactive({
    dat <- filtered_data()
    req(nrow(dat) > 0)

    # Temperature plot: x = end_temp
    # Include 8-min timepoint for heat-stressed samples,
    # and timepoint=0 for baseline (where start_temp == end_temp)
    dat <- dat %>%
      filter(
        timepoint == 8 | (start_temp == end_temp & timepoint == 0)
      )

    if (nrow(dat) == 0) return(dat)

    if (is_category_mode()) {
      dat <- add_categories(dat, gene_list())
    }

    dat
  })

  # ── Data preparation: time plot ──────────────────────────────────────────
  time_plot_data <- reactive({
    req(input$end_temp)
    dat <- filtered_data()
    req(nrow(dat) > 0)
    et <- as.numeric(input$end_temp)

    # Time plot: x = timepoint
    # Include the specific start_temp -> end_temp transition,
    # plus baseline (timepoint=0, start_temp==end_temp)
    dat <- dat %>%
      filter(
        (end_temp == et) |
          (end_temp == start_temp & timepoint == 0)
      )

    if (nrow(dat) == 0) return(dat)

    if (is_category_mode()) {
      dat <- add_categories(dat, gene_list())
    }

    dat
  })

  # ── Summarize data: compute means and optionally keep biorep points ──────
  summarize_plot_data <- function(dat, x_var) {
    label_col <- if (is_category_mode() && "category" %in% names(dat)) "category" else "display_label"

    # Always compute the mean line
    mean_dat <- dat %>%
      group_by(.data[[label_col]], species, .data[[x_var]]) %>%
      summarise(
        sd = sd(pSup, na.rm = TRUE),
        pSup = mean(pSup, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        pSup_lo = pmax(0, pSup - sd),
        pSup_hi = pmin(1, pSup + sd),
        display_label = .data[[label_col]]
      )

    if (isTRUE(input$show_bioreps)) {
      # Return both mean and raw biorep points
      biorep_dat <- dat %>%
        mutate(display_label = .data[[label_col]])
      list(mean = mean_dat, bioreps = biorep_dat)
    } else {
      list(mean = mean_dat, bioreps = NULL)
    }
  }

  # ── Category assignment helper ───────────────────────────────────────────
  add_categories <- function(dat, cats) {
    gene_to_cat <- bind_rows(lapply(names(cats), function(cat_name) {
      tibble(query_gene = cats[[cat_name]], category = cat_name)
    }))

    dat %>%
      mutate(
        query_match = case_when(
          toupper(gene) %in% gene_to_cat$query_gene ~ toupper(gene),
          toupper(scer_ortholog) %in% gene_to_cat$query_gene ~ toupper(scer_ortholog),
          toupper(orf) %in% gene_to_cat$query_gene ~ toupper(orf),
          TRUE ~ NA_character_
        )
      ) %>%
      left_join(gene_to_cat, by = c("query_match" = "query_gene")) %>%
      filter(!is.na(category)) %>%
      select(-query_match)
  }

  # ── Plot rendering ──────────────────────────────────────────────────────
  output$plot_ui <- renderUI({
    w <- paste0(input$plot_width, "px")
    h <- paste0(input$plot_height, "px")
    plotOutput("plot", width = w, height = h)
  })

  current_plot <- reactive({
    if (input$plotType == "temperature") {
      build_temp_plot()
    } else {
      build_time_plot()
    }
  })

  output$plot <- renderPlot({
    current_plot()
  })

  # ── Shared plot builder ──────────────────────────────────────────────────
  build_plot <- function(dat, x_var, x_label, x_min, x_max, x_breaks) {
    req(nrow(dat) > 0)

    label_col <- if (is_category_mode() && "category" %in% names(dat)) "category" else "display_label"
    plot_data <- summarize_plot_data(dat, x_var)
    mean_dat <- plot_data$mean
    biorep_dat <- plot_data$bioreps

    # Base plot: mean lines with color = gene/category, shape = species
    p <- ggplot(mean_dat, aes(
      x = .data[[x_var]], y = pSup,
      colour = display_label,
      shape = species,
      group = interaction(display_label, species)
    )) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 3)

    # Error bars
    if (!isTRUE(input$show_bioreps) && isTRUE(input$show_errorbars)) {
      p <- p + geom_errorbar(
        aes(ymin = pSup_lo, ymax = pSup_hi),
        width = (x_max - x_min) * 0.015,
        linewidth = 0.4
      )
    }

    # Biorep points: overlay individual data on the mean lines
    if (!is.null(biorep_dat)) {
      p <- p + geom_point(
        data = biorep_dat,
        aes(
          x = .data[[x_var]], y = pSup,
          colour = display_label,
          shape = species
        ),
        size = 2, alpha = 0.4,
        inherit.aes = FALSE
      )
    }

    # Number of distinct gene/category levels for viridis palette
    n_labels <- n_distinct(mean_dat$display_label)

    p +
      scale_colour_viridis_d(
        end = 0.9,
        option = "viridis"
      ) +
      scale_shape_manual(
        values = species_shapes,
        labels = species_labels_md,
        name = "Species"
      ) +
      labs(colour = NULL) +
      xlab(x_label) +
      scale_y_pSup() +
      scale_x_continuous(
        breaks = x_breaks,
        expand = expansion(mult = 0.02)
      ) +
      coord_cartesian(xlim = c(x_min, x_max)) +
      theme(legend.text = element_markdown())
  }

  build_temp_plot <- function() {
    dat <- temp_plot_data()
    req(nrow(dat) > 0)

    t_min <- if (!is.null(input$temp_min)) input$temp_min else 23
    t_max <- if (!is.null(input$temp_max)) input$temp_max else 50

    build_plot(
      dat,
      x_var = "end_temp",
      x_label = "Temperature (\u00b0C) after 8 min.",
      x_min = t_min,
      x_max = t_max,
      x_breaks = temp_axis_breaks(t_min, t_max)
    )
  }

  build_time_plot <- function() {
    dat <- time_plot_data()
    req(nrow(dat) > 0)

    et <- as.numeric(input$end_temp)
    ti_min <- if (!is.null(input$time_min)) input$time_min else 0
    ti_max <- if (!is.null(input$time_max)) input$time_max else 20

    # Time axis breaks: include all timepoints that exist in data
    time_breaks <- sort(unique(c(
      seq(ti_min, ti_max, by = 4),
      unique(dat$timepoint[dat$timepoint >= ti_min & dat$timepoint <= ti_max])
    )))

    build_plot(
      dat,
      x_var = "timepoint",
      x_label = paste0("Minutes at ", et, "\u00b0C"),
      x_min = ti_min,
      x_max = ti_max,
      x_breaks = time_breaks
    )
  }

  # ── Data info panel ────────────────────────────────────────────────────
  output$data_info <- renderUI({
    genes <- all_genes()
    found <- ps_dt %>%
      match_genes(genes) %>%
      distinct(gene, scer_ortholog, orf)

    found_upper <- c(
      toupper(found$gene),
      toupper(found$scer_ortholog),
      toupper(found$orf)
    )
    missing <- setdiff(genes, found_upper)

    info_parts <- c()
    if (length(missing) > 0) {
      info_parts <- c(
        info_parts,
        paste0("<b>Not found:</b> ", paste(missing, collapse = ", "))
      )
    }

    if (length(info_parts) > 0) {
      HTML(paste(info_parts, collapse = "<br>"))
    }
  })

  # ── Example links ──────────────────────────────────────────────────────
  observeEvent(input$flat_examples,
    updateTextAreaInput(session, "ids", value = "PGK1;PMA1;PAB1")
  )
  observeEvent(input$category_examples,
    updateTextAreaInput(session, "ids",
      value = "glycolytic=PGK1,FBA1,TDH3;superaggregator=DED1,NUG1,OLA1;ribosomal=RPL4A,RPL19A,RPS9B"
    )
  )

  # ── Figure download ────────────────────────────────────────────────────
  output$dl_plot <- downloadHandler(
    filename = function() {
      paste0("psup_dia_", Sys.Date(), ".", input$dl_format)
    },
    content = function(file) {
      # Convert pixel dimensions to inches at 150 dpi
      w_in <- input$plot_width / 150
      h_in <- input$plot_height / 150
      p <- current_plot()

      if (input$dl_format == "png") {
        ggsave(file, plot = p, device = "png",
               width = w_in, height = h_in, dpi = 300)
      } else if (input$dl_format == "svg") {
        ggsave(file, plot = p, device = "svg",
               width = w_in, height = h_in)
      } else if (input$dl_format == "pdf") {
        ggsave(file, plot = p, device = "pdf",
               width = w_in, height = h_in)
      }
    }
  )
}

# ── Run ──────────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
