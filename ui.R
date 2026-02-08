library(shiny)
# install.packages("devtools")
# devtools::install_github("aoles/shinyURL")
library(shinyURL)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("Proportion in supernatant (pSup) for yeast genes"),

    # Sidebar with a slider input for the number of bins
    sidebarLayout(
        sidebarPanel(
            textAreaInput("ids",
                        "Gene identifiers separated by semicolons:",
                        value = "PGK1;PAB1;PMA1", rows=4),
            helpText("Examples (click to update):"),
            actionLink("flat_examples", "Individual proteins"),
            actionLink("category_examples", "Protein categories"),
            helpText(" "),
            selectInput("idType", "Identify by:", c("gene","orf"), selected="gene"),
            checkboxInput("interval", "Show 95% intervals", FALSE),
            selectInput("plotType", "Plot against:", c("time","temperature"), selected="time"),
            hr(),
            h5("Plot settings"),
            fluidRow(
                column(6, numericInput("plot_width", "Width (px):", value = 700, min = 300, max = 1600, step = 50)),
                column(6, numericInput("plot_height", "Height (px):", value = 500, min = 200, max = 1200, step = 50))
            )
        ),

        # Show a plot of the generated distribution
        mainPanel(
            uiOutput("plot_ui"),
            fluidRow(
                column(3, selectInput("dl_format", NULL,
                    choices = c("PNG" = "png", "SVG" = "svg", "PDF" = "pdf"),
                    selected = "png", width = "100%")),
                column(3, downloadButton("dl_plot", "Download figure"))
            )
        )
    ),

    # Subtitle with paper reference
    mainPanel(
        p("Data from: Wallace EWJ, Kear-Scott JL, Pilipenko EV, Schwartz MH, Laskowski PR, Rojek AE,",
          "Katanski CD, Riback JA, Dion MF, Franks AM, Airoldi EM, Pan T, Budnik BA, Drummond DA.",
          "\u201cReversible, specific, active aggregates of endogenous proteins assemble upon heat stress.\u201d",
          tags$em("Cell"),
          "162(6), 2015. ",
          a("doi:10.1016/j.cell.2015.08.041", href = "https://doi.org/10.1016/j.cell.2015.08.041"),
          " | ",
          a("Paper page", href = "https://drummondlab.org/papers/paper/endogenous-aggregates")
        )
    ),

    # Include button/text box for URL.
    shinyURL.ui()
))