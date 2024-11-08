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
                        "Enter gene identifiers separated by commas:",
                        value = "PGK1,PAB1,PMA1", rows=4),
            selectInput("idType", "identify by", c("gene","orf"), selected="gene"),
            checkboxInput("interval", "Show 95% intervals", FALSE),
            selectInput("plotType", "plot against", c("time","temperature"), selected="time")
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("plot") # ,plotOutput("tempPlot")
        )
    ),
    
    # Subtitle with paper reference
    mainPanel(
        p("Data from: 
          Reversible, specific, active aggregates of endogenous proteins assemble upon heat stress,
          Wallace et al., Cell 162 (6), 2015, http://drummondlab.org/endogenous-aggregates")
    ),
    
    # Include button/text box for URL.
    shinyURL.ui()
))