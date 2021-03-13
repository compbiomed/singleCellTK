scDblFinderHelpModal <- function() {
  modalDialog(
    tags$style(HTML("
      div {
        word-wrap: break-word;
      }
      ")),
    tags$div(
      h3("scDblFinder Parameters"),
      fluidRow(
        column(4, tags$b("Parameter Name")),
        column(8, tags$b("Description"))
      ),
      tags$hr(),
      fluidRow(
        column(4, "nNeighbors"),
        column(8, "Integer. Number of nearest neighbors used to calculate density for doublet detection. Default 50.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "simDoublets"),
        column(8, 'Integer. Number of simulated doublets created for doublet detection. Default 10000.')
      )
    )
  )
}

