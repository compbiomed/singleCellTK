doubletCellsHelpModal <- function() {
  modalDialog(
    tags$style(HTML("
      div {
        word-wrap: break-word;
      }
      ")),
    tags$div(
      h3("Doublet Cells Parameters"),
      fluidRow(
        column(4, tags$b("Parameter Name")),
        column(8, tags$b("Description"))
      ),
      tags$hr(),
      fluidRow(
        column(4, "nNeighbors"),
        column(8, "Number of nearest neighbors used to calculate density for doublet detection. Default 50.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "simDoublets"),
        column(8, 'Number of simulated doublets created for doublet detection. Default 10000.')
      ),
    )
  )
}

QCMHelpModal <- function() {
  modalDialog(
    tags$style(HTML("
      div {
        word-wrap: break-word;
      }
      ")),
    tags$div(
      h3("QC Metrics Parameters"),
      fluidRow(
        column(4, tags$b("Parameter Name")),
        column(8, tags$b("Description"))
      ),
      tags$hr(),
      fluidRow(
        column(4, "collectionName"),
        column(8, "Character. Name of a GeneSetCollection obtained by using one of the importGeneSet* functions. Default NULL.")
      ),
    )
  )
}

