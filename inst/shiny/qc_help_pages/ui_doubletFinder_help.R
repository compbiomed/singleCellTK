doubletFinderHelpModal <- function() {
  modalDialog(
    tags$style(HTML("
      div {
        word-wrap: break-word;
      }
      ")),
    tags$div(
      h3("Doublet Finder Parameters"),
      fluidRow(
        column(4, tags$b("Parameter Name")),
        column(8, tags$b("Description"))
      ),
      tags$hr(),
      fluidRow(
        column(4, "seuratNfeatures"),
        column(8, "Integer. Number of highly variable genes to use. Default 2000.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "seuratRes"),
        column(8, 'Numeric vector. The resolution parameter used in seurat, which adjusts the number of 
               clusters determined via the algorithm. Default 1.5.')
      ),
      tags$hr(),
      fluidRow(
        column(4, "formationRate"),
        column(8, "Doublet formation rate used within algorithm. Default 0.075.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "seuratPcs"),
        column(8, "Numeric vector. The PCs used in seurat function to determine number of clusters. Default 1:15.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "verbose"),
        column(8, "Check off to print messages from Seurat and DoubletFinder in the console.")
      ),
    )
  )
}

