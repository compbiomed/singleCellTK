decontXHelpModal <- function() {
  modalDialog(
    tags$style(HTML("
      div {
        word-wrap: break-word;
      }
      ")),
    tags$div(
      h3("DecontX Parameters"),
      fluidRow(
        column(4, tags$b("Parameter Name")),
        column(8, tags$b("Description"))
      ),
      tags$hr(),
      fluidRow(
        column(4, "maxIter"),
        column(8, "Integer. Maximum iterations of the EM algorithm. Default 500.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "nativePrior"),
        column(8, "Integer. The first element of a vector containing the concentration parameters 
               for the Dirichlet prior for the contamination in each cell.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "contaminationPrior"),
        column(8, "Integer. The second element of a vector containing the concentration parameters 
               for the Dirichlet prior for the contamination in each cell.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "convergence"),
        column(8, "Numeric. The EM algorithm will be stopped if the maximum difference in the contamination 
               estimates between the previous and current iterations is less than this. Default 0.001.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "iterLogLik"),
        column(8, "Integer. Calculate log likelihood every iterLogLik iteration. Default 10.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "varGenes"),
        column(8, "Integer. The number of variable genes to use in dimensionality reduction before clustering. 
               Variability is calcualted using modelGeneVar function from the 'scran' package. Used only when z is not provided. 
               Default 5000.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "dbscanEps"),
        column(8, "Numeric. The clustering resolution parameter used in 'dbscan' to estimate broad cell clusters. 
               Used only when z is not provided. Default 1.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "estimateDelta"),
        column(8, "Check off to update delta at each iteration.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "verbose"),
        column(8, "Check off to print log messages. Default TRUE.")
      ),
    )
  )
}

