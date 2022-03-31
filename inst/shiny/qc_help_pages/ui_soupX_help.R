soupXHelpModal <- function() {
    modalDialog(
        tags$style(HTML("
      div {
        word-wrap: break-word;
      }
      ")),
        tags$div(
            h3("SoupX Parameters"),
            fluidRow(
                column(4, tags$b("Parameter Name")),
                column(8, tags$b("Description"))
            ),
            tags$hr(),
            fluidRow(
                column(4, "cluster"),
                column(8, "Prior knowledge of clustering labels on cells. Can 
                       be specified with preloaded cell annotations. When not 
                       supplied, scran::quickCluster method will be applied.")
            ),
            tags$hr(),
            fluidRow(
                column(4, "tfidfMin"),
                column(8, "Minimum value of tfidf to accept for a marker gene. 
                       Default 1. Parameter for SoupX::autoEstCont()")
            ),
            tags$hr(),
            fluidRow(
                column(4, "soupQuantile"),
                column(8, "Only use genes that are at or above this expression 
                       quantile in the soup. This prevents inaccurate estimates 
                       due to using genes with poorly constrained contribution 
                       to the background. Default 0.9. Parameter for 
                       SoupX::autoEstCont()")
            ),
            tags$hr(),
            fluidRow(
                column(4, "maxMarkers"),
                column(8, "If we have heaps of good markers, keep only the best 
                       maxMarkers of them. Default 100. Parameter for 
                       SoupX::autoEstCont()")
            ),
            tags$hr(),
            fluidRow(
                column(4, "contaminationRange"),
                column(8, "This constrains the contamination fraction to lie 
                       within this range. Must be between 0 and 1. The high end 
                       of this range is passed to 
                       SoupX::estimateNonExpressingCells as 
                       `maximumContamination`. Use 'Lower range' and 'Higher 
                       range' for the two values. Default 0.01 - 0.8. Parameter 
                       for SoupX::autoEstCont()")
            ),
            tags$hr(),
            fluidRow(
                column(4, "rhoMaxFDR"),
                column(8, "Integer. False discovery rate passed to 
                       SoupX::estimateNonExpressingCells, to test if rho is 
                       less than `maximumContamination`. Default 0.2. Parameter 
                       for SoupX::autoEstCont()")
            ),
            tags$hr(),
            fluidRow(
                column(4, "priorRho"),
                column(8, "Mode of gamma distribution prior on contamination 
                       fraction. Default 0.05. Parameter for 
                       SoupX::autoEstCont()")
            ),
            tags$hr(),
            fluidRow(
                column(4, "priorRhoStdDev"),
                column(8, "Standard deviation of gamma distribution prior on 
                       contamination fraction. Default 0.1. Parameter for 
                       SoupX::autoEstCont()")
            ),
            tags$hr(),
            fluidRow(
                column(4, "forceAccept"),
                column(8, "Should we allow very high contamination fractions to 
                       be used. Default FALSE (not checked). Parameter for 
                       SoupX::autoEstCont()")
            ),
            tags$hr(),
            fluidRow(
                column(4, "adjustMethod"),
                column(8, "Method to use for correction. One of 'subtraction', 
                       'soupOnly', or 'multinomial'. Default 'subtraction'. 
                       Parameter for SoupX::adjustCounts")
            ),
            tags$hr(),
            fluidRow(
                column(4, "roundToInt"),
                column(8, "Should the resulting matrix be rounded to integers? 
                       Default FALSE (not checked). Parameter for 
                       SoupX::adjustCounts")
            ),
            tags$hr(),
            fluidRow(
                column(4, "tol"),
                column(8, "Allowed deviation from expected number of soup 
                       counts. Don't change this. Default 0.001. Parameter for 
                       SoupX::adjustCounts")
            ),
            tags$hr(),
            fluidRow(
                column(4, "pCut"),
                column(8, "The p-value cut-off used when method = 'soupOnly'. 
                       Default 0.01. Parameter for SoupX::adjustCounts")
            )
        )
    )
}

