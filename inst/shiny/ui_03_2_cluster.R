shinyPanelCluster <- fluidPage(
  tags$div(
    class = "container",
    h3("Clustering"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v05-tab03_Dimensionality-Reduction-and-Clustering.html",
              "(help)", target = "_blank")),
    wellPanel(
      # CLUSTERING --> VISUALIZATION
      fluidRow(
        column(
          6,
          selectInput("clustAlgo", "Select Algorithm",
                      c("Scran SNN"=1, "K-Means"=2),
                      selected = "Scran SNN")
        )
      ),
      h4("Input Parameters:"),
      fluidRow(
        # Scran SNN ####
        conditionalPanel(
          "input.clustAlgo == 1",
          column(
            6,
            selectInput("clustScranSNNInType", "Select Input Matrix Type:",
                        c("Assay", "ReducedDim", "AltExp")),
            conditionalPanel(
              "input.clustScranSNNInType == 'AltExp'",
              uiOutput("clustScranSNNAltExpAssayUI")
            )
          ),
          column(
            6,
            conditionalPanel(
              "input.clustScranSNNInType == 'Assay'",
              selectInput("clustScranSNNAssay", "Select An Assay:", currassays)
            ),
            conditionalPanel(
              "input.clustScranSNNInType == 'ReducedDim'",
              selectInput("clustScranSNNReddim", "Select A ReducedDim:", currreddim)
            ),
            conditionalPanel(
              "input.clustScranSNNInType == 'AltExp'",
              selectInput("clustScranSNNAltExp", "Select An AltExp:", curraltExps)
            )
          ),
          column(
            12,
            conditionalPanel(
              "input.clustScranSNNInType == 'Assay'",
              helpText("An 'assay' contains full sized data matrix with all cells and features.")
            ),
            conditionalPanel(
              "input.clustScranSNNInType == 'ReducedDim'",
              helpText("A 'reducedDim' contains low-dimension representation of an assay.\n Dimension reduction has to be run in advance.")
            ),
            conditionalPanel(
              "input.clustScranSNNInType == 'AltExp'",
              helpText("An 'altExp' contains an assay with subseted features.")
            )
          ),
          column(
            4,
            numericInput("clustScranSNNK", "K value:", 10, min = 1, step = 1),
            selectInput("clustScranSNNAlgo", "Graph Clustering Algorithm:",
                        c("walktrap", "louvain", "infomap", "fastGreedy",
                          "labelProp", "leadingEigen"), selected = "walktrap")
          ),
          conditionalPanel(
            "input.clustScranSNNInType != 'ReducedDim'",
            column(
              4,
              numericInput("clustScranSNNd", "Number of Components:",
                           50, min = 2, step = 5)
            )
          ),
          column(
            4,
            selectInput("clustScranSNNType", "Edge Weight Type:",
                        c("rank", "number", "jaccard"), selected = "rank")
          )
        ),

        # K-Means ####
        conditionalPanel(
          "input.clustAlgo == 2",
          column(
            6,
            selectInput("clustKMeansReddim", "Select A ReducedDim:", currreddim)
          ),
          column(6),
          column(
            12,
            helpText("A 'reducedDim' contains low-dimension representation of an assay.\n Dimension reduction has to be run in advance.")
          ),
          column(
            4,
            numericInput("clustKMeansN", "Number of Centers (Clusters):", value = NULL),
            selectInput("clustKMeansAlgo", "Algorithm",
                        c("Hartigan-Wong", "Lloyd", "MacQueen"),
                        selected = "Hartigan-Wong")
          ),
          column(
            4,
            numericInput("clustKMeansNIter", "Max Number of Iterations:",
                         10, min = 2, step = 1)
          ),
          column(
            4,
            numericInput("clustKMeansNStart", "Number of Random Sets:",
                         1, min = 1, step = 1)
          )
        )
      ), # fuildRow ends here
      uiOutput("clustNameUI"),
      withBusyIndicatorUI(actionButton("clustRun", "Run"))
    )
  )
)
