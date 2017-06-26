shiny_panel_cluster <- fluidPage(
  tabsetPanel(
    tabPanel(
      "Clustering",
      fluidRow(
        column(
          4,
          wellPanel(
            ###  VISUALIZATION (e.g. PCA, tSNE)
            selectInput("dimRedPlotMethod", "Visualization Method:", c("PCA", "tSNE", "Dendrogram")),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'PCA'", "dimRedPlotMethod"),
              selectInput("pcX", "X axis:", pcComponents),
              selectInput("pcY", "Y axis:", pcComponents, selected = "PC2")
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'PCA' || input['%s'] == 'tSNE'", "dimRedPlotMethod", "dimRedPlotMethod"),
              selectInput("colorBy", "Color points by:", c("No Color", "Gene Expression", clusterChoice)),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'Gene Expression'", "colorBy"),
                radioButtons("colorGeneBy", "Gene list:", c("Manual Input", "Biomarker (from DE tab)")),
                conditionalPanel(
                  condition = sprintf("input['%s'] == 'Manual Input'", "colorGeneBy"),
                  selectInput(
                    "colorGenes", "Select Gene(s):", geneChoice, multiple = TRUE
                  )
                ),
                conditionalPanel(
                  condition = sprintf("input['%s'] == 'Biomarker (from DE tab)'", "colorGeneBy"),
                  textInput("colorGenesBiomarker", "Enter Name of Gene List:", "")
                ),
                radioButtons("colorBinary", "Color scale:", c("Binary", "Continuous"))
              ),
              selectInput("shapeBy", "Shape points by:", c("No Shape", clusterChoice)),
              withBusyIndicatorUI(actionButton("reRunTSNE", "Re-run tSNE"))
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'PCA' || input['%s'] == 'tSNE'", "dimRedPlotMethod", "dimRedPlotMethod"),
              radioButtons("booleanCluster", "Cluster Data?", c("Yes", "No"), selected = "No")
            ),
            ### CLUSTERING --> VISUALIZATION
            conditionalPanel(
              condition = sprintf("input['%s'] == 'Yes'", "booleanCluster"),
              selectInput("selectClusterInputData", "Data to Cluster:", c("Raw Data", "PCA Components", "tSNE Components")),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'PCA' || input['%s'] == 'tSNE' ", "dimRedPlotMethod", "dimRedPlotMethod"),
                radioButtons("clusteringAlgorithm", "Select Clustering Algorithm:", c("K-Means", "Clara"))
              ),
              ##----------------------------------#
              ## K-Means
              conditionalPanel(
                condition = sprintf("input['%s'] == 'K-Means'  && input['%s'] != 'Dendrogram'", "clusteringAlgorithm", "dimRedPlotMethod"),
                selectInput("Knumber", "Number of Clusters (k):", numClusters)
              ),
              ##----------------------------------#
              ## Clara
              conditionalPanel(
                condition = sprintf("input['%s'] == 'Clara' && input['%s'] != 'Dendrogram'", "clusteringAlgorithm", "dimRedPlotMethod"),
                selectInput("Cnumber", "Number of Clusters:", numClusters)
              ),
              ##----------------------------------#
              ## K-Means & Clara
              conditionalPanel(
                condition = sprintf("input['%s'] != 'Dendrogram'", "dimRedPlotMethod"),
                conditionalPanel(
                  condition = sprintf("input['%s'] == 'Clara' || input['%s'] == 'K-Means'", "clusteringAlgorithm", "clusteringAlgorithm"),
                  textInput("clusterName", "Name of Clusters:", value = "")
                )
              ),
              ##----------------------------------#
              ## Input other clustering algorithms here
              ##----------------------------------#
              conditionalPanel(
                condition = sprintf("input['%s'] != 'Dendrogram' && input['%s'] == 'K-Means' || input['%s'] == 'Clara'", "dimRedPlotMethod", "clusteringAlgorithm", "clusteringAlgorithm"),
                withBusyIndicatorUI(actionButton("clusterData", "Cluster Data"))
              )
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'Dendrogram'", "dimRedPlotMethod"),
              radioButtons("clusteringAlgorithmD", "Select Clustering Algorithm:", c("Hierarchical", "Phylogenetic Tree"), selected = "Hierarchical"),
              selectInput("dendroDistanceMetric", "Select Distance Metric:", c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))
            )
          )
        ),
        column(8,
          conditionalPanel(
            condition = sprintf("input['%s'] == 'Dendrogram'", "dimRedPlotMethod"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'Hierarchical' || input['%s'] == 'Phylogenetic Tree'", "clusteringAlgorithmD", "clusteringAlgorithmD"),
              plotOutput("treePlot")
            )
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'PCA' || input['%s'] == 'tSNE'", "dimRedPlotMethod", "dimRedPlotMethod"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'No' || input['%s'] == 'K-Means' || input['%s'] == 'Clara'", "booleanCluster", "clusteringAlgorithm", "clusteringAlgorithm"),
              conditionalPanel(
                condition = sprintf("input['%s'] != 'Gene Expression'", "colorBy"),
                plotlyOutput("clusterPlot")
              )
            )
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'PCA' || input['%s'] == 'tSNE'", "dimRedPlotMethod", "dimRedPlotMethod"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'Gene Expression'", "colorBy"),
              plotOutput("geneExpressionPlot")
            )
          )
        ),
        column(6, ""),
        column(6, tableOutput("pctable"))
      ),
      mainPanel(
        h1("Instructions"),
        h3("Visualizing Samples by PCA or tSNE"),
        p("1. (For PCA) Choose which prinicipal components to plot"),
        p("2. Choose what option to color and/or shape your data by"),
        h4("Color by Gene Expression"),
        p("a. Choose to color by a manually inputted list of gene(s) (up to 9) or a biomarker lit from the differential expression tab"),
        p("b. Select either a continuous or binary color scale"),
        h3("Cluster with K-Means or Clara"),
        p("1. Select cluster data"),
        p("2. Choose which data to cluster (raw or dimension reduced data)"),
        p("3. Choose the clustering algorithm"),
        p("4. Choose the number of cluster centers"),
        p("5. Assign a name to your clusters"),
        p("6. Press Cluster Data"),
        h3("Visualizing Samples by Dendrogram"),
        p("1. Select visualization method to be Dendrogram)"),
        p("2. Choose which data to cluster (raw or dimension reduced data)"),
        p("3. Choose the clustering algorithm"),
        p("4. Choose the distance metric")
      )
    ),
    tabPanel(
      "PCA (MAST example)",
      fluidRow(
        column(
          4,
          wellPanel(
            selectInput(inputId = "pcaAlgorithm", label = "Algorithm", choices = c("regular PCA", "randomized PCA", "robust PCA", "RR PCA")),
            checkboxGroupInput(inputId = "pcaCheckbox", label = "Involving variables", choices = c("PC1", "PC2", "PC3"), selected = c("PC1", "PC2", "PC3")),
            selectInput(inputId = "selectAdditionalVariables", label = "Additional variables:", choices = clusterChoice, multiple = TRUE),
            selectInput("plotTypeId", "Plot Type", c("Paired Plot", "Single Plot"), "Paired Plot"),
            selectInput("colorClusters_MAST", "Color Clusters By", clusterChoice),
            withBusyIndicatorUI(actionButton("plotPCA", "Plot PCA Data"))
          )#wellPanel
        ),
        column(
          8,
          plotOutput("pcaPlot")
        )
      ),
      mainPanel(
        h1("Instructions"),
        p(""), strong("Principle Component Analysis:"),
        p("1. Choose algorithm (one of the PCA algorithms)"),
        p("2. Choose components to be involved in paired plot"),
        p("3. Choose 2 components in Checkbox to produce single plot(if more than 2 components are selected, only the first 2 will be used)"),
        p("3. Choose feature to color data by"),
        p("4. Visualize your data")
      )
    )
  ),
  includeHTML("www/footer.html")
)
