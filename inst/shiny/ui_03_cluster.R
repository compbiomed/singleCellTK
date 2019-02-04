shinyPanelCluster <- fluidPage(
  tags$div(
    class = "container",
    h1("Visualization & Clustering"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v05-tab03_Dimensionality-Reduction-and-Clustering.html",
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        ###  VISUALIZATION (e.g. PCA, tSNE)
        selectInput("dimRedAssaySelect", "Select Assay:", currassays),
        selectInput("dimRedPlotMethod", "Visualization Method:", c("PCA", "tSNE", "UMAP", "Dendrogram")),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'PCA'", "dimRedPlotMethod"),
          selectInput("pcX", "X axis:", pcComponents),
          selectInput("pcY", "Y axis:", pcComponents, selected = pcComponentsSelectedY)
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'PCA' || input['%s'] == 'tSNE' || input['%s'] == 'UMAP'", "dimRedPlotMethod", "dimRedPlotMethod", "dimRedPlotMethod"),
          selectInput("colorBy", "Color points by:", c("No Color", "Gene Expression", clusterChoice)),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'Gene Expression'", "colorBy"),
            #radioButtons("colorGeneBy", "Gene list:", c("Manual Input", "Biomarker (from DE tab)")), #TODO: implement biomarker color by
            radioButtons("colorGeneBy", "Gene list:", "Manual Input"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'Manual Input'", "colorGeneBy"),
              selectizeInput(
                "colorGenes", label = "Select Gene(s):", NULL, multiple = TRUE
              )
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'Biomarker (from DE tab)'", "colorGeneBy"),
              textInput("colorGenesBiomarker", "Enter Name of Gene List:", "")
            ),
            radioButtons("colorBinary", "Color scale:", c("Binary", "Continuous"), selected = "Continuous")
          ),
          selectInput("shapeBy", "Shape points by:", c("No Shape", clusterChoice)),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'tSNE'", "dimRedPlotMethod"),
            withBusyIndicatorUI(actionButton("reRunTSNE", "Re-run tSNE"))
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'PCA'", "dimRedPlotMethod"),
            withBusyIndicatorUI(actionButton("reRunPCA", "Re-run PCA"))
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'UMAP'", "dimRedPlotMethod"),
            withBusyIndicatorUI(actionButton("reRunUMAP", "Re-run UMAP"))
          )
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'PCA' || input['%s'] == 'tSNE' || input['%s'] == 'UMAP'", "dimRedPlotMethod", "dimRedPlotMethod", "dimRedPlotMethod"),
          radioButtons("booleanCluster", "Cluster Data?", c("Yes", "No"), selected = "No")
        ),
        # CLUSTERING --> VISUALIZATION
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Yes'", "booleanCluster"),
          selectInput("selectClusterInputData", "Data to Cluster:", c("PCA Components", "tSNE Components", "UMAP Components")),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'PCA' || input['%s'] == 'tSNE' || input['%s'] == 'UMAP'", "dimRedPlotMethod", "dimRedPlotMethod", "dimRedPlotMethod"),
            radioButtons("clusteringAlgorithm", "Select Clustering Algorithm:", c("K-Means", "Clara"))
          ),
          ##----------------------------------#
          # K-Means
          conditionalPanel(
            condition = sprintf("input['%s'] == 'K-Means'  && input['%s'] != 'Dendrogram'", "clusteringAlgorithm", "dimRedPlotMethod"),
            selectInput("Knumber", "Number of Clusters (k):", numClusters)
          ),
          ##----------------------------------#
          # Clara
          conditionalPanel(
            condition = sprintf("input['%s'] == 'Clara' && input['%s'] != 'Dendrogram'", "clusteringAlgorithm", "dimRedPlotMethod"),
            selectInput("Cnumber", "Number of Clusters:", numClusters)
          ),
          ##----------------------------------#
          # K-Means and Clara
          conditionalPanel(
            condition = sprintf("input['%s'] != 'Dendrogram'", "dimRedPlotMethod"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'Clara' || input['%s'] == 'K-Means'", "clusteringAlgorithm", "clusteringAlgorithm"),
              textInput("clusterName", "Name of Clusters:", value = "clusters")
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
      ),
      mainPanel(
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Dendrogram'", "dimRedPlotMethod"),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'Hierarchical' || input['%s'] == 'Phylogenetic Tree'", "clusteringAlgorithmD", "clusteringAlgorithmD"),
            plotOutput("treePlot")
          )
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'PCA' || input['%s'] == 'tSNE' || input['%s'] == 'UMAP'", "dimRedPlotMethod", "dimRedPlotMethod", "dimRedPlotMethod"),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'No' || input['%s'] == 'K-Means' || input['%s'] == 'Clara'", "booleanCluster", "clusteringAlgorithm", "clusteringAlgorithm"),
            conditionalPanel(
              condition = sprintf("input['%s'] != 'Gene Expression'", "colorBy"),
              plotlyOutput("clusterPlot", height = "600px")
            )
          )
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'PCA' || input['%s'] == 'tSNE' || input['%s'] == 'UMAP'", "dimRedPlotMethod", "dimRedPlotMethod", "dimRedPlotMethod"),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'Gene Expression'", "colorBy"),
            plotOutput("geneExpressionPlot", height = "600px")
          )
        ),
        tableOutput("pctable")
      )
    )
  )
)
