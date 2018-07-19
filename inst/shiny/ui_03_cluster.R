shinyPanelCluster <- fluidPage(
  tags$div(
    class = "container",
    h1("Visualization & Clustering"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v05-tab03_Dimensionality-Reduction-and-Clustering.html",
              "(help)", target = "_blank")),
    
    tabsetPanel(
      tabPanel(
        "Dimensionality Reduction",
        wellPanel(
          
          # ---- Run New Dimensional Reduction
          HTML('<button type="button" class="btn btn-default btn-block" 
            data-toggle="collapse" data-target="#c-collapse-run">
            Run New Dimensional Reduction</button>'),
          tags$div( 
            id="c-collapse-run", class="collapse",
            wellPanel(
              fluidRow(
                column(8,
                  wellPanel(  
                    fluidRow(
                      column(6,
                        tags$h4("Select:"),
                        selectInput("dimRedAssaySelect", "Assay:", currassays),
                        # Note: Removed "Dendrogram" option from method select to disable conditionalPanels.
                        selectInput("dimRedPlotMethod", "Method:", c("PCA", "tSNE")),
                        tags$br(),
                        ## BUTTONS NEED REPLACING:
                        ## these are just the old "re-run" buttons moved up & re-labeled to look like 1 button
                        tags$br(),
                        conditionalPanel(
                          condition = sprintf("input['%s'] == 'tSNE'", "dimRedPlotMethod"),
                          withBusyIndicatorUI(actionButton("reRunTSNE", "Run"))
                        ),
                        conditionalPanel(
                          condition = sprintf("input['%s'] == 'PCA'", "dimRedPlotMethod"),
                          withBusyIndicatorUI(actionButton("reRunPCA", "Run"))
                        )
                      ),
                      column(6,
                        tags$h4("DR Options:"),
                        ## NOT LINKED UP
                        textInput("dimRedNameInput", "reducedDim Name:", ""),
                        tags$br(),
                        HTML('<button type="button" class="btn btn-default btn-block" 
                          data-toggle="collapse" data-target="#c-collapse-run-options">
                          View More Options</button>'
                        ),
                        tags$div(
                          id="c-collapse-run-options", class="collapse",
                          tags$p("Content coming soon.")
                        )
                      )
                    )
                  )
                ),
                column(4,
                  wellPanel(
                    h4("Available Reduced Dims:"),
                    tableOutput("reducedDimsList"),
                    tags$hr(), 
                    h4("Remove a reducedDim:"),
                    fluidRow(
                      column(8,
                        selectInput("delRedDimType", label=NULL, currreddim)
                      ),
                      column(4,
                        withBusyIndicatorUI(actionButton("delRedDim", "Delete"))
                      )
                    )
                  )
                )
              )
            )
          ),
          
          # ---- Visualize
          HTML('<button type="button" class="btn btn-default btn-block" 
            data-toggle="collapse" data-target="#c-collapse-visualize">
            Visualize</button>'),
          tags$div(
            id="c-collapse-visualize", class="collapse",
            wellPanel(
              sidebarLayout(
                sidebarPanel(
                  tags$h3("Visualization Settings:"),  
                  ## NOT LINKED UP
                  selectInput("usingReducedDims", "Select Reduced Dimension Data:", ""),
                  tags$h4("Axis Settings"),
                  selectInput("pcX", "X Axis:", pcComponents),
                  selectInput("pcY", "Y Axis:", pcComponents, selected = pcComponentsSelectedY),
                
                  tags$h4("Style Settings"),
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
                  ## NOT LINKED UP
                  selectInput("sizeBy", "Size points by:", c("No Sizing", clusterChoice)),
                  tags$hr(),
                  ## NOT LINKED UP
                  withBusyIndicatorUI(actionButton("cUpdatePlot", "Update Plot"))
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
                    condition = sprintf("input['%s'] == 'PCA' || input['%s'] == 'tSNE'", "dimRedPlotMethod", "dimRedPlotMethod"),
                    conditionalPanel(
                      condition = sprintf("input['%s'] == 'No' || input['%s'] == 'K-Means' || input['%s'] == 'Clara'", "booleanCluster", "clusteringAlgorithm", "clusteringAlgorithm"),
                      conditionalPanel(
                        condition = sprintf("input['%s'] != 'Gene Expression'", "colorBy"),
                        plotlyOutput("clusterPlot", height = "600px")
                      )
                    )
                  ),
                  conditionalPanel(
                    condition = sprintf("input['%s'] == 'PCA' || input['%s'] == 'tSNE'", "dimRedPlotMethod", "dimRedPlotMethod"),
                    conditionalPanel(
                      condition = sprintf("input['%s'] == 'Gene Expression'", "colorBy"),
                      plotOutput("geneExpressionPlot", height = "600px")
                    )
                  ),
                  tags$hr(),
                  tags$h4("PC Table:"),
                  tableOutput("pctable")
                ) 
              ) 
            )
          ), 
          
          # ---- Clustering
          HTML('<button type="button" class="btn btn-default btn-block" 
            data-toggle="collapse" data-target="#c-collapse-cluster">
            Clustering (optional)</button>'),
          tags$div(
            id="c-collapse-cluster", class="collapse",
            wellPanel(
              # CLUSTERING --> VISUALIZATION
              fluidRow(
                column(6,
                  tags$h4("Data to Cluster:"),
                  selectInput("selectClusterInputData", label=NULL, c("PCA Components", "tSNE Components")),
                  conditionalPanel(
                    condition = sprintf("input['%s'] == 'PCA' || input['%s'] == 'tSNE' ", "dimRedPlotMethod", "dimRedPlotMethod"),
                    tags$h4("Select Clustering Algorithm:"),
                    radioButtons("clusteringAlgorithm", label=NULL, c("K-Means", "Clara"))
                  )
                ),
                column(6,
                  ##----------------------------------#
                  # K-Means
                  conditionalPanel(
                    condition = sprintf("input['%s'] == 'K-Means'  && input['%s'] != 'Dendrogram'", "clusteringAlgorithm", "dimRedPlotMethod"),
                    tags$h4("Number of Clusters (k):"),
                    selectInput("Knumber", label=NULL, numClusters)
                  ),
                  ##----------------------------------#
                  # Clara
                  conditionalPanel(
                    condition = sprintf("input['%s'] == 'Clara' && input['%s'] != 'Dendrogram'", "clusteringAlgorithm", "dimRedPlotMethod"),
                    tags$h4("Number of Clusters:"),
                    selectInput("Cnumber", label=NULL, numClusters)
                  ),
                  ##----------------------------------#
                  # K-Means and Clara
                  conditionalPanel(
                    condition = sprintf("input['%s'] != 'Dendrogram'", "dimRedPlotMethod"),
                    conditionalPanel(
                      condition = sprintf("input['%s'] == 'Clara' || input['%s'] == 'K-Means'", "clusteringAlgorithm", "clusteringAlgorithm"),
                      tags$h4("Name of Clusters:"),
                      textInput("clusterName", label=NULL, value = "clusters")
                    )
                  )
                  ##----------------------------------#
                  ## Input other clustering algorithms here
                  ##----------------------------------#
                )
              ),
              conditionalPanel(
                condition = sprintf("input['%s'] != 'Dendrogram' && input['%s'] == 'K-Means' || input['%s'] == 'Clara'", "dimRedPlotMethod", "clusteringAlgorithm", "clusteringAlgorithm"),
                withBusyIndicatorUI(actionButton("clusterData", "Cluster Data"))
              ),
                
              conditionalPanel(
                condition = sprintf("input['%s'] == 'Dendrogram'", "dimRedPlotMethod"),
                radioButtons("clusteringAlgorithmD", "Select Clustering Algorithm:", c("Hierarchical", "Phylogenetic Tree"), selected = "Hierarchical"),
                selectInput("dendroDistanceMetric", "Select Distance Metric:", c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))
              )
            )
          )
        )
      ),
      
      tabPanel(
        '"Visualize"',
        wellPanel(
          sidebarLayout(
            sidebarPanel(
              h4("Visualization Options:"),
              selectInput("visAssaySelect", "Select Assay:", currassays),
              selectInput("visPlotMethod", "Visualization Method:", c("boxplot", "scatterplot", "barplot", "heatmap")),
              selectInput("visCondn", "Condition:", c(clusterChoice)),
              selectizeInput("selectvisGenes", label = "Select Gene(s):", NULL, multiple = TRUE),
              withBusyIndicatorUI(actionButton("plotvis", "Plot"))
            ),
            mainPanel(
              fluidRow(
                plotOutput("visPlot")
              )
            )
          )
        )
      ),
      
      tabPanel(
        "Dendrogram",
        wellPanel(
          tags$br(),
          tags$p("Work in progress - Need to move dendrogram functionality over from DR tab.")
        )
      )
    )
  )
)
