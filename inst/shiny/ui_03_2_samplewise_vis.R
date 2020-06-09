shinyPanelCluster <- fluidPage(
  tags$div(
    class = "container",
    h3("Samplewise Visualization and Clustering"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v05-tab03_Dimensionality-Reduction-and-Clustering.html",
              "(help)", target = "_blank")),
      # SHINYJS COLLAPSE --------------------------
      # Section 1 - Assay Settings
      actionButton("c_button1", "Run New Dimensional Reduction"),
      # open by default
      tags$div(id = "c_collapse1",
               wellPanel(
                 fluidRow(
                   column(8,
                          wellPanel(
                            fluidRow(
                              column(6,
                                     tags$h4("Select:"),
                                     selectInput("dimRedAssaySelect", "Assay:", currassays),
                                     # Note: Removed "Dendrogram" option from method select to disable conditionalPanels.
                                     selectInput("dimRedPlotMethod", "Method:", c("PCA", "tSNE", "UMAP")),
                                     tags$br(),
                                     withBusyIndicatorUI(actionButton("runDimred", "Run"))
                              ),
                              column(6,
                                     tags$h4("DR Options:"),
                                     textInput("dimRedNameInput", "reducedDim Name:", ""),
                                     tags$br(),
                                     HTML('<button type="button" class="btn btn-default btn-block"
                                          data-toggle="collapse" data-target="#c-collapse-run-options">
                                          View More Options</button>'
                                     ),
                                     tags$div(
                                       id = "c-collapse-run-options", class = "collapse",
                                       conditionalPanel(
                                         condition = sprintf("input['%s'] == 'UMAP'", "dimRedPlotMethod"),
                                         sliderInput("iterUMAP", "# of iterations", min = 50, max = 500, value = 200),
                                         sliderInput("neighborsUMAP", "# of nearest neighbors", min = 2, max = 100, value = 5),
                                         numericInput("alphaUMAP", "learning rate(alpha)", value = 1)
                                       ),
                                       conditionalPanel(
                                         condition = sprintf("input['%s'] == 'tSNE'", "dimRedPlotMethod"),
                                         sliderInput("iterTSNE", "# of iterations", min = 100, max = 2000, value = 1000),
                                         sliderInput("perplexityTSNE", "Perplexity paramter", min = 5, max = 50, value = 5)
                                       )
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
                                   selectInput("delRedDimType", label = NULL, currreddim)
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
      # Section 2 - Visualize
      actionButton("c_button2", "Visualize"),
      # open by default
      tags$div(id = "c_collapse2",
               wellPanel(
                 sidebarLayout(
                   sidebarPanel(
                     tags$h3("Visualization Settings:"),
                     ## NOT LINKED UP
                     uiOutput("usingReducedDims"),
                     uiOutput("dimRedAxisSettings"),
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
                       # conditionalPanel(
                       #   condition = sprintf("input['%s'] == 'Biomarker (from DE tab)'", "colorGeneBy"),
                       #   textInput("colorGenesBiomarker", "Enter Name of Gene List:", "")
                       # ),
                       radioButtons("colorBinary", "Color scale:", c("Binary", "Continuous"), selected = "Continuous")
                     ),
                     selectInput("shapeBy", "Shape points by:", c("No Shape", clusterChoice)),
                     ## NOT LINKED UP
                     ##selectInput("sizeBy", "Size points by:", c("No Sizing", clusterChoice)),
                     checkboxInput("axisNames", "Customise axis labels?"),
                     conditionalPanel(
                       condition = "input.axisNames == true",
                       textInput("dimRedAxis1", "Axis Name 1:", ""),
                       textInput("dimRedAxis2", "Axis Name 2:", "")
                    ),
                     tags$hr(),
                     withBusyIndicatorUI(actionButton("cUpdatePlot", "Update Plot"))
                   ),
                   mainPanel(
                     conditionalPanel(
                       condition = sprintf("input['%s'] != 'Gene Expression'", "colorBy"),
                       plotlyOutput("clusterPlot", height = "600px")
                     ),
                     conditionalPanel(
                       condition = sprintf("input['%s'] == 'Gene Expression'", "colorBy"),
                       plotOutput("geneExpPlot", height = "600px")
                     ),
                     #plotlyOutput("clusterPlot", height = "600px"),
                     tags$hr(),
                     uiOutput("pctable")
                   )
                 )
               )
      ),
      # Section 3 - Clustering
      actionButton("c_button3", "Clustering (optional)"),
      # open by default
      tags$div(id = "c_collapse3",
               wellPanel(
                 # CLUSTERING --> VISUALIZATION
                 fluidRow(
                   column(6,
                          tags$h4("Data to Cluster:"),
                          selectInput("selectClusterInputData", label = NULL, c("PCA Components", "tSNE Components", "UMAP Components")),
                          conditionalPanel(
                            condition = sprintf("input['%s'] == 'PCA' || input['%s'] == 'tSNE' || input['%s'] == 'UMAP' ", "dimRedPlotMethod", "dimRedPlotMethod", "dimRedPlotMethod"),
                            tags$h4("Select Clustering Algorithm:"),
                            radioButtons("clusteringAlgorithm", label = NULL, c("K-Means", "Clara"))
                          )
                   ),
                   column(6,
                          ##----------------------------------#
                          # K-Means
                          conditionalPanel(
                            condition = sprintf("input['%s'] == 'K-Means'", "clusteringAlgorithm"),
                            tags$h4("Number of Clusters (k):"),
                            selectInput("Knumber", label = NULL, numClusters)
                          ),
                          ##----------------------------------#
                          # Clara
                          conditionalPanel(
                            condition = sprintf("input['%s'] == 'Clara'", "clusteringAlgorithm"),
                            tags$h4("Number of Clusters:"),
                            selectInput("Cnumber", label = NULL, numClusters)
                          ),
                          ##----------------------------------#
                          # K-Means and Clara
                          conditionalPanel(
                            condition = sprintf("input['%s'] == 'Clara' || input['%s'] == 'K-Means'", "clusteringAlgorithm", "clusteringAlgorithm"),
                            tags$h4("Name of Clusters:"),
                            textInput("clusterName", label = NULL, value = "clusters")
                          )
                          ##----------------------------------#
                          ## Input other clustering algorithms here
                          ##----------------------------------#
                   )
                 ),
                 conditionalPanel(
                   condition = sprintf("input['%s'] == 'K-Means' || input['%s'] == 'Clara'", "clusteringAlgorithm", "clusteringAlgorithm"),
                   withBusyIndicatorUI(actionButton("clusterData", "Cluster Data")),
                   helpText("To visualize this, run any clustering method here -> scroll up (visualize) and select your saved cluster name in plot by shape/color.")
                 )
               )
              )
    )
)
