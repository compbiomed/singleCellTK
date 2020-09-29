shinyPanelCluster <- fluidPage(
  tags$div(
    class = "container",
    h1("Visualization & Clustering"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v05-tab03_Dimensionality-Reduction-and-Clustering.html",
              "(help)", target = "_blank")),
    tabsetPanel(
      tabPanel(
        "Gene Visualization",
        wellPanel(
          br(),
          fluidRow(
            sidebarLayout(
              sidebarPanel(
                h4("Visualization Options:"),
                selectInput("visAssaySelect", "Select Assay:", currassays),
                selectInput("visPlotMethod", "Visualization Method:", c("boxplot", "scatterplot", "barplot", "heatmap")),
                selectInput("visCondn", "Condition:", c("none", clusterChoice)),
                h3("Choose data source:"),
                radioButtons(
                  "visGeneList", label = NULL, c("Select Gene(s)" = "selVisRadioGenes",
                                                 "Saved top genes" = "visBiomarker")
                ),
                conditionalPanel(
                  condition = sprintf("input['%s'] == 'selVisRadioGenes'", "visGeneList"),
                  selectizeInput("selectvisGenes", label = "Select Gene(s):", NULL, multiple = TRUE)
                ),
                conditionalPanel(
                  helpText("To use this, first run Differential expression and save top genes."),
                  condition = sprintf("input['%s'] == 'visBiomarker'", "visGeneList"),
                  uiOutput("visBioGenes")
                ),
                uiOutput("visOptions"),
                withBusyIndicatorUI(actionButton("plotvis", "Plot"))
              ),
              mainPanel(
                fluidRow(
                  plotOutput("visPlot", height = '600px')
                )
              )
            )
          )
        )
      ),
      tabPanel(
        "Dimensionality Reduction",
        wellPanel(
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
                        ## NOT LINKED UP -- LINKED NOW
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
                            sliderInput("iterUMAP", "# of iterations", min = 50, max = 500, value = 100),
                            sliderInput("neighborsUMAP", "# of nearest neighbors", min = 2, max = 100, value = 5),
                            sliderInput("mindistUMAP", "minimum distance between points", min = 0.001, max = 0.1, value = 0.01),
                            numericInput("alphaUMAP", "learning rate(alpha)", value = 1)
                          ),
                          conditionalPanel(
                            condition = sprintf("input['%s'] == 'tSNE'", "dimRedPlotMethod"),
                            sliderInput("iterTSNE", "# of iterations", min = 100, max = 2000, value = 500),
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
                  #selectInput("usingReducedDims", "Select Reduced Dimension Data:", currreddim),
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
                  ##selectInput("sizeBy", "Size points by:", c("No Sizing", clusterChoice)),
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
                withBusyIndicatorUI(actionButton("clusterData", "Cluster Data"))
              )
            )
          )
        )
      ),
      tabPanel(
        title = "Dendrogram",
        wellPanel(
          fluidRow(
            sidebarLayout(
              sidebarPanel(
                tags$h4("Options:"),
                uiOutput("dendroRedDim"),
                radioButtons("clusteringAlgorithmD", "Select Clustering Algorithm:", c("Hierarchical", "Phylogenetic Tree"), selected = "Hierarchical"),
                selectInput("dendroDistanceMetric", "Select Distance Metric:", c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")),
                withBusyIndicatorUI(actionButton("dendroPlot", "Plot"))
              ),
              mainPanel(
                fluidRow(
                  plotOutput("treePlot", height = '600px')
                )
              )
            )
          )
        )
      )
    )
  )
)



# #Server Code (Commenting it out for now) #Irzam
# #Dendrogram
# output$dendroRedDim <- renderUI({
#   selectInput("dendroRedDim", "Select Reduced Dimension Data:", names(reducedDims(vals$counts)))
# })
# 
# observeEvent(input$dendroPlot, {
#   withBusyIndicatorServer(input$dendroPlot, {
#     if (is.null(vals$counts)){
#       shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
#     } else {
#       if (grepl(pattern = "PCA_", x = input$usingReducedDims)) {
#         comp <- "PCA Components"
#       } else  if (grepl(pattern = "TSNE_", x = input$usingReducedDims)) {
#         comp <- "TSNE Components"
#       } else {
#         comp <- "UMAP Components"
#       }
#       data <- getClusterInputData(inSCE = vals$counts,
#                                   inputData = comp,
#                                   useAssay = input$dimRedAssaySelect,
#                                   reducedDimName = input$usingReducedDims)
#       d <- stats::dist(data)
#       h <- stats::hclust(d, input$dendroDistanceMetric)
#       if (input$clusteringAlgorithmD == "Phylogenetic Tree") {
#         vals$dendrogram <- ggtree::ggtree(as.phylo(h), layout = "circular", open.angle = 360) + ggtree::geom_tiplab2(size = 2)
#       } else if (input$clusteringAlgorithmD == "Hierarchical") {
#         vals$dendrogram <- ggtree::ggtree(as.phylo(h)) + ggtree::theme_tree2() + ggtree::geom_tiplab(size = 2)
#       } else {
#         stop("Input clustering algorithm not found ", input$clusteringAlgorithmD)
#       }
#       vals$dendrogram
#     }
#   })
# })
# output$treePlot <- renderPlot({
#   req(vals$dendrogram)
#   vals$dendrogram
# }, height = 600)
# 
