shiny_panel_cluster <- fluidPage(
  tabsetPanel(
    tabPanel(
      "PCA (MAST example)",
      fluidRow(
        column(
          4,
          wellPanel(
            selectInput(inputId ="pcaAlgorithm",label = "Algorithm",choices = c("regular PCA","randomized PCA","robust PCA","RR PCA")),
            checkboxGroupInput(inputId = "pcaCheckbox",label = "Involving variables", choices = c("PC1","PC2","PC3"),selected = c("PC1","PC2","PC3")),
            selectInput(inputId = "selectAdditionalVariables",label = "Additional variables:",choices = clusterChoice,multiple = TRUE),
            selectInput("plotTypeId","Plot Type",c("Paired Plot","Single Plot"),"Paired Plot"),
            selectInput("colorClusters_MAST","Color Clusters By",clusterChoice),
            withBusyIndicatorUI(actionButton("plotPCA", "Plot PCA Data"))
          )#wellPanel
        ),
        column(
          8,
          plotOutput("pcaPlot")
        )
      ),#fluidRow
      mainPanel(
        h1("Instructions"),
        p(""), strong("Principle Component Analysis:"), 
        p("1. Choose algorithm (one of the PCA algorithms)"),
        p("2. Choose components to be involved in paired plot"),
        p("3. Choose 2 components in Checkbox to produce single plot(if more than 2 components are selected, only the first 2 will be used)"),
        p("3. Choose feature to color data by"),
        p("4. Visualize your data")
      )
    ),#subtab "PCA (MAST example)" end
    
    tabPanel("Dimensionality Reduction", 
             fluidRow(
               column(4,
                      wellPanel(
                        selectInput("selectDimRed","Algorithm",c("PCA","tSNE")),
                        conditionalPanel(
                          condition = sprintf("input['%s'] == 'PCA'", "selectDimRed"),
                          selectInput("pcX", "X axis:", pcComponents),
                          selectInput("pcY", "Y axis:", pcComponents, selected = "PC2")
                        ),
                        conditionalPanel(
                          condition = sprintf("input['%s'] == 'tSNE'", "selectDimRed")
                        ),
                        selectInput("colorDims","Color Points By",clusterChoice),
                        withBusyIndicatorUI(actionButton("plotData", "Plot Data"))
                      )),
               column(8,
                      plotlyOutput("dimredPlot"))
             ),
             mainPanel(
               h1("Instructions"),
               p(""), strong("tSNE:"), 
               p("1. Choose algorithm (PCA or tSNE)"),
               p("2. Choose which components to use for the x and y axes"),
               p("3. Choose feature to color data by"),
               p("4. Visualize your data")
             )
    ),#Dimensionality Reduction (PCA/tSNE) end
    
    tabPanel("Clustering",
             fluidRow(
               column(4,
                      wellPanel(
                                p("Select data type:"),
                                selectInput("selectDataC", "Data", c("", "Raw Data", "PCA Components", "tSNE Components")),
                                conditionalPanel(
                                  condition = sprintf("input['%s'] == 'PCA Components'", "selectDataC"),
                                  selectInput("pcX_Clustering_Data", "X axis:", pcComponents),
                                  selectInput("pcY_Clustering_Data", "Y axis:", pcComponents, selected = "PC2")
                                ),
                            
                                p("Select clustering algorithm:"),
                                selectInput("selectCluster", "Algorithm", c("","K-Means", "Hierarchical", "Phylogenetic Tree")),
                                conditionalPanel(
                                  condition = sprintf("input['%s'] == 'K-Means'", "selectCluster"),
                                  selectInput("numberKClusters","Number of Clusters",numClusters)
                                ),
                                
                                conditionalPanel(
                                  condition = sprintf("input['%s'] == 'Hierarchical'", "selectCluster"),
                                  selectInput("numberHClusters","Number of Clusters",numClusters),
                                  withBusyIndicatorUI(actionButton("plotClustersDendogram", "Plot Dendogram"))
                                ),
                                
                                conditionalPanel(
                                  condition = sprintf("input['%s'] == 'Phylogenetic Tree'", "selectCluster"),
                                  withBusyIndicatorUI(actionButton("plotClustersPTree", "Plot Phylogenetic Tree"))
                                ),
                                
                                conditionalPanel(
                                  condition = sprintf("input['%s'] == 'K-Means'", "selectCluster"),
                                  p("Select Visualization"),
                                  selectInput("selectCOutput", "Output", c("", "PCA Plot", "Scatter Plot", "Table")),
                                  conditionalPanel(
                                    condition = sprintf("input['%s'] == 'PCA Plot'", "selectCOutput"),
                                    selectInput("pcX_Clustering_Plot", "X axis:", pcComponents),
                                    selectInput("pcY_Clustering_Plot", "Y axis:", pcComponents, selected = "PC2"),
                                    selectInput("colorClusters_Plot", "Color Clusters By", clusterChoice),
                                    selectInput("shapeClusters_Plot", "Shape Clusters By", clusterChoice),
                                    withBusyIndicatorUI(actionButton("plotClustersPCA", "Plot Clusters"))
                                  ), 
                                  conditionalPanel(
                                    condition = sprintf("input['%s'] == 'Scatter Plot'", "selectCOutput"),
                                    withBusyIndicatorUI(actionButton("plotClustersScatter", "Plot Clusters"))
                                  ),
                                  conditionalPanel(
                                    condition = sprintf("input['%s'] == 'Table'", "selectCOutput"),
                                    withBusyIndicatorUI(actionButton("makeCTable", "Generate Table"))
                                  )
                                )
                      )
               ),
               column(8,
                      conditionalPanel(
                        condition = sprintf("input['%s'] == 'Hierarchical'", "selectCluster"),
                        plotOutput("dendoPlot")
                      ),
                      conditionalPanel(
                        condition = sprintf("input['%s'] == 'Phylogenetic Tree'", "selectCluster"),
                        plotOutput("phyloPlot")
                      ),
                      plotlyOutput("clusterPlot"))
             ),
             mainPanel(
               h1("Instructions"),
               p(""), strong("Clustering:"),
               p("1. Choose clustering algorithm (K-Means, ...)"),
               p("2. Choose which data to use (raw data, principal component values, tSNE values)"), 
               p("3. Choose which components to use for the x and y axes"),
               p("4. Choose feature to color data by"),
               p("5. Visualize your data and clusters")
             )
    )
  ),#Clustering end
  includeHTML("www/footer.html")
)
