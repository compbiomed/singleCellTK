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
        p("1. Choose algorithm (one of the PCA algorithm)"),
        p("2. Choose components to be involved in paired plot"),
        p("3. Choose 2 components in Checkbox to produce single plot(if more than 2 components are selected, only the first 2 will be used)"),
        p("3. Choose feature to color data by"),
        p("4. Visualize your data")
      )
    ),#subtab "Principle Component Analysis" end
    
    tabPanel("Dimensionality Reduction", 
             fluidRow(
               column(4,
                      wellPanel(
                        selectInput("selectDimRed","Algorithm",c("PCA","tSNE")),
                        selectInput("pcX", "X axis:", pcComponents),
                        selectInput("pcY", "Y axis:", pcComponents, selected = "PC2"),
                        selectInput("colorClusters","Color Clusters By",clusterChoice),
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
    ),#tSNE end
    tabPanel("Clustering",
             fluidRow(
               column(4,
                      wellPanel("Select clustering algorithm:",
                                selectInput("selectCluster", "Algorithm", c("K-Means")),
                                selectInput("selectK", "Cluster Centers", c(1, 2, 3, 4, 5)),
                                p("Select data to cluster:"),
                                selectInput("selectDataC", "Data", c("Raw Data", "PCA Components", "tSNE Components")),
                                selectInput("pcX", "X axis:", c("1"="PC1","2"="PC2","3"="PC3")),
                                selectInput("pcY", "Y axis:", c("1"="PC1","2"="PC2","3"="PC3")),
                                selectInput("colorClusters","Color Clusters By",c("Tissue")),
                                actionButton("plotClusters", "Plot Clusters")
                      )
               ),
               column(8,
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
