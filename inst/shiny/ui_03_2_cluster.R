shinyPanelCluster <- fluidPage(
  tags$script("Shiny.addCustomMessageHandler('close_dropDownClust', function(x){
                  $('html').click();
                });"),
  tags$div(
    class = "container",
    h1("Clustering"),
    h5(tags$a(href = paste0(docs.artPath, "clustering.html"),
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel = sidebarPanel(
        # CLUSTERING --> VISUALIZATION
        selectInput("clustAlgo", "Select Algorithm",
                    list("Scran SNN" = c("louvain" = 1, "leiden" = 2, 
                                         "walktrap" = 3, "infomap" = 4, 
                                         "fast greedy" = 5, "label prop" = 6, 
                                         "leading eigen" = 7),
                         "K-Means" = c("Hartigan-Wong" = 8, "Lloyd" = 9,
                                       "MacQueen" = 10),
                         "Seurat" = c("louvain" = 11, "multilevel" = 12,
                                      "SLM" = 13),
                         "Scanpy" = c("louvain" = 14,
                                      "leiden" = 15)),
        ),
        # Scran SNN ####
        conditionalPanel(
          condition = "input.clustAlgo >=1 && input.clustAlgo <= 7",
          selectizeInput(
            inputId = "clustScranSNNMat", 
            label = "Select input matrix:", 
            choices = NULL, 
            selected = NULL, 
            multiple = FALSE,
            options = NULL),
          numericInput("clustScranSNNK", "K value:", 8, min = 1, step = 1),
          numericInput("clustScranSNNd", "Number of Components:",
                       10, min = 1),
          selectInput("clustScranSNNType", "Edge Weight Type:",
                      c("rank", "number", "jaccard"), selected = "rank"),
          conditionalPanel(
            condition = 'input.clustAlgo == 2',
            numericInput("clustScranSNNLeidenReso", "Resolution:",
                         value = 1, min = 0, step = 0.1),
            selectInput("clusterScranSNNLeidenObjFunc", "Objective Function:",
                        c("Constant Potts Model (CPM)" = "CPM",
                          "Modularity" = "modularity"), selected = "CPM")
          ),
          conditionalPanel(
            condition = 'input.clustAlgo == 3',
            numericInput("clustScranSNNWalktrapStep", "Steps:",
                         value = 4, min = 1, step = 1)
          )
        ),
        # K-Means ####
        conditionalPanel(
          "(input.clustAlgo >= 8 && input.clustAlgo <= 10) || input.clustAlgo == 14 || input.clustAlgo == 15",
          selectInput("clustKMeansReddim", "Select A ReducedDim:", currreddim),
          conditionalPanel("(input.clustAlgo >= 8 && input.clustAlgo <= 10)",
                           numericInput("clustKMeansN", "Number of Centers (Clusters):",
                                        value = NULL),
                           numericInput("clustKMeansNIter", "Max Number of Iterations:",
                                        10, min = 2, step = 1),
                           numericInput("clustKMeansNStart", "Number of Random Sets:",
                                        1, min = 1, step = 1)
                           )
        ),
        # Seurat ####
        conditionalPanel(
          "input.clustAlgo >= 11 && input.clustAlgo <= 13",
          numericInput("clustSeuratDims", "How Many Dimensions to Use:", 10,
                       min = 2, step = 1),
          checkboxInput("clustSeuratGrpSgltn", "Group Singletons",
                        value = TRUE),
          numericInput("clustSeuratRes", "Resolution", 0.8, step = 0.05)
        ),
        # Scanpy ###
        conditionalPanel(
          "input.clustAlgo == 14 || input.clustAlgo == 15",
          numericInput("clustScanpyDims", "No. of dimensions to use:", 10,
                       min = 2, step = 1),
          numericInput("clustScanpyNeighbors", "No. of neigbors:", 15,
                       min = 2, step = 1),
          numericInput("clustScanpyRes", "Resolution:", 0.8, step = 0.05),
          numericInput("clustScanpyIter", "No. of Iterations:", -1, step = 1),
          checkboxInput("clustScanpyWeights", "Use weights from knn graph:",
                        value = FALSE),
          selectInput("clustScanpyCorrMethod", "Correlation method:", c('pearson', 'kendall', 'spearman'))
        ),
        useShinyjs(),
        textInput("clustName", "Name of Clustering Result:",
                  ""),
        withBusyIndicatorUI(actionButton("clustRun", "Run"))
      ),
      mainPanel = mainPanel(
        panel(
          
          fluidRow(
            column(
              width = 3,
              dropdown(
                fluidRow(
                  column(
                    width = 12,
                    fluidRow(actionBttn(inputId = "closeDropDownClust", label = NULL, style = "simple", color = "danger", icon = icon("times"), size = "xs"), align = "right"),
                    radioButtons("clustVisChoicesType", NULL,
                                 c("Select from Current Results:" = 1,
                                   "Select from All Present Annotation:" = 2),
                                 selected = 1, inline = TRUE, ),
                    conditionalPanel(
                      "input.clustVisChoicesType == 1",
                      selectInput("clustVisRes", NULL, "")
                    ),
                    conditionalPanel(
                      "input.clustVisChoicesType == 2",
                      selectInput("clustVisCol", NULL, clusterChoice)
                    ),
                    selectInput("clustVisReddim", "Use Reduction:", currreddim),
                    withBusyIndicatorUI(
                      actionBttn(
                        inputId = "clustPlot",
                        label = "Update",
                        style = "bordered",
                        color = "primary",
                        size = "sm"
                      )
                    )
                  )
                ),
                inputId = "dropDownClust",
                icon = icon("cog"),
                status = "primary",
                circle = FALSE,
                inline = TRUE
              )
            ),
            column(
              width = 9,
              fluidRow(
                column(
                  width = 12,
                  h6(
                    "A scatter plot of the selected low-dimensionality representation of the dataset will be generated, with the selected cluster labeling colored on each dot (cell)."
                  ),
                  h6(
                    "If you are using an expression matrix or subset for calculation, please click on the cog icon on the left to specify the dimensions to plot."
                  )
                ),
                align="center"
              )
            )
          ),
          hr(),
          br(),
          shinyjqui::jqui_resizable(
            plotlyOutput("clustVisPlot")
          )
        )
      )
    )
  ),
  nonLinearWorkflowUI(id = "nlw-cl")
)

