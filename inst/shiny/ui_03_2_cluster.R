shinyPanelCluster <- fluidPage(
  tags$script("Shiny.addCustomMessageHandler('close_dropDownClust', function(x){
                  $('html').click();
                });"),
  tags$div(
    class = "container",
    h1("Clustering"),
    h5(tags$a(href = paste0(docs.artPath, "clustering.html"),
              "(help)", target = "_blank")),
    wellPanel(
      # CLUSTERING --> VISUALIZATION
      fluidRow(
        column(
          6,
          selectInput("clustAlgo", "Select Algorithm",
                      list("Scran SNN" = c("louvain" = 1, "leiden" = 2, 
                                           "walktrap" = 3, "infomap" = 4, 
                                           "fast greedy" = 5, "label prop" = 6, 
                                           "leading eigen" = 7),
                           "K-Means" = c("Hartigan-Wong" = 8, "Lloyd" = 9,
                                         "MacQueen" = 10),
                           "Seurat" = c("louvain" = 11, "multilevel" = 12,
                                        "SLM" = 13)),
                      )
        )
      ),
      h4("Input Parameters:"),
      fluidRow(
        # Scran SNN ####
        conditionalPanel(
          "input.clustAlgo >=1 && input.clustAlgo <= 7",
          column(
            6,
            selectizeInput(
              inputId = "clustScranSNNMat", 
              label = "Select input matrix:", 
              choices = NULL, 
              selected = NULL, 
              multiple = FALSE,
              options = NULL),
            #uiOutput('clustScranSNNMat'),
          ),
          column(
            4,
            numericInput("clustScranSNNK", "K value:", 8, min = 1, step = 1),
          ),
          column(
            4,
            numericInput("clustScranSNNd", "Number of Components:",
                         50, min = 2, step = 5)
          ),
          column(
            4,
            selectInput("clustScranSNNType", "Edge Weight Type:",
                        c("rank", "number", "jaccard"), selected = "rank")
          ),
          conditionalPanel(
            condition = 'input.clustAlgo == 2',
            column(
              4,
              numericInput("clustScranSNNLeidenReso", "Resolution:",
                           value = 1, min = 0, step = 0.1)
            ),
            column(
              6,
              selectInput("clusterScranSNNLeidenObjFunc", "Objective Function:",
                          c("Constant Potts Model (CPM)" = "CPM",
                            "Modularity" = "modularity"), selected = "CPM")
            )
          ),
          conditionalPanel(
            condition = 'input.clustAlgo == 3',
            column(
              4,
              numericInput("clustScranSNNWalktrapStep", "Steps:",
                           value = 4, min = 1, step = 1)
            )
          )
        ),

        # K-Means ####
        conditionalPanel(
          "input.clustAlgo >= 8 && input.clustAlgo <= 10",
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
            numericInput("clustKMeansN", "Number of Centers (Clusters):",
                         value = NULL),
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
        ),

        # Seurat ####
        conditionalPanel(
          "input.clustAlgo >= 11 && input.clustAlgo <= 13",
          column(
            6,
            selectInput("clustSeuratReddim", "Select A ReducedDim:", currreddim)
          ),
          column(6),
          column(
            12,
            helpText("A 'reducedDim' contains low-dimension representation of an assay.\n Dimension reduction has to be run in advance.")
          ),
          column(
            4,
            numericInput("clustSeuratDims", "How Many Dimensions to Use:", 10,
                         min = 2, step = 1),
          ),
          column(
            4,
            checkboxInput("clustSeuratGrpSgltn", "Group Singletons",
                          value = TRUE)
          ),
          column(
            4,
            numericInput("clustSeuratRes", "Resolution", 0.8, step = 0.05)
          )
        )
      ), # fuildRow ends here
      useShinyjs(),
      textInput("clustName", "Name of Clustering Result:",
                ""),
      #uiOutput("clustNameUI"),
      withBusyIndicatorUI(actionButton("clustRun", "Run"))
    ),
    h3("Visualization"),
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
  ),
  nonLinearWorkflowUI(id = "nlw-cl")
)

