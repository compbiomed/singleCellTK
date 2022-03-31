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
                      list("Scran SNN" = c("walktrap" = 1, "louvain" = 2,
                                           "infomap" = 3, "fastGreedy" = 4,
                                           "labelProp" = 5, "leadingEigen" = 6),
                           "K-Means" = c("Hartigan-Wong" = 7, "Lloyd" = 8,
                                         "MacQueen" = 9),
                           "Seurat" = c("louvain" = 10, "multilevel" = 11,
                                        "SLM" = 12)),
                      )
        )
      ),
      h4("Input Parameters:"),
      fluidRow(
        # Scran SNN ####
        conditionalPanel(
          "input.clustAlgo >=1 && input.clustAlgo <= 6",
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
            numericInput("clustScranSNNK", "K value:", 10, min = 1, step = 1),
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
          "input.clustAlgo >= 7 && input.clustAlgo <= 9",
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
          "input.clustAlgo >= 10 && input.clustAlgo <= 12",
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
      uiOutput("clustNameUI"),
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

