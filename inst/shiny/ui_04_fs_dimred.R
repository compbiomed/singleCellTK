shinyPanelFS_DimRed <- fluidPage(
  tags$script("Shiny.addCustomMessageHandler('close_dropDownDimRedEmbedding', function(x){
                  $('html').click();
                });"),
  tags$script("Shiny.addCustomMessageHandler('close_dropDownFS', function(x){
                  $('html').click();
                });"),
  h1("Feature Selection & Dimensionality Reduction"),
  tabsetPanel(
    id = "FSDimRedTabsetPanel",
    tabPanel(
      "Feature Selection",
      fluidRow(column(
        4,
        panel(
          heading = "Compute HVG",
          h5(tags$a(
            href = paste0(docs.artPath, "cnsl_feature_selection.html"),
            "(help)",
            target = "_blank"
          )),
          selectInput(
            inputId = "hvgMethodFS",
            label = "Select HVG method: ",
            choices = c(
              "Seurat - vst" = "vst",
              "Seurat - mean.var.plot" = "mean.var.plot",
              "Seurat - dispersion" = "dispersion",
              "Scran - modelGeneVar" = "modelGeneVar"
            )
          ),
          selectizeInput(
            inputId = "assaySelectFS_Norm", 
            label = "Select input matrix:", 
            choices = NULL, 
            selected = NULL, 
            multiple = FALSE,
            options = NULL),
          #uiOutput("assaySelectFS_Norm"),
          withBusyIndicatorUI(actionButton("findHvgButtonFS",
                                           "Compute Variability"))
        ),
        panel(
          heading = "Select and Subset",
          fluidRow(h6("selection of features will be based on the latest computation above"), align="center"),
          numericInput("hvgNumberSelect", "Number of HVG to select",
                       2000, step = 100),
          selectizeInput(
            inputId = "hvgSubsetAssay", 
            label = "Select input matrix:", 
            choices = NULL, 
            selected = NULL, 
            multiple = FALSE,
            options = NULL),
          #uiOutput("hvgSubsetAssay"),
          textInput("hvgAltExpName", "Name for the subset",
                    "featureSubset"),
          withBusyIndicatorUI(actionButton("hvgSubsetRun", "Select"))
        )
      ),
      column(
        8,
        panel(
          heading = "Plot",
          fluidRow(
            column(4, dropdown(
              fluidRow(
                column(12,
                       fluidRow(actionBttn(inputId = "closeDropDownFS", label = NULL, style = "simple", color = "danger", icon = icon("times"), size = "xs"), align = "right"),
                       numericInput(
                         inputId = "hvgNoFeaturesViewFS",
                         label = "Select number of features to display: ",
                         value = 100
                       ),
                       actionBttn(
                         inputId = "updatePlotFS",
                         label = "Update",
                         style = "bordered",
                         color = "primary",
                         size = "sm"
                       )
                )
              ),
              inputId = "dropDownFS",
              icon = icon("cog"),
              status = "primary",
              circle = FALSE,
              inline = TRUE
            )),
            column(7, fluidRow(h6("Scatterplot of features showing the variability versus average expression"), align="center"))
          ),
          hr(),
          br(),
          shinyjqui::jqui_resizable(plotOutput(outputId = "plotFS")),
          fluidRow(h6("highlighted features"), align="center"),
          div(
            style = "margin-top: 2px;",
            verbatimTextOutput(outputId = "hvgOutputFS",
                               placeholder = TRUE)
          )
        )
      )),
      nonLinearWorkflowUI(id = "nlw-fs")
    ),
    tabPanel(
      "Dimensionality Reduction",
      h5(tags$a(
        href = paste0(docs.artPath, "cnsl_dimensionality_reduction.html"),
        "(help)",
        target = "_blank"
      )),
      fluidRow(column(4,
                      fluidRow(column(
                        12,
                        panel(
                          heading = "Options",
                          selectizeInput(
                            inputId = "dimRedAssaySelect", 
                            label = "Select input matrix:", 
                            choices = NULL, 
                            selected = NULL, 
                            multiple = FALSE,
                            options = NULL),
                          #uiOutput("dimRedAssaySelect"),
                          selectInput(
                            "dimRedPlotMethod",
                            "Select method:",
                            c(
                              "Scater - PCA" = "scaterPCA",
                              "Seurat - PCA" = "seuratPCA",
                              "Seurat - ICA" = "seuratICA"
                            )
                          ),
                          uiOutput("dimRedNameUI"),
                          numericInput(
                            inputId = "dimRedNumberDims",
                            label = "Number of dimensions:",
                            value = 10,
                            step = 1
                          ),
                          conditionalPanel(
                            condition = "input.dimRedPlotMethod != 'seuratICA'",
                            materialSwitch(
                              inputId = "computeElbowPlot",
                              label = "Compute ElbowPlot?",
                              value = TRUE
                            ),
                            materialSwitch(
                              inputId = "computeJackstrawPlot",
                              label = "Compute JackstrawPlot?",
                              value = FALSE
                            )
                          ),
                          materialSwitch(
                            inputId = "computeHeatmapPlot",
                            label = "Compute HeatmapPlot?",
                            value = TRUE
                          ),
                          conditionalPanel(
                            condition = "input.computeHeatmapPlot == true",
                            numericInput(
                              inputId = "dimRedNFeaturesHeatmap",
                              label = "Select number of features for heatmap plot:",
                              value = 20,
                              min = 2,
                              step = 1
                            )
                          ),
                          withBusyIndicatorUI(actionButton("runDimred", "Run"))
                        )
                      ))),
               column(8,
                      fluidRow(column(
                        12,
                        hidden(
                          tags$div(
                            class = "dimRedPCAICA_plotTabset_class",
                            tabsetPanel(id = "dimRedPCAICA_plotTabset",
                                        type = "tabs")
                          )
                        )
                      )))),
      nonLinearWorkflowUI(id = "nlw-dr")
    ),

    tabPanel(title = "Embedding",
             h5(
               tags$a(
                 href = paste0(docs.artPath, "cnsl_2d_embedding.html"),
                 "(help)",
                 target = "_blank"
               )
             ),
             fluidRow(
               column(4,
                      fluidRow(column(
                        12,
                        panel(
                          heading = "Options",
                          selectizeInput(
                            inputId = "dimRedAssaySelect_tsneUmap", 
                            label = "Select input matrix:", 
                            choices = NULL, 
                            selected = NULL, 
                            multiple = FALSE,
                            options = NULL),
                          #uiOutput("dimRedAssaySelect_tsneUmap"),
                          selectInput("dimRedPlotMethod_tsneUmap", "Select method:",
                                      c("rtSNE" = "rTSNE",
                                        "scaterUMAP" = "scaterUMAP",
                                        "seuratTSNE" = "seuratTSNE",
                                        "seuratUMAP" = "seuratUMAP")),
                          uiOutput("dimRedNameUI_tsneUmap"),
                          conditionalPanel(
                            condition = "input.dimRedPlotMethod_tsneUmap == 'scaterUMAP'",
                            checkboxInput("logNormUMAP", " Log Normalize the data",
                                          FALSE),
                            numericInput("iterUMAP", "# of iterations", min = 50,
                                         max = 500, value = 200),
                            numericInput("neighborsUMAP", "# of nearest neighbors",
                                         min = 2, max = 100, value = 30),
                            numericInput("mindistUMAP",
                                         "minimum distance between points",
                                         min = 0.001, max = 0.1, value = 0.5),
                            numericInput("alphaUMAP", "learning rate(alpha)",
                                         value = 1),
                            numericInput("spreadUMAP", "spread", min = 0.001, value = 5)
                          ),
                          conditionalPanel(
                            condition = "input.dimRedPlotMethod_tsneUmap == 'rTSNE'",
                            numericInput("iterTSNE", "No. of iterations:", min = 100,
                                         max = 2000, value = 1000),
                            numericInput("perplexityTSNE", "Set perplexity:",
                                         min = 5, max = 50, value = 5)
                          ),
                          conditionalPanel(
                            condition = "input.dimRedPlotMethod_tsneUmap == 'seuratTSNE'
                                           || input.dimRedPlotMethod_tsneUmap == 'seuratUMAP'",
                            selectInput(inputId = "reductionMethodUMAPTSNEDimRed",
                                        label = "Select reduction method: ",
                                        choices = c("pca", "ica")),
                            numericInput(
                              inputId = "dimRedNumberDims_tsneUmap",
                              label = "Set number of dimensions:",
                              value = 10),
                          ),
                          conditionalPanel(
                            condition = "input.dimRedPlotMethod_tsneUmap == 'seuratTSNE'",
                            numericInput(inputId = "perplexityTSNEDimRed",
                                         label = "Set perplexity:",
                                         value = 30)
                          ),
                          conditionalPanel(
                            condition = "input.dimRedPlotMethod_tsneUmap == 'seuratUMAP'",
                            numericInput(inputId = "minDistUMAPDimRed",
                                         label = "Set min.dist:",
                                         value = 0.3),
                            numericInput(inputId = "nNeighboursUMAPDimRed",
                                         label = "Set n.neighbors:",
                                         value = 30,
                                         step = 1),
                            numericInput(inputId = "spreadUMAPDimRed",
                                         label = "Set spread:",
                                         value = 1),
                          ),
                          withBusyIndicatorUI(actionButton("runDimred_tsneUmap", "Run"))
                        )
                      ))),
               column(8,
                      fluidRow(column(
                        12,
                        tags$div(class = "dimRedTSNEUMAP_plotTabset_class",
                                 tabPanel(
                                   title = "Plot",
                                   panel(heading = "Plot",
                                         fluidRow(
                                           column(4, dropdown(
                                             fluidRow(
                                               column(12,
                                                      fluidRow(actionBttn(inputId = "closeDropDownDimRedEmbedding", label = NULL, style = "simple", color = "danger", icon = icon("times"), size = "xs"), align = "right"),
                                                      selectizeInput(
                                                        inputId = "selectRedDimPlot_tsneUmap",
                                                        label = "Select reducedDim:",
                                                        choices = NULL
                                                      ),
                                                      actionBttn(
                                                        inputId = "updateRedDimPlot_tsneUmap",
                                                        label = "Update",
                                                        style = "bordered",
                                                        color = "primary",
                                                        size = "sm"
                                                      )
                                               )
                                             ),
                                             inputId = "dropDownDimRedEmbedding",
                                             icon = icon("cog"),
                                             status = "primary",
                                             circle = FALSE,
                                             inline = TRUE
                                           )),
                                           column(7, fluidRow(h6("Scatterplot of cells on a 2D embedding"), align="center"))
                                         ),
                                         hr(),
                                         br(),
                                         plotlyOutput(outputId = "plotDimRed_tsneUmap"))
                                 ))
                      )))
             ))
  )
)
