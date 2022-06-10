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
          heading = "1. Compute variability metric",
          h5(tags$a(
            href = paste0(docs.artPath, "cnsl_feature_selection.html"),
            "(help)",
            target = "_blank"
          )),
          selectInput(
            inputId = "hvgMethodFS",
            label = "Select method for calculating variability of features:",
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
          actionButton("findHvgButtonFS","Run")
        ),
        panel(
          heading = "2. Select number of variable features",
          disabled(selectInput("hvgMetricSelect", "Select variability metric:",
                               choices = NULL)),
          disabled(numericInput(
            "hvgNumberSelect", "Select the number of highly variable features:",
            value = 2000, step = 100, min = 0)),
          disabled(textInput(
            inputId = "hvgSubsetName",
            label = "Name of the feature subset:",
            value = NULL
          )),
          disabled(actionButton("hvgSubsetRun", "Run"))
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
                       fluidRow(actionBttn(inputId = "closeDropDownFS", 
                                           label = NULL, style = "simple", 
                                           color = "danger", 
                                           icon = icon("times"), size = "xs"), 
                                align = "right"),
                       disabled(
                         selectInput(
                           inputId = "hvgPlotMethod",
                           label = "Select variability metric (created in step 1): ",
                           choices = NULL
                         )
                       ),
                       disabled(
                         selectInput(
                           inputId = "hvgPlotSubsetSelect",
                           label = "Select feature subset (created in step 2): ",
                           choices = "None"
                         )
                       ),
                       disabled(
                         numericInput(
                           inputId = "hvgPlotNLabel",
                           label = "Number of features to label:",
                           value = 20, min = 0, step = 10
                         )
                       ),
                       selectInput(
                         inputId = "hvgPlotFeatureDisplay",
                         label = "Choose feature label:",
                         choices = c("Rownames (Default)",
                                     featureChoice)
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
            column(7, fluidRow(h6("Scatterplot showing the variability of each feature versus its average expression across all cells"), align="center"))
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
                          selectInput(
                            "dimRedPlotMethod",
                            "Select method:",
                            c(
                              "Scater - PCA" = "scaterPCA",
                              "Seurat - PCA" = "seuratPCA",
                              "Seurat - ICA" = "seuratICA"
                            )
                          ),
                          selectInput(
                            inputId = "dimRedHVGSelect",
                            label = "Select HVG list:",
                            choices = "None"
                          ),
                          checkboxInput(
                            inputId = "dimRedScale",
                            label = "Scale",
                            value = TRUE
                          ),
                          uiOutput("dimRedNameUI"),
                          numericInput(
                            inputId = "dimRedNumberDims",
                            label = "Number of dimensions:",
                            value = 50,
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
                          numericInput(inputId = "seed_dimRed",
                                       label = "Seed value for reproducibility of result:",
                                       value = 12345,
                                       step = 1),
                          actionButton("runDimred", "Run")
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
                          selectInput("dimRedPlotMethod_tsneUmap", "Select method:",
                                      c("scaterUMAP" = "scaterUMAP",
                                        "rtSNE" = "rTSNE",
                                        "seuratUMAP" = "seuratUMAP",
                                        "seuratTSNE" = "seuratTSNE")),
                          uiOutput("dimRedNameUI_tsneUmap"),
                          conditionalPanel(
                            condition = "input.dimRedPlotMethod_tsneUmap == 'scaterUMAP' ||
                                         input.dimRedPlotMethod_tsneUmap == 'rTSNE'",
                            checkboxInput("logNorm_tsneUmap", 
                                          " Log Normalize the data",
                                          FALSE),
                            selectInput("hvg_tsneUmap", "Use HVG list", 
                                        choices = "None"),
                            checkboxInput("scale_tsneUmap", "Scale assay data", 
                                          TRUE),
                            checkboxInput("pca_tsneUmap", 
                                          "Run PCA on assay data", 
                                          TRUE),
                          ),
                          numericInput(
                            inputId = "dimRedNumberDims_tsneUmap",
                            label = "Number of dimensions to use:",
                            value = 10),
                          conditionalPanel(
                            condition = "input.dimRedPlotMethod_tsneUmap == 'scaterUMAP'",
                            numericInput("iterUMAP", "# of iterations", 
                                         min = 50, max = 500, value = 200),
                            numericInput("neighborsUMAP", "# of nearest neighbors",
                                         min = 2, max = 100, value = 30),
                            numericInput("mindistUMAP",
                                         "minimum distance between points",
                                         min = 0.001, max = 0.1, value = 0.01),
                            numericInput("alphaUMAP", "learning rate(alpha)",
                                         value = 1),
                            numericInput("spreadUMAP", "spread", min = 0.001, 
                                         value = 1)
                          ),
                          conditionalPanel(
                            condition = "input.dimRedPlotMethod_tsneUmap == 'rTSNE'",
                            numericInput("iterTSNE", "No. of iterations:", 
                                         min = 100, max = 2000, value = 1000),
                            numericInput("thetaTSNE", "Set theta value:", 0.5, 
                                         min = 0, step = 0.1),
                            numericInput("perplexityTSNE", "Set perplexity:",
                                         min = 5, max = 50, value = 30)
                          ),
                          conditionalPanel(
                            condition = "input.dimRedPlotMethod_tsneUmap == 'seuratTSNE'
                                           || input.dimRedPlotMethod_tsneUmap == 'seuratUMAP'",
                            selectInput(inputId = "reductionMethodUMAPTSNEDimRed",
                                        label = "Select reduction method: ",
                                        choices = c("pca", "ica")),
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
                          numericInput(inputId = "seed__tsneUmap",
                                       label = "Seed value for reproducibility of result:",
                                       value = 12345,
                                       step = 1),
                          actionButton("runDimred_tsneUmap", "Run")
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
