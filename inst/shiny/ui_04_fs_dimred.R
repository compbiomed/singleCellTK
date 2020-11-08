shinyPanelFS_DimRed <- fluidPage(
  h1("Feature Selection & Dimensionality Reduction"),
  h5(tags$a(href = "https://www.sctk.science/articles/tab04_fs-dimred",
            "(help)", target = "_blank")),
  tabsetPanel(
    tabPanel("Feature Selection",
             fluidRow(
               column(
                 4,
                 panel(
                   heading = "Compute HVG",
                   selectInput(
                     inputId = "hvgMethodFS",
                     label = "Select HVG method: ",
                     choices = c(
                       "Seurat - vst" = "vst",
                       "Seurat - mean.var.plot" = "mean.var.plot",
                       "Seurat - dispersion" = "dispersion",
                       "Scran - modelGeneVar" = "modelGeneVar")),
                   selectInput(
                     inputId = "assaySelectFS_Norm",
                     label = "Select normalized assay:",
                     choices = currassays),
                   withBusyIndicatorUI(actionButton("findHvgButtonFS",
                                                    "Compute Variability"))
                 ),
                 panel(
                   heading = "Select and Subset",
                   p("Selection will be based on the latest computation above."),
                   numericInput("hvgNumberSelect", "Number of HVG to select",
                                2000, step = 100),
                   textInput("hvgAltExpName", "Name for the subset",
                             "featureSubset"),
                   withBusyIndicatorUI(actionButton("hvgSubsetRun", "Select"))
                 )
               ),
               column(
                 8,
                 panel(
                   heading = "Plot",
                   p("Visualization will be based on the latest computation on the left."),
                   numericInput(
                     inputId = "hvgNoFeaturesViewFS",
                     label = "Select number of features to display: ",
                     value = 100),
                   withBusyIndicatorUI(actionButton("showHVG", "Show")),
                   div(
                     style = "margin-top: 15px;",
                     verbatimTextOutput(
                       outputId = "hvgOutputFS",
                       placeholder = TRUE)
                   ),
                   plotOutput(
                     outputId = "plotFS",
                     width = 400,
                     height = 400
                   )
                 )
               )
             )
    ),
    tabPanel(
      "Dimensionality Reduction",
      tabsetPanel(      
        tabPanel("PCA/ICA",
                            fluidRow(
                              column(4,
                                     fluidRow(
                                       column(12,
                                              panel(heading = "Options",
                                                    h6("Select assay type:"),
                                                    radioButtons(inputId = "dimRedAssayType", label = NULL,
                                                                 choices = c("Use full sized assay" = 1,
                                                                             "Use subset" = 2),
                                                                 selected = 1, inline = TRUE),
                                                    conditionalPanel(
                                                      condition = "input.dimRedAssayType == 1",
                                                      selectInput("dimRedAssaySelect", "Select assay:", currassays),
                                                    ),
                                                    conditionalPanel(
                                                      condition = "input.dimRedAssayType == 2",
                                                      selectInput("dimRedAltExpSelect", "Select subset:", curraltExps),
                                                      uiOutput("dimRedAltExpAssayUI")
                                                    ),
                                                    selectInput("dimRedPlotMethod", "Select method:",
                                                                c("Scran - PCA" = "PCA",
                                                                  "Seurat - PCA" = "PCASeurat",
                                                                  "Seurat - ICA" = "ICASeurat")),
                                                    uiOutput("dimRedNameUI"),
                                                    textInput(
                                                      inputId = "dimRedNumberDims", 
                                                      label = "Number of Dimensions:",
                                                      value = 50),
                                                    materialSwitch(
                                                      inputId = "computeElbowPlot", 
                                                      label = "Compute ElbowPlot?", 
                                                      value = TRUE),
                                                    materialSwitch(
                                                      inputId = "computeHeatmapPlot", 
                                                      label = "Compute HeatmapPlot?", 
                                                      value = TRUE),
                                                    materialSwitch(
                                                      inputId = "computeJackstrawPlot", 
                                                      label = "Compute JackstrawPlot?", 
                                                      value = FALSE),
                                                    withBusyIndicatorUI(actionButton("runDimred", "Run"))
                                                    )
                                              )
                                     )),
                              column(8,
                                     fluidRow(
                                       column(12,
                                              #hidden(
                                                tags$div(
                                                  class = "seurat_pca_plots", 
                                                  tabsetPanel(id = "dimRedPCAICA_plotTabset", 
                                                              type = "tabs"
                                                              )
                                                )#)
                                              )
                                     )
                                     )
                            )
                 ),
        tabPanel("tSNE/UMAP",
                 fluidRow(
                   column(4,
                          fluidRow(
                            column(12,
                                   panel(heading = "Options",
                                         h6("Select assay type:"),
                                         radioButtons(inputId = "dimRedAssayType_tsneUmap", label = NULL,
                                                      choices = c("Use full sized assay" = 1,
                                                                  "Use subset" = 2),
                                                      selected = 1, inline = TRUE),
                                         conditionalPanel(
                                           condition = "input.dimRedAssayType_tsneUmap == 1",
                                           selectInput("dimRedAssaySelect_tsneUmap", "Select assay:", currassays),
                                         ),
                                         conditionalPanel(
                                           condition = "input.dimRedAssayType_tsneUmap == 2",
                                           selectInput("dimRedAltExpSelect_tsneUmap", "Select subset:", curraltExps),
                                           uiOutput("dimRedAltExpAssayUI_tsneUmap")
                                         ),
                                         selectInput("dimRedPlotMethod_tsneUmap", "Select method:",
                                                     c("tSNE", "UMAP")),
                                         uiOutput("dimRedNameUI_tsneUmap"),
                                         textInput(
                                           inputId = "dimRedNumberDims_tsneUmap", 
                                           label = "Number of Dimensions:",
                                           value = 50),
                                         conditionalPanel(
                                           condition = "input.dimRedPlotMethod_tsneUmap == 'UMAP'",
                                           checkboxInput("logNormUMAP", " Log Normalize the data",
                                                         TRUE),
                                           sliderInput("iterUMAP", "# of iterations", min = 50,
                                                       max = 500, value = 200),
                                           sliderInput("neighborsUMAP", "# of nearest neighbors",
                                                       min = 2, max = 100, value = 5),
                                           sliderInput("mindistUMAP",
                                                       "minimum distance between points",
                                                       min = 0.001, max = 0.1, value = 0.01),
                                           numericInput("alphaUMAP", "learning rate(alpha)",
                                                        value = 1)
                                         ),
                                         conditionalPanel(
                                           condition = "input.dimRedPlotMethod_tsneUmap == 'tSNE'",
                                           sliderInput("iterTSNE", "# of iterations", min = 100,
                                                       max = 2000, value = 1000),
                                           sliderInput("perplexityTSNE", "Perplexity paramter",
                                                       min = 5, max = 50, value = 5)
                                         ),
                                         withBusyIndicatorUI(actionButton("runDimred_tsneUmap", "Run"))
                                   )
                            )
                          )),
                   column(8,
                          fluidRow(
                            column(12,
                                   #hidden(
                                   tags$div(
                                     class = "seurat_pca_plots", 
                                     tabsetPanel(id = "dimRedTSNEUMAP_plotTabset", 
                                                 type = "tabs"
                                     )
                                   )#)
                            )
                          )
                   )
                 )
                 )
        )
    )
  )
)
