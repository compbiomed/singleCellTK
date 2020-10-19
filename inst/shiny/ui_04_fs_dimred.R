shinyPanelFS_DimRed <- fluidPage(
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
      "Run New Dimensional Reduction",
      # SHINYJS COLLAPSE --------------------------
      # Section 1 - Assay Settings
      # open by default
      tags$div(
        id = "dimred",
        wellPanel(
          fluidRow(
            column(
              8,
              wellPanel(
                fluidRow(
                  column(
                    6,
                    tags$h4("Select:"),
                    radioButtons(inputId = "dimRedAssayType", label = NULL,
                                 choices = c("Use full sized assay" = 1,
                                             "Use subset" = 2),
                                 selected = 1, inline = TRUE),
                    conditionalPanel(
                      condition = "input.dimRedAssayType == 1",
                      selectInput("dimRedAssaySelect", "Assay:", currassays),
                    ),
                    conditionalPanel(
                      condition = "input.dimRedAssayType == 2",
                      selectInput("dimRedAltExpSelect", "Subset:", curraltExps),
                      uiOutput("dimRedAltExpAssayUI")
                    ),
                    # Note: Removed "Dendrogram" option from method select
                    # to disable conditionalPanels.
                    selectInput("dimRedPlotMethod", "Method:",
                                c("PCA", "tSNE", "UMAP")),
                    tags$br(),
                    withBusyIndicatorUI(actionButton("runDimred", "Run"))
                  ),
                  column(
                    6,
                    tags$h4("DR Options:"),
                    uiOutput("dimRedNameUI"),
                    #textInput("dimRedNameInput", "reducedDim Name:", ""),
                    tags$br(),
                    HTML('<button type="button"
                            class="btn btn-default btn-block"
                            data-toggle="collapse"
                            data-target="#c-collapse-run-options">
                              View More Options
                          </button>'
                    ),
                    tags$div(
                      id = "c-collapse-run-options", class = "collapse",
                      conditionalPanel(
                        condition = "input.dimRedPlotMethod == 'PCA'",
                        HTML('<p style="color:rgb(255,0,0);">
                                No parameters available for PCA
                              </p>')
                      ),
                      conditionalPanel(
                        condition = "input.dimRedPlotMethod == 'UMAP'",
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
                        condition = "input.dimRedPlotMethod == 'tSNE'",
                        sliderInput("iterTSNE", "# of iterations", min = 100,
                                    max = 2000, value = 1000),
                        sliderInput("perplexityTSNE", "Perplexity paramter",
                                    min = 5, max = 50, value = 5)
                      )
                    )
                  )
                )
              )
            ),
            column(
              4,
              wellPanel(
                h4("Available Reduced Dims:"),
                tableOutput("reducedDimsList"),
                tags$hr(),
                h4("Remove a reducedDim:"),
                fluidRow(
                  column(
                    8,
                    selectInput("delRedDimType", label = NULL, currreddim)
                  ),
                  column(
                    4,
                    withBusyIndicatorUI(actionButton("delRedDim", "Delete"))
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)

