shinyPanelBatchcorrect <- fluidPage(
  includeCSS('styles.css'),
  tabsetPanel(
  tabPanel(
    "Normalization", fluid = TRUE,
        fluidRow(
            column(
                width = 4,
                fluidRow(
                    column(
                        width = 12,
                        panel(
                            heading = "Normalization Options",
                            selectInput(
                                inputId = "normalizeAssayMethodSelect",
                                label = "Select normalization:",
                                choices = c(
                                    "Seurat - LogNormalize" = "LogNormalize",
                                    "Seurat - CLR" = "CLR",
                                    "Seurat - RC" = "RC",
                                    "Seurat - SCTransform" = "SCT",
                                    "Scater - logNormCounts" = "LNC",
                                    "CPM" = "CPM")
                                ),
                            selectInput("normalizeAssaySelect", "Select Assay:", currassays),
                            conditionalPanel(
                                condition = "input.normalizeAssayMethodSelect == 'LogNormalize'
                                              || input.normalizeAssayMethodSelect == 'CLR'
                                              || input.normalizeAssayMethodSelect == 'RC'",
                                textInput(
                                inputId = "normalizationScaleFactor",
                                label = "Set scaling factor: ",
                                value = "10000"
                                    )
                            ),
                            textInput("normalizeAssayOutname", "Assay Name", "",
                            placeholder = "What should the assay be called?"),
                            withBusyIndicatorUI(actionButton("normalizeAssay", "Normalize"))
                        )
                    )
                ),
                fluidRow(
                    column(
                        width = 12,
                        panel(
                            heading = "Transformation, Scaling, and Trimming",
                            selectInput(
                                "assayModifyAction",
                                "Options:",
                                c(
                                "Log Transform" = "log",
                                "log1p" = "log1p",
                                "Z-Score" = "z.score",
                                "Trim" = "trim"
                                )
                            ),
                            selectInput("modifyAssaySelect", "Select Assay:", currassays),
                            textInput("modifyAssayOutname", "Assay Name", "", placeholder = "What should the assay be called?"),
                            conditionalPanel(
                                condition = "input.assayModifyAction == 'seurat.scale'",
                                selectInput(
                                    inputId = "scaleSeuratModel",
                                    label = "Select model for scaling: ",
                                    choices = c("linear", "poisson", "negbinom")
                                    ),
                                materialSwitch(
                                    inputId = "scaleSeuratDoScale",
                                    label = "Scale data?",
                                    value = TRUE
                                    ),
                                materialSwitch(
                                    inputId = "scaleSeuratDoCenter",
                                    label = "Center data?",
                                    value = TRUE
                                    ),
                                textInput(
                                    inputId = "scaleSeuratMaxValue",
                                    label = "Max value for scaled data: ",
                                    value = "10"
                                    )
                            ),
                            conditionalPanel(
                              condition = "input.assayModifyAction != 'trim'",
                              materialSwitch(
                                inputId = "trimAssayCheckbox",
                                label = "Trim Assay",
                                value = FALSE
                              )
                            ),
                            conditionalPanel(
                                condition = "input.trimAssayCheckbox == true || input.assayModifyAction == 'trim'",
                                numericInput(
                                    inputId = "trimUpperValueAssay",
                                    label = "Specify upper trim value",
                                    value = 10
                                ),
                                numericInput(
                                  inputId = "trimLowerValueAssay",
                                  label = "Specify lower trim value",
                                  value = -10
                                )
                            ),
                            withBusyIndicatorUI(actionButton("modifyAssay", "Run")),
                        )
                    )
                )
            ),
            column(
                width = 8,
                fluidRow(
                    column(
                        width = 12,
                        panel(
                            heading = "Available Assays",
                            tableOutput("assayList")
                        )
                    )
                )
            )
        )
    ),
  tabPanel(
    "Batch Correction", 
    sidebarLayout(
      sidebarPanel(
        selectInput("combatAssay", "Select Assay:", currassays),
        tags$hr(),
        h4("Plot Batch Effect:"),
        selectInput("batchVarPlot", "Select Batch Annotation:", c("none", clusterChoice)),
        selectInput("conditionVarPlot", "Select Condition Annotation:", c("none", clusterChoice)),
        tags$hr(),
        h4("Run Batch Correction:"),
        selectInput("batchMethod", "Select Method:", "ComBat"),
        selectInput("combatBatchVar", "Select Batch Condition:", clusterChoice),
        selectInput("combatConditionVar", "Select Additional Covariates:",
                    clusterChoice, multiple = TRUE),
        radioButtons("combatParametric", "Adjustments:", c("Parametric",
                                                           "Non-parametric"),
                     selected = "Parametric"),
        checkboxInput("combatMeanOnly", "Correct mean of the batch effect only",
                      value = FALSE),
        checkboxInput("combatRef", "Run reference batch combat:",
                      value = FALSE),
        uiOutput("selectCombatRefBatchUI"),
        textInput("combatSaveAssay", "Assay Name to Use:", value = "combat"),
        withBusyIndicatorUI(actionButton("combatRun", "Run"))
      ),
      mainPanel(
        uiOutput("combatStatus"),
        plotOutput("combatBoxplot", height = "600px")
      )
    )
  )
))