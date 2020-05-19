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
                                inputId = "normalizeLibrarySelect",
                                label = "Select normalization:",
                                choices = c(
                                    "Seurat" = "seurat",
                                    "CPM" = "cpm")
                                ),
                            selectInput("normalizeAssaySelect", "Select Assay:", currassays),
                            conditionalPanel(
                                condition = "input.normalizeLibrarySelect == 'seurat'",
                                selectInput(
                                    inputId = "normalizeAssayMethodSelect",
                                    label = "Select normalization method: ",
                                    choices = c("LogNormalize", "CLR", "RC")
                                    ),
                                textInput(
                                inputId = "normalizationScaleFactor",
                                label = "Set scaling factor: ",
                                value = "10000"
                                    )
                            ),
                            conditionalPanel(
                                condition = "input.normalizeLibrarySelect == 'cpm'",
                                conditionalPanel(
                                    condition = "input.assayModifyAction != 'delete'",
                                    textInput("normalizeAssayOutname", "Assay Name", "",
                                    placeholder = "What should the assay be called?")
                                    ),
                                ),
                            withBusyIndicatorUI(actionButton("normalizeAssay", "Normalize"))
                        )
                    )
                ),
                fluidRow(
                    column(
                        width = 12,
                        panel(
                            heading = "Assay Options",
                            selectInput(
                                "assayModifyAction",
                                "Assay Actions:",
                                c(
                                "Log Transform" = "log",
                                "log1p" = "log1p",
                                "Z-Score" = "z.score"
                                )
                            ),
                            selectInput("modifyAssaySelect", "Select Assay:", currassays),
                            conditionalPanel(
                                condition = "input.assayModifyAction != 'delete'
                                    && input.assayModifyAction != 'seurat.scale'",
                                textInput("modifyAssayOutname", "Assay Name", "",
                                placeholder = "What should the assay be called?")
                            ),
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
                            materialSwitch(
                                inputId = "trimAssayCheckbox",
                                label = "Trim Assay",
                                value = FALSE
                            ),
                            conditionalPanel(
                                condition = "input.trimAssayCheckbox == true",
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
  tabPanel("Batch Correction",
    fluidRow(
      panel(
        heading = "Basic input option:",
        div(style="display: inline-block;vertical-align:top; width: 150px;",
            selectInput("batchCorrAssay", "Select Assay:", currassays)),
        div(style="display: inline-block;vertical-align:top; width: 200px;",
            selectInput("batchCorrVar", "Select Batch Annotation:", clusterChoice)),
        div(style="display: inline-block;vertical-align:top; width: 150px;",
            selectInput("batchCorrCond", "Select Condition:", clusterChoice))
      )
    ),
    sidebarLayout(
      sidebarPanel(
        h3("Parameters"),
        selectInput('batchCorrMethods', "Select Batch Correction Method:",
                    c("ComBat", "FastMNN", "LIGER", "Limma", "MNN", "scMerge")),
        # ComBat ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'ComBat'",
          radioButtons("combatParametric", "Adjustments:",
                       c("Parametric", "Non-parametric"),
                       selected = "Parametric"),
          checkboxInput("combatMeanOnly",
                        "Correct mean of the batch effect only",
                        value = FALSE),
          checkboxInput("combatRef", "Run reference batch combat:",
                        value = FALSE),
          uiOutput("selectCombatRefBatchUI"),
          textInput("combatSaveAssay", "Assay Name to Use:", value = "combat"),
          withBusyIndicatorUI(actionButton("combatRun", "Run"))
        ),
        # FastMNN ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'FastMNN'",
          checkboxInput('FastMNNPcInput', "Use low-dimension input instead",
                        value = FALSE),
          conditionalPanel(
            condition = 'input.FastMNNPcInput == true',
            selectInput('FastMNNReddim', "Select Reduced dimension:", currreddim)
          ),
          textInput("FastMNNSaveReddim", "ReducedDim Name to Use:", value = "FastMNN"),
          withBusyIndicatorUI(actionButton("FastMNNRun", "Run"))
        ),
        # LIGER ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'LIGER'",
          numericInput("ligerNComp", label = "Number of output dimension:",
                       value = 20L, min = 2, max = 100000, step = 1),
          numericInput("ligerLambda", label = "Lambda:",
                       value = 5.0, min = 0, max = 100000, step = 0.1),
          numericInput("ligerResolution", label = "Resolution:",
                       value = 1.0, min = 0, max = 100000, step = 0.1),
          textInput("ligerSaveReddim", "ReducedDim Name to Save:",
                    value = "LIGER"),
          withBusyIndicatorUI(actionButton("ligerRun", "Run"))
        ),
        # Limma ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'Limma'",
          textInput("limmaSaveAssay", "Assay Name to Use:", value = "LIMMA"),
          withBusyIndicatorUI(actionButton("limmaRun", "Run"))
        ),
        # MNN ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'MNN'",
          numericInput('MNNK', 'K value',
                       value = 20, min = 1, max = 100000, step = 1),
          textInput("MNNSigma", 'Sigma value', value = 0.1,
                    placeholder = 'Type a number.'),
          textInput("MNNSaveAssay", "Assay Name to Use:", value = "MNN"),
          withBusyIndicatorUI(actionButton("MNNRun", "Run"))
        ),
        # scMerge ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'scMerge'",
          ## SEG options
          radioButtons('scMergeSEGOpt', "Choose Stable Expressed Gene (SEG) seg",
                       choiceNames = c('Automatically identify',
                                       'Use pre-computed dataset',
                                       'Customized list'),
                       choiceValues = c(1, 2, 3)),
          conditionalPanel(
            condition = 'input.scMergeSEGOpt == 2',
            radioButtons('scMergeSEGSpecies', "The species of SEG:",
                         choices = c('human', 'mouse'))
          ),
          conditionalPanel(
            condition = 'input.scMergeSEGOpt == 3',
            textAreaInput("scMergeSEGCustom", "Use customized SEG:",
                          placeholder = "paste them here.\nOne per line with no symbol separator.")
          ),
          checkboxInput("scMergeAutoKmk", "Automatically detect Kmeans-K", value = TRUE),
          conditionalPanel(
            condition = "input.scMergeAutoKmk == false",
            uiOutput('scMergeNBatch'),
            textInput('scMergeUserKmk', label = NULL)
          ),
          textInput("scMergeSaveAssay", "Assay Name to Use:", value = "scMerge"),
          withBusyIndicatorUI(actionButton("scMergeRun", "Run"))
        ),
        uiOutput("batchCorrStatus")
      ),
      mainPanel(
        h3("Variance explained by batch and condition"),
        fluidRow(class = "batchVarPlotRow",
          column(
            width = 6,
            style='border-right: 1px solid #CCCCCC',
            h4('Original Variance'),
            plotOutput('batchVarPlot',
                       height = "500px", width = "500px"),
          ),
          column(
            width = 6,
            h4("Corrected Variance"),
            plotOutput('batchVarCorrectedPlot',
                       height = "500px", width = "500px")
          )
        ), tags$head(tags$style("
           .batchVarPlotRow{height:560px;
                            width:1050px;
                            background-color: #ECF0F1;
           }"
           )
        )
      )
    )
  )
  )
)
