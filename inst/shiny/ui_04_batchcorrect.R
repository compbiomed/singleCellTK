shinyPanelBatchcorrect <- fluidPage(
  includeCSS('styles.CSS'),
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
                              label = "Select normalization method: ",
                              choices = c("Seurat - LogNormalize" = "LogNormalize",
                                          "Seurat - CLR" = "CLR",
                                          "Seurat - RC" = "RC",
                                          "Seurat - SCTransform" = "SCT",
                                          "Scater - LogNormCounts" = "LNC",
                                          "Scater - CPM" = "CPM")
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
                            textInput(
                              inputId = "normalizeAssayOutname",
                              label = "Assay Name:",
                              value = "SeuratLogNormalize"
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
                            textInput("modifyAssayOutname", "Assay Name",
                                      value = "countsLog"),
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
  tabPanel(
    "Batch Correction",
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/batch_correction.html#ui-usage-1",
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        h3("Parameters"),
        #uiOutput("batchCorrAssayUI"),
        selectInput("batchCorrAssay", "Select Assay:", currassays),
        selectInput("batchCorrVar", "Select Batch Annotation:", clusterChoice),
        selectInput('batchCorrMethods', "Select Batch Correction Method:",
                    c("ComBat", "BBKNN", "FastMNN", "Harmony", "LIGER", "Limma",
                      "MNN", "scanorama", "scMerge", "Seurat3 Integration",
                      "ZINBWaVE")),
        # BBKNN ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'BBKNN'",
          h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/references/runBBKNN.html",
                    "(help for BBKNN)", target = "_blank")),
          numericInput("BBKNNNComp", label = "Number of output dimension:",
                       value = 50L, min = 2, max = 100000, step = 1),
          textInput("BBKNNSaveReddim", "ReducedDim Name to Use:",
                    value = "BBKNN"),
          withBusyIndicatorUI(actionButton("BBKNNRun", "Run"))
        ),
        # ComBat ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'ComBat'",
          h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/references/runComBat.html",
                    "(help for ComBat)", target = "_blank")),
          selectInput("combatCond", "Select Condition of Covariance:",
                      clusterChoice),
          radioButtons("combatParametric", "Adjustments:",
                       c("Parametric", "Non-parametric"),
                       selected = "Parametric"),
          checkboxInput("combatMeanOnly",
                        "Correct mean of the batch effect only",
                        value = FALSE),
          checkboxInput("combatRef", "Run reference batch combat:",
                        value = FALSE),
          uiOutput("selectCombatRefBatchUI"),
          textInput("combatSaveAssay", "Assay Name to Use:", value = "ComBat"),
          withBusyIndicatorUI(actionButton("combatRun", "Run"))
        ),
        # FastMNN ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'FastMNN'",
          h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/references/runFastMNN.html",
                    "(help for fastMNN)", target = "_blank")),
          checkboxInput('FastMNNPcInput', "Use low-dimension input instead",
                        value = FALSE),
          conditionalPanel(
            condition = 'input.FastMNNPcInput == true',
            selectInput('FastMNNReddim', "Select Reduced dimension:",
                        currreddim)
          ),
          textInput("FastMNNSaveReddim", "ReducedDim Name to Use:",
                    value = "FastMNN"),
          withBusyIndicatorUI(actionButton("FastMNNRun", "Run"))
        ),
        # Harmony ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'Harmony'",
          h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/references/runHarmony.html",
                    "(help for Harmony)", target = "_blank")),
          checkboxInput('HarmonyPcInput', "Use low-dimension input instead",
                        value = FALSE),
          conditionalPanel(
            condition = 'input.HarmonyPcInput == true',
            selectInput('HarmonyReddim', "Select Reduced dimension:",
                        currreddim)
          ),
          numericInput("HarmonyNComp", label = "Number of output dimension:",
                       value = 50L, min = 2, max = 100000, step = 1),
          textInput("HarmonyTheta", "Theta value", value = '5',
                    placeholder = "Type a number"),
          numericInput("HarmonyNIter", "Number of iteration",
                       value = 10L, min = 1, step = 1),
          textInput("HarmonySaveReddim", "ReducedDim Name to Use:",
                    value = "Harmony"),
          withBusyIndicatorUI(actionButton("HarmonyRun", "Run"))
        ),
        # LIGER ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'LIGER'",
          h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/references/runLIGER.html",
                    "(help for LIGER)", target = "_blank")),
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
          h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/references/runLimmaBC.html",
                    "(help for Limma)", target = "_blank")),
          textInput("limmaSaveAssay", "Assay Name to Use:", value = "LIMMA"),
          withBusyIndicatorUI(actionButton("limmaRun", "Run"))
        ),
        # MNN ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'MNN'",
          h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/references/runMNNCorrect.html",
                    "(help for MNN)", target = "_blank")),
          numericInput('MNNK', 'K value',
                       value = 20, min = 1, step = 1),
          numericInput("MNNSigma", 'Sigma value', value = 0.1),
          textInput("MNNSaveAssay", "Assay Name to Use:", value = "MNN"),
          withBusyIndicatorUI(actionButton("MNNRun", "Run"))
        ),
        # scanorama ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'scanorama'",
          h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/references/runSCANORAMA.html",
                    "(help for Scanorama)", target = "_blank")),
          numericInput('scnrmSIGMA', 'Sigma value', 15),
          numericInput('scnrmALPHA', 'Alpha value', 0.1),
          numericInput('scnrmKNN', 'KNN value', 20, min = 1, step = 1),
          textInput("scnrmSaveAssay", "Assay Name to Use:",
                    value = "SCANORAMA"),
          withBusyIndicatorUI(actionButton("scnrmRun", "Run"))
        ),
        # scMerge ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'scMerge'",
          h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/references/runSCMerge.html",
                    "(help for scMerge)", target = "_blank")),
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
          checkboxInput("scMergeAutoKmk", "Automatically detect Kmeans-K",
                        value = TRUE),
          conditionalPanel(
            condition = "input.scMergeAutoKmk == false",
            uiOutput('scMergeNBatch'),
            textInput('scMergeUserKmk', label = NULL)
          ),
          selectInput("scMergeCT", "Cell Type Annotation:", clusterChoice),
          textInput("scMergeSaveAssay", "Assay Name to Use:",
                    value = "scMerge"),
          withBusyIndicatorUI(actionButton("scMergeRun", "Run"))
        ),
        # Seurat3 Integration ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'Seurat3 Integration'",
          h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/references/runSeurat3Integration.html",
                    "(help for Seurat Integration)", target = "_blank")),
          uiOutput('Srt3IntNAnchUI'),
          textInput("Srt3IntSaveAssay", "Assay Name to Use:",
                    value = "Seurat3Int"),
          withBusyIndicatorUI(actionButton("Srt3IntRun", "Run"))
        ),
        # ZINBWaVE ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'ZINBWaVE'",
          h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/references/runZINVWaVE.html",
                    "(help for ZINBWaVE)", target = "_blank")),
          span("Test on small example not passed yet, don't run.",
               style = 'color:red;'),
          uiOutput('zinbwaveNHvgUI'),
          uiOutput('zinbwaveEpsUI'),
          numericInput('zinbwaveNIter', 'Max number of iteration:',
                       value = 10, step = 1),
          numericInput('zinbwaveNComp', 'Number of dimensions to output:',
                       value = 50, min = 2, step = 1),
          textInput("zinbwaveSaveReddim", "ReducedDim Name to Use:",
                    value = "ZINBWaVE"),
          withBusyIndicatorUI(actionButton("zinbwaveRun", "Run"))
        ),
        uiOutput("batchCorrStatus")
      ),
      mainPanel(
        fluidRow(class = "batchVarPlotRow",
          panel(heading = "Visualization",
            column(
              width = 4,
              style='border-right: 1px solid #CCCCCC',
              h4('Original Status'),
              plotOutput('batchOriVars',
                height = "300px", width = "300px"),
              plotOutput('batchOriPCA',
                height = "300px", width = "300px")
            ),
            column(
              width = 4,
              h4("Corrected Status"),
              plotOutput('batchCorrVars',
                height = "300px", width = "300px"),
              plotOutput('batchCorrReddim',
                height = "300px", width = "300px")
            ),
            column(
              width = 4,
              h3("Visualization Setting"),
              h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/batch_correction.html#visualization",
                        "(What are plotted?)", target = "_blank")),
              selectInput("batchCheckOrigAssay", "Original Assay:", currassays),
              selectInput("batchCheckVar", "Batch Annotation:", clusterChoice),
              selectInput("batchCheckCond", "Additional Condition (optional)",
                clusterChoice),
              uiOutput("batchCheckResUI"),
              withBusyIndicatorUI(actionButton("plotBatchCheck", "Plot"))
            )
          )
        )
      )
    )
  )
  )
)

