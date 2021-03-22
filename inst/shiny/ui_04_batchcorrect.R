shinyPanelBatchcorrect <- fluidPage(
  includeCSS('styles.CSS'),
  h1("Normalization & Batch Correction"),
  h5(tags$a(href = "https://www.sctk.science/articles/tab04_batch-correction",
            "(help)", target = "_blank")),
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
                                          "Scater - CPM" = "CPM",
                                          "Custom Normalization" = "custom")
                            )
                            
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
                            heading = "Options",
                            fluidRow(
                              column(6, style='border-right:1px solid; border-color:#EDEDED;',
                                     conditionalPanel(
                                       condition = "input.normalizeAssayMethodSelect != 'custom'",
                                       uiOutput("normalizeAssaySelect"),
                                       #selectInput("normalizeAssaySelect", "Select Assay:", currassays),
                                       #uiOutput("about"),
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
                                       awesomeCheckbox(
                                         inputId = "normalizationScale",
                                         label = "Scale data after normalization?",
                                         value = FALSE
                                       )
                                     ),
                                     conditionalPanel(
                                       condition = "input.normalizeAssayMethodSelect == 'custom'",
                                       h5("Assay Options:"),
                                       uiOutput("modifyAssaySelect"),
                                       textInput("modifyAssayOutname", "Assay Name",
                                                 value = "customNormalizedAssay"),
                                       tags$hr(),
                                       h5("Select Options:"),
                                       awesomeCheckbox(
                                         inputId = "customNormalizeOptionsNormalize",
                                         label = "Normalize",
                                         value = FALSE
                                       ),
                                       conditionalPanel(
                                         condition = "input.customNormalizeOptionsNormalize == true",
                                         h5("Normalize Options:"),
                                         selectInput(
                                           inputId = "customNormalizeAssayMethodSelect",
                                           label = "Select normalization method: ",
                                           choices = c("Seurat - LogNormalize" = "LogNormalize",
                                                       "Seurat - CLR" = "CLR",
                                                       "Seurat - RC" = "RC",
                                                       "Seurat - SCTransform" = "SCT",
                                                       "Scater - LogNormCounts" = "LNC",
                                                       "Scater - CPM" = "CPM")
                                         )
                                       ),
                                       awesomeCheckbox(
                                         inputId = "customNormalizeOptionsTransform",
                                         label = "Transform",
                                         value = FALSE
                                       ),
                                       conditionalPanel(
                                         condition = "input.customNormalizeOptionsTransform == true",
                                         h5("Transformation Options:"),
                                         selectInput(
                                           inputId = "customNormalizeTransformOptions",
                                           label = "Select transformation method:",
                                           choices = c("Log2" = "log2",
                                                       "Log1p" = "log1p",
                                                       "Sqrt" = "sqrt")
                                         )
                                       ),
                                       awesomeCheckbox(
                                         inputId = "customNormalizeOptionsPsuedocounts",
                                         label = "Psuedocounts",
                                         value = FALSE
                                       ),
                                       conditionalPanel(
                                         condition = "input.customNormalizeOptionsPsuedocounts == true",
                                         h5("Pseudocounts Options:"),
                                         awesomeRadio(
                                           inputId = "customNormalizePseudoOptions",
                                           label = "Select when to add a pseudo value:",
                                           choices = c("before normalization",
                                                       "before transformation",
                                                       "after transformation")
                                         ),
                                         numericInput(
                                           inputId = "customNormalizePseudoValue",
                                           label = "Enter a pseudo value to add:",
                                           value = 1,
                                           min = 1
                                         )
                                       ),
                                       awesomeCheckbox(
                                         inputId = "customNormalizeOptionsScale",
                                         label = "Scale",
                                         value = FALSE
                                       ),
                                       conditionalPanel(
                                         condition = "input.customNormalizeOptionsScale == true",
                                         h5("Scale Options:"),
                                         selectInput(
                                           inputId = "customNormalizeScaleOptions",
                                           label = "Select scaling method:",
                                           choices = c("Z.Score" = "zscore")
                                         )
                                       ),
                                       awesomeCheckbox(
                                         inputId = "customNormalizeOptionsTrim",
                                         label = "Trim",
                                         value = FALSE
                                       ),
                                       conditionalPanel(
                                         condition = "input.customNormalizeOptionsTrim == true",
                                         h5("Trim Options:"),
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
                                       )
                                     )
                                     ),
                              column(6,
                                     h4("Description:"),
                                     textOutput("normalizeTabDescription"),
                                     tags$hr(),
                                     h4("Output Data Type:"),
                                     uiOutput("normalizationDataTagUI"),
                                     tags$hr(),
                                     h4("Selected Options:"),
                                     conditionalPanel(
                                       condition = "input.normalizeAssayMethodSelect != 'custom'",
                                       uiOutput("normalizationNormalizeSelectedMethodUI")
                                     ),
                                     conditionalPanel(
                                       condition = "(input.customNormalizeOptionsNormalize == true
                                       && input.normalizeAssayMethodSelect == 'custom')",
                                       h5("Normalize")
                                     ),
                                     conditionalPanel(
                                       condition = "input.customNormalizeOptionsTransform == true
                                       && input.normalizeAssayMethodSelect == 'custom'",
                                       h5("Transform")
                                     ),
                                     conditionalPanel(
                                       condition = "input.customNormalizeOptionsScale == true
                                       && input.normalizeAssayMethodSelect == 'custom'",
                                       h5("Scale")
                                     ),
                                     conditionalPanel(
                                       condition = "input.customNormalizeOptionsPsuedocounts == true
                                       && input.normalizeAssayMethodSelect == 'custom'",
                                       h5("Pseudocounts")
                                     ),
                                     conditionalPanel(
                                       condition = "input.customNormalizeOptionsTrim == true
                                       && input.normalizeAssayMethodSelect == 'custom'",
                                       h5("Trim")
                                     )
                                     )
                            ),
                            fluidRow(
                              tags$hr(),
                              column(12,
                                     conditionalPanel(
                                       condition = "input.normalizeAssayMethodSelect != 'custom'",
                                       div(style = "display:inline-block; float:right", withBusyIndicatorUI(actionButton("normalizeAssay", "Run"))),
                                      ),
                                     conditionalPanel(
                                       condition = "input.normalizeAssayMethodSelect == 'custom'",
                                       div(style = "display:inline-block; float:right", withBusyIndicatorUI(actionButton("modifyAssay", "Run")))                                     
                                       )
                                     )
                            )
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
        h3("Parameters"),
        selectInput('batchCorrMethods', "Select Batch Correction Method:",
                    c("ComBatSeq", "BBKNN", "FastMNN", "Limma", #"Harmony", "LIGER",
                      "MNN", "scanorama", "scMerge", "Seurat3 Integration",
                      "ZINBWaVE")),
        uiOutput("batchCorrAssay"),
        selectInput("batchCorrVar", "Select Batch Annotation:", clusterChoice),
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
        # ComBatSeq ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'ComBatSeq'",
          h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/references/runComBatSeq.html",
                    "(help for ComBatSeq)", target = "_blank")),
          radioButtons("combatKnownCT", "Have known cell type variable?",
                       choices = c("Yes", "No")),
          conditionalPanel(
            condition = "input.combatKnownCT == 'Yes'",
            selectInput("combatCond", "Select Condition of Covariance:",
                        clusterChoice),
          ),
          conditionalPanel(
            condition = "input.combatKnownCT == 'No'",
            radioButtons("combatCTBalance", "Are cell types balanced?",
                         choices = c("Yes", "No")),
            conditionalPanel(
              condition = "input.combatCTBalance == 'No'",
              p("Will estimate surrogate variables and use them as an empirical control"),
            )
          ),
          selectInput("combatBioCond", "Select Biology Condition:",
                      clusterChoice),
          checkboxInput("combatShrink",
                        "Apply shrinkage on parameter estimation",
                        value = FALSE),
          conditionalPanel(
            condition = "input.combatShrink == true",
            numericInput("combatNGene",
                         "Number of random genes to use in empirical Bayes estimation",
                         value = NULL)
          ),
          checkboxInput("combatShrinkDisp",
                        "Apply shrinkage on dispersion",
                        value = FALSE),
          textInput("combatSaveAssay", "Assay Name to Use:", value = "ComBatSeq"),
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
        # # Harmony ####
        # conditionalPanel(
        #   condition = "input.batchCorrMethods == 'Harmony'",
        #   h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/references/runHarmony.html",
        #             "(help for Harmony)", target = "_blank")),
        #   checkboxInput('HarmonyPcInput', "Use low-dimension input instead",
        #                 value = FALSE),
        #   conditionalPanel(
        #     condition = 'input.HarmonyPcInput == true',
        #     selectInput('HarmonyReddim', "Select Reduced dimension:",
        #                 currreddim)
        #   ),
        #   numericInput("HarmonyNComp", label = "Number of output dimension:",
        #                value = 50L, min = 2, max = 100000, step = 1),
        #   textInput("HarmonyTheta", "Theta value", value = '5',
        #             placeholder = "Type a number"),
        #   numericInput("HarmonyNIter", "Number of iteration",
        #                value = 10L, min = 1, step = 1),
        #   textInput("HarmonySaveReddim", "ReducedDim Name to Use:",
        #             value = "Harmony"),
        #   withBusyIndicatorUI(actionButton("HarmonyRun", "Run"))
        # ),
        # # LIGER ####
        # conditionalPanel(
        #   condition = "input.batchCorrMethods == 'LIGER'",
        #   h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/references/runLIGER.html",
        #             "(help for LIGER)", target = "_blank")),
        #   numericInput("ligerNComp", label = "Number of output dimension:",
        #                value = 20L, min = 2, max = 100000, step = 1),
        #   numericInput("ligerLambda", label = "Lambda:",
        #                value = 5.0, min = 0, max = 100000, step = 0.1),
        #   numericInput("ligerResolution", label = "Resolution:",
        #                value = 1.0, min = 0, max = 100000, step = 0.1),
        #   textInput("ligerSaveReddim", "ReducedDim Name to Save:",
        #             value = "LIGER"),
        #   withBusyIndicatorUI(actionButton("ligerRun", "Run"))
        # ),
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
              uiOutput("batchCheckOrigAssay"),
              #selectInput("batchCheckOrigAssay", "Original Assay:", currassays),
              selectInput("batchCheckVar", "Batch Annotation:", clusterChoice),
              selectInput("batchCheckCond", "Additional Condition (optional)",
                clusterChoice),
              p("Only result generated in the current session will be presented. ",
                style = "color:grey;"),
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

