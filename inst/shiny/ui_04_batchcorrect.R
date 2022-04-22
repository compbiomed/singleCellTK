shinyPanelBatchcorrect <- fluidPage(
  includeCSS('styles.CSS'),
  tags$script("Shiny.addCustomMessageHandler('close_dropDownBC', function(x){
                  $('html').click();
                });"),
  h1("Normalization & Batch Correction"),
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
                            h5(tags$a(href = paste0(docs.artPath, "cnsl_normalization.html"),
                                      "(help)", target = "_blank")),
                            selectInput(
                              inputId = "normalizeAssayMethodSelect",
                              label = "Select normalization method: ",
                              choices = c("Seurat - LogNormalize" = "LogNormalize",
                                          "Seurat - CLR" = "CLR",
                                          "Seurat - RC" = "RC",
                                          "Seurat - SCTransform" = "SCTransform",
                                          "Scater - LogNormCounts" = "logNormCounts",
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
                                       selectizeInput(
                                         inputId = "normalizeAssaySelect", 
                                         label = "Select input matrix:", 
                                         choices = NULL, 
                                         selected = NULL, 
                                         multiple = FALSE,
                                         options = NULL),
                                       #uiOutput("normalizeAssaySelect"),
                                       conditionalPanel(
                                         condition = "input.normalizeAssayMethodSelect == 'LogNormalize'
                              || input.normalizeAssayMethodSelect == 'CLR'
                              || input.normalizeAssayMethodSelect == 'RC'",
                                         numericInput(
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
                                       ),
                                       conditionalPanel(
                                         condition = "input.normalizationScale == true",
                                         awesomeCheckbox(
                                           inputId = "normalizationTrim",
                                           label = "Trim data after scale?",
                                           value = TRUE
                                         ),
                                         conditionalPanel(
                                           condition = "input.normalizationTrim == true",
                                           numericInput(
                                             inputId = "normalizationTrimUpper",
                                             label = "Upper trim value:",
                                             value = 10
                                           ),
                                           numericInput(
                                             inputId = "normalizationTrimLower",
                                             label = "Lower trim value:",
                                             value = -10
                                           )
                                         )
                                       )
                                     ),
                                     conditionalPanel(
                                       condition = "input.normalizeAssayMethodSelect == 'custom'",
                                       h5("Assay Options:"),
                                       selectizeInput(
                                         inputId = "modifyAssaySelect", 
                                         label = "Select input matrix:", 
                                         choices = NULL, 
                                         selected = NULL, 
                                         multiple = FALSE,
                                         options = NULL),
                                       #uiOutput("modifyAssaySelect"),
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
                                                       "Seurat - SCTransform" = "SCTransform",
                                                       "Scater - LogNormCounts" = "logNormCounts",
                                                       "Scater - CPM" = "CPM")
                                         )
                                       ),
                                       conditionalPanel(
                                         condition = "!(input.customNormalizeOptionsNormalize == true
                                         && (input.customNormalizeAssayMethodSelect == 'LogNormalize'
                                         || input.customNormalizeAssayMethodSelect == 'CLR'
                                         || input.customNormalizeAssayMethodSelect == 'SCTransform'
                                         || input.customNormalizeAssayMethodSelect == 'logNormCounts'))",
                                         awesomeCheckbox(
                                           inputId = "customNormalizeOptionsTransform",
                                           label = "Transform",
                                           value = FALSE
                                         )
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
                                         conditionalPanel(
                                           condition = "input.customNormalizeOptionsNormalize == true",
                                           awesomeCheckbox(
                                             inputId = "customNormalizePseudoOptionsBefore",
                                             label = "before normalization",
                                             value = FALSE
                                           ),
                                           conditionalPanel(
                                             condition = "input.customNormalizePseudoOptionsBefore == true",
                                             numericInput(
                                               inputId = "customNormalizePseudoValueBefore",
                                               label = "Enter a pseudovalue to add:",
                                               value = 1,
                                               min = 1
                                             )
                                           )
                                         ),
                                         conditionalPanel(
                                           condition = "input.customNormalizeOptionsTransform == true",
                                           awesomeCheckbox(
                                             inputId = "customNormalizePseudoOptionsAfter",
                                             label =  "before transformation",
                                             value = FALSE
                                           ),
                                           conditionalPanel(
                                             condition = "input.customNormalizePseudoOptionsAfter == true",
                                             numericInput(
                                               inputId = "customNormalizePseudoValueAfter",
                                               label = "Enter a pseudovalue to add:",
                                               value = 1,
                                               min = 1
                                             )
                                           )
                                         ),
                                         conditionalPanel(
                                           condition = "input.customNormalizeOptionsNormalize == false
                                           && input.customNormalizeOptionsTransform == false",
                                           h6("To add a pseudovalue, must select 'Normalize' or 'Transform' options!")
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
                                     ),
                                     tags$hr(),
                                     h4("Output Data Type:"),
                                     uiOutput("normalizationDataTagUI")
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
        h5(tags$a(href = paste0(docs.artPath, "batch_correction.html"),
                  "(help)", target = "_blank")),
        selectInput('batchCorrMethods', "Select Batch Correction Method:",
                    c("ComBatSeq", "BBKNN", "FastMNN", "Limma", #"Harmony", "LIGER",
                      "MNN", "scanorama", "scMerge", "Seurat3 Integration",
                      "ZINBWaVE")),
        selectizeInput(
          inputId = "batchCorrAssay", 
          label = "Select input matrix:", 
          choices = NULL, 
          selected = NULL, 
          multiple = FALSE,
          options = NULL),
        #uiOutput("batchCorrAssay"),
        selectInput("batchCorrVar", "Select Batch Annotation:", clusterChoice),
        # BBKNN ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'BBKNN'",
          h5(tags$a(href = paste0(docs.base, "reference/runBBKNN.html"),
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
          h5(tags$a(href = paste0(docs.base, "reference/runComBatSeq.html"),
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
          h5(tags$a(href = paste0(docs.base, "reference/runFastMNN.html"),
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
        #   h5(tags$a(href = paste0(docs.base, "reference/runHarmony.html"),
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
        #   h5(tags$a(href = paste0(docs.base, "reference/runLIGER.html"),
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
          h5(tags$a(href = paste0(docs.base, "reference/runLimmaBC.html"),
                    "(help for Limma)", target = "_blank")),
          textInput("limmaSaveAssay", "Assay Name to Use:", value = "LIMMA"),
          withBusyIndicatorUI(actionButton("limmaRun", "Run"))
        ),
        # MNN ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'MNN'",
          h5(tags$a(href = paste0(docs.base, "reference/runMNNCorrect.html"),
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
          h5(tags$a(href = paste0(docs.base, "reference/runSCANORAMA.html"),
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
          h5(tags$a(href = paste0(docs.base, "reference/runSCMerge.html"),
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
          h5(tags$a(href = paste0(docs.base, "reference/runSeuratIntegration.html"),
                    "(help for Seurat Integration)", target = "_blank")),
          uiOutput('Srt3IntNAnchUI'),
          textInput("Srt3IntSaveAssay", "Assay Name to Use:",
                    value = "Seurat3Int"),
          withBusyIndicatorUI(actionButton("Srt3IntRun", "Run"))
        ),
        # ZINBWaVE ####
        conditionalPanel(
          condition = "input.batchCorrMethods == 'ZINBWaVE'",
          h5(tags$a(href = paste0(docs.base, "reference/runZINBWaVE.html"),
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
        fluidRow(
          class = "batchVarPlotRow",
          panel(
            heading = "Visualization",

            fluidRow(
              column(
                width = 4,
                dropdown(
                  fluidRow(
                    column(
                      width = 12,
                      fluidRow(actionBttn(inputId = "closeDropDownBC", label = NULL, style = "simple", color = "danger", icon = icon("times"), size = "xs"), align = "right"),
                      h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/batch_correction.html#visualization",
                                "(What are plotted?)", target = "_blank")),
                      selectizeInput(
                        inputId = "batchCheckOrigAssay", 
                        label = "Select input matrix:", 
                        choices = NULL, 
                        selected = NULL, 
                        multiple = FALSE,
                        options = NULL),
                      #uiOutput("batchCheckOrigAssay"),
                      #selectInput("batchCheckOrigAssay", "Original Assay:", currassays),
                      selectInput("batchCheckVar", "Batch Annotation:", clusterChoice),
                      selectInput("batchCheckCond", "Additional Condition (optional)",
                                  clusterChoice),
                      p("Only result generated in the current session will be presented. ",
                        style = "color:grey;"),
                      uiOutput("batchCheckResUI"),
                      withBusyIndicatorUI(
                        actionBttn(
                        inputId = "plotBatchCheck",
                        label = "Update",
                        style = "bordered",
                        color = "primary",
                        size = "sm"
                      )
                      )
                    )
                  ),
                  inputId = "dropDownBC",
                  icon = icon("cog"),
                  status = "primary",
                  circle = FALSE,
                  inline = TRUE
                )
              ),
              column(
                width = 7,
                fluidRow(
                  h6(
                    "The top two plots shows the variance explained by the grouping of batches and user specified conditions, and the bottom two plots present the low dimension representation of the datasets. Please refer to our documentation for detatil."
                  ),
                  align="center"
                )
              )
            ),
            hr(),
            br(),


            column(
              width = 6,
              style='border-right: 1px solid #CCCCCC',
              h4('Original Status'),
              plotOutput('batchOriVars',
                height = "300px", width = "300px"),
              plotOutput('batchOriPCA',
                height = "300px", width = "400px")
            ),
            column(
              width = 6,
              h4("Corrected Status"),
              plotOutput('batchCorrVars',
                height = "300px", width = "300px"),
              plotOutput('batchCorrReddim',
                height = "300px", width = "400px")
            ),
          )
        )
      )
    )
  )
  ),
  nonLinearWorkflowUI(id = "nlw-nbc")
)

