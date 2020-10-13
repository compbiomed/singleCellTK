shinyPanelCelda2 <- fluidPage(
  useShinyalert(),
  tags$div(
    class = "container",
    h1("celda: CEllular Latent Dirichlet Allocation"),
    h5(tags$a(href = "https://github.com/campbio/celda",
      "(help)", target = "_blank")),
    tabsetPanel(
      tabPanel(
        "Celda Clustering",
        wellPanel(
          br(),
          fluidRow(
            sidebarLayout(
              sidebarPanel(
                # SHINYJS ACCORDION --------------------------
                # Section 1 - Basic Settings
                actionButton("celdaBasicSet", "Basic Settings"),
                # open by default
                tags$div(id = "celdaCollapse1",
                  wellPanel(
                    selectInput("celdaAssay", "Select Assay:",
                      currassays),
                    selectInput("celdaModel",
                      "Select Celda Model:",
                      c("celda_C", "celda_G", "celda_CG"),
                      selected = "celda_CG"),
                    # c("celda_C"),
                    # selected = "celda_C"),
                    conditionalPanel(
                      condition = sprintf("input['%s'] == 'celda_C'",
                        "celdaModel"),
                      numericInput("cellClusterC",
                        label = "Number of Cell Clusters (K):",
                        value = 15,
                        min = 1,
                        max = 100,
                        step = 1)
                    ),
                    conditionalPanel(
                      condition = sprintf("input['%s'] == 'celda_G'",
                        "celdaModel"),
                      numericInput("geneModuleG",
                        label = "Number of Gene Modules (L):",
                        value = 50,
                        min = 1,
                        max = 200,
                        step = 1)
                    ),
                    conditionalPanel(
                      condition = sprintf("input['%s'] == 'celda_CG'",
                        "celdaModel"),
                      numericInput("cellClusterCG",
                        label = "Number of Cell Clusters (K):",
                        value = 15,
                        min = 1,
                        max = 100,
                        step = 1),
                      numericInput("geneModuleCG",
                        label = "Number of Gene Modules (L):",
                        value = 50,
                        min = 1,
                        max = 200,
                        step = 1))
                  )
                ),
                # Section 2 - Advanced Settings
                actionButton("celdaAdvSet", "Advanced Settings"),
                shinyjs::hidden(
                  tags$div(id = "celdaCollapse2",
                    wellPanel(
                      conditionalPanel(
                        condition = sprintf("input['%s'] == 'celda_C' ||
                  input['%s'] == 'celda_CG'", "celdaModel", "celdaModel"),
                        selectInput("celdaAlgorithm",
                          "Select Algorithm:",
                          list("Expectation Maximization" = "EM",
                            "Gibbs Sampling" = "Gibbs"),
                          selected = "Expectation Maximization")
                      ),
                      conditionalPanel(
                        condition = sprintf("input['%s'] == 'celda_C' ||
                  input['%s'] == 'celda_CG'", "celdaModel", "celdaModel"),
                        numericInput("celdaAlpha",
                          label = "Alpha",
                          value = 1,
                          min = 0.00000001,
                          max = 100000)
                      ),
                      numericInput("celdaBeta",
                        label = "Beta",
                        value = 1,
                        min = 0.00000001,
                        max = 100000),
                      conditionalPanel(
                        condition = sprintf("input['%s'] == 'celda_G' ||
                  input['%s'] == 'celda_CG'", "celdaModel", "celdaModel"),
                        numericInput("celdaDelta",
                          label = "Delta:",
                          value = 1,
                          min = 0.00000001,
                          max = 100000),
                        numericInput("celdaGamma",
                          label = "Gamma:",
                          value = 1,
                          min = 0.00000001,
                          max = 100000)
                      ),
                      numericInput("celdaMaxIter",
                        label = "Maximum Number of Iterations:",
                        value = 200,
                        min = 1,
                        max = 100000,
                        step = 1),
                      numericInput("celdaStopIter",
                        label = "Number of Converging Iterations for Gibbs
                          sampler to stop:",
                        value = 10,
                        min = 1,
                        max = 100000,
                        step = 1),
                      numericInput("celdaSplitIter",
                        label = "Split on Every This Number of Iteration:",
                        value = 10,
                        min = 1,
                        max = 100000,
                        step = 1),
                      numericInput("celdaNChains",
                        label = "Number of random cluster initializations
                          for every K/L combination:",
                        value = 3,
                        min = 1,
                        max = 100000,
                        step = 1),
                      # numericInput("celdaCores",
                      #   label = "Number of Cores used for parallel computing:",
                      #   value = 1,
                      #   min = 1,
                      #   max = 100000,
                      #   step = 1),
                      numericInput("celdaSeed",
                        label = "Base Seed For Random Number Generation:",
                        value = 12345,
                        min = 1,
                        max = 100000,
                        step = 1)
                    )
                  )
                ),
                tags$hr(),
                withBusyIndicatorUI(actionButton(inputId = "runCelda",
                  label = "Run Celda")),
                withBusyIndicatorUI(actionButton(inputId = "renderHeatmap",
                  label = "Render Heatmap")),
                tags$hr(),
                downloadButton("downloadSCECelda", "Download SingleCellExperiment object (.rds)")
              ),
              mainPanel(
                conditionalPanel(
                  condition = sprintf("input['%s'] ==
                  'celda_C' || input['%s'] ==
                    'celda_G' || input['%s'] == 'celda_CG'",
                    "celdaModel",
                    "celdaModel",
                    "celdaModel"),
                  plotOutput("celdaHeatmap", height = "600px")
                )
              )
            )
          )
        )
      ),
      tabPanel(
        "Celda Grid Search",
        wellPanel(
          br(),
          fluidRow(
            sidebarLayout(
              sidebarPanel(
                # SHINYJS ACCORDION --------------------------
                # Section 1 - Basic Settings
                actionButton("celdaBasicSetGS", "Basic Settings"),
                # open by default
                tags$div(id = "celdaCollapseGS1",
                  wellPanel(
                    selectInput("celdaAssayGS", "Select Assay:",
                      currassays),
                    selectInput("celdaModelGS",
                      "Select Celda Model:",
                      c("celda_C", "celda_G", "celda_CG"),
                      selected = "celda_CG"),
                    conditionalPanel(
                      condition = sprintf("input['%s'] == 'celda_C'",
                        "celdaModelGS"),

                      h4("Range of Cell Clusters (K):"),

                      numericInput("GSRangeKlow",
                        "Lower bound:",
                        value = 2,
                        min = 2,
                        step = 1),

                      numericInput("GSRangeKup",
                        "Upper bound:",
                        value = 4,
                        min = 2,
                        step = 1),

                      numericInput("interK",
                        label = "Cell Cluster Increment Step Size:",
                        value = 1,
                        min = 1,
                        step = 1)
                    ),
                    conditionalPanel(
                      condition = sprintf("input['%s'] == 'celda_G'",
                        "celdaModelGS"),

                      h4("Range of Gene Modules (L):"),

                      numericInput("GSRangeLlow",
                        "Lower bound:",
                        value = 2,
                        min = 2,
                        step = 1),

                      numericInput("GSRangeLup",
                        "Upper bound:",
                        value = 4,
                        min = 2,
                        step = 1),
                      numericInput("interL",
                        label = "Gene Module Search Increment Step Size:",
                        value = 1,
                        min = 1,
                        step = 1)
                    ),
                    conditionalPanel(
                      condition = sprintf("input['%s'] == 'celda_CG'",
                        "celdaModelGS"),
                      h4("Range of Cell Clusters (K):"),

                      numericInput("GSRangeKCGlow",
                        "Lower bound:",
                        value = 2,
                        min = 2,
                        step = 1),

                      numericInput("GSRangeKCGup",
                        "Upper bound:",
                        value = 4,
                        min = 2,
                        step = 1),

                      numericInput("interKCG",
                        label = "Cell Cluster Search Increment Step Size:",
                        value = 1,
                        min = 1,
                        step = 1),

                      h4("Range of Gene Modules (L):"),

                      numericInput("GSRangeLCGlow",
                        "Lower bound:",
                        value = 2,
                        min = 2,
                        step = 1),

                      numericInput("GSRangeLCGup",
                        "Upper bound:",
                        value = 4,
                        min = 2,
                        step = 1),

                      numericInput("interLCG",
                        label = "Gene module Search Increment Step Size:",
                        value = 1,
                        min = 1,
                        step = 1)
                    )
                    # selectInput("celdaGSVerbose", "Verbose:",
                    #   c(TRUE, FALSE),
                    #   selected = FALSE
                    # )
                  )
                ),
                # Section 2 - Advanced Settings
                actionButton("celdaAdvSetGS", "Advanced Settings"),
                shinyjs::hidden(
                  tags$div(id = "celdaCollapseGS2",
                    wellPanel(
                      numericInput("celdaMaxIterGS",
                        label = "Maximum Number of Iterations:",
                        value = 200,
                        min = 1,
                        max = 100000,
                        step = 1),
                      numericInput("celdaNChainsGS",
                        label = "Number of Random Cluster Initializations
                          For Every K/L Combination:",
                        value = 3,
                        min = 1,
                        max = 100000,
                        step = 1),
                      numericInput("celdaCoresGS",
                        label = "Number of Cores To Use For
                          Parallel Estimation of Chains:",
                        value = 1,
                        min = 1,
                        max = 100000,
                        step = 1),
                      numericInput("celdaSeedGS",
                        label = "Base Seed For Random Number Generation:",
                        value = 12345,
                        min = 1,
                        max = 100000,
                        step = 1)
                    )
                  )
                ),
                tags$hr(),
                withBusyIndicatorUI(actionButton(inputId = "runCeldaGS",
                  label = "Run Celda Grid Search")),
                withBusyIndicatorUI(actionButton(
                  inputId = "renderPerplexityPlot",
                  label = "Render Perplexity Plot")),
                tags$hr(),
                selectInput("celdaSelectGSList",
                  "Select Celda Grid Search List:",
                  NULL),
                selectInput("celdaSelectGSMod",
                  "Select Celda Model:",
                  NULL),
                withBusyIndicatorUI(actionButton(
                  inputId = "confirmCeldaModel",
                  label = "Confirm Selection")),
                tags$hr(),
                downloadButton("downloadAllCeldaLists",
                  "Download All Celda Lists")
              ),
              mainPanel(
                conditionalPanel(
                  condition = sprintf("input['%s'] ==
                  'celda_C' || input['%s'] ==
                    'celda_G' || input['%s'] == 'celda_CG'",
                    "celdaModel",
                    "celdaModel",
                    "celdaModel"),
                  plotOutput("celdaPerplexityPlot", height = "600px")
                )
              )
            )
          )
        )
      ),
      tabPanel(
        "Visualize",
        wellPanel(
          fluidRow(
            tabsetPanel(
              tabPanel("t-SNE",
                wellPanel(
                  fluidRow(
                    sidebarLayout(
                      sidebarPanel(
                        # SHINYJS ACCORDION --------------------------
                        # Settings
                        #actionButton("celdatSNESet", "Settings"),
                        # closed by default
                        #shinyjs::hidden(
                        #tags$div(id = "celdaCollapsetSNE",
                        #wellPanel(
                        selectInput("celdaAssaytSNE", "Select Assay:",
                          currassays),

                        numericInput("celdatSNEmaxCells",
                          label =
                            "Max.cells: Maximum number of cells to
                                plot",
                          value = 25000,
                          min = 1,
                          step = 1),

                        numericInput("celdatSNEminClusterSize",
                          label =
                            "Min.cluster.size: Do not subsample cell
                              clusters below this threshold",
                          value = 100,
                          min = 1,
                          step = 1),

                        numericInput("celdatSNEPerplexity",
                          label =
                            "Perplexity: ",
                          value = 20),

                        numericInput("celdatSNEmaxIter",
                          label =
                            "Max.iter: Maximum number of iterations in
                              tSNE generation",
                          value = 2500),

                        numericInput("celdatSNESeed",
                          label =
                            "Seed: ",
                          value = 12345),
                        #)
                        #),
                        #)
                        tags$hr(),
                        withBusyIndicatorUI(actionButton(
                          inputId = "runCeldatSNE",
                          label = "Run Celda t-SNE")),
                        withBusyIndicatorUI(actionButton(
                          inputId = "renderCeldatSNEByCellCluster",
                          label = "Render cell cluster t-SNE plot")),
                        withBusyIndicatorUI(actionButton(
                          inputId = "renderCeldatSNEModule",
                          label = "Render module probability t-SNE plot")),
                        br(),
                        selectInput("celdatSNEFeature", "Select Feature:",
                          NULL, multiple = TRUE),
                        withBusyIndicatorUI(actionButton(
                          inputId = "renderCeldatSNEFeature",
                          label = "Render Gene expression t-SNE plot"))
                      ), mainPanel(
                        plotOutput("celdatSNECellClusterPlot",
                          height = "600px"),
                        plotOutput("celdatSNEModulePlot",
                          height = "600px"),
                        plotOutput("celdatSNEFeaturePlot",
                          height = "600px")
                      )
                    )
                  )
                )
              ), tabPanel("Probability Map",
                wellPanel(
                  fluidRow(
                    sidebarLayout(
                      sidebarPanel(
                        selectInput("celdaAssayProbabilityMap",
                          "Select Assay:",
                          currassays),
                        tags$hr(),
                        withBusyIndicatorUI(actionButton(
                          inputId = "renderCeldaProbabilityMap",
                          label = "Render Celda Probability Map"))
                      ), mainPanel(
                        plotOutput("celdaProbabilityMapPlot",
                          height = "600px")
                      )
                    )
                  )
                )
              ),
              tabPanel("Module Heatmap",
                wellPanel(
                  fluidRow(
                    sidebarLayout(
                      sidebarPanel(
                        selectInput("celdaAssayModuleHeatmap",
                          "Select Assay:",
                          currassays),
                        selectInput("celdaFeatureModule",
                          "Select Feature Modules:",
                          NULL,
                          multiple = TRUE),
                        numericInput("celdaModuleTopCells",
                          "top.cells: Number of cells with the highest and
                          lowest probabilities for modules to include in
                          the heatmap",
                          value = 100,
                          min = 1,
                          step = 1),
                        # selectInput("celdaFeatureModuleNormalize",
                        #   "Normalize: Whether to normalize the columns of
                        #   'counts'",
                        #   choices = c(TRUE, FALSE),
                        #   selected = TRUE),
                        selectInput("celdaModuleFeatureNames",
                          "Show Feature Names:",
                          choices = c(TRUE, FALSE),
                          selected = TRUE),
                        tags$hr(),
                        withBusyIndicatorUI(actionButton(
                          inputId = "renderCeldaModuleHeatmap",
                          label = "Render Celda Module Heatmap"))
                      ), mainPanel(
                        plotOutput("celdaModuleHeatmapPlot",
                          height = "600px")
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
  )
)

