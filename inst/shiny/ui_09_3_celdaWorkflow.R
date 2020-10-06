# User Interface for Celda Workflow ---
shinyPanelCelda <- fluidPage(
    inlineCSS(list(".panel-danger>.panel-heading" = "background-color:#dcdcdc; color:#000000", ".panel-primary>.panel-heading" = "background-color:#f5f5f5; color:#000000; border-color:#dddddd", ".panel-primary" = "border-color:#dddddd;", ".panel-primary>.panel-heading+.panel-collapse>.panel-body" = "border-color:#dddddd;")),
    bsCollapse(id = "CeldaUI", open = "Data Input",
        bsCollapsePanel("Module Splitting",
            sidebarPanel(
                numericInput("celdaLinit", "Select Number of Initial Feature Modules:", min = 1, max = 25, value = 10),
                numericInput("celdaLmax", "Select Number of Maximum Feature Modules:", min = 15, max = 200, value = 100),
                actionButton("celdamodsplit", "Recursive Module Split"),
                #actionButton("celdamodsplitdiff", "Plot Perplexity Difference"),
                hidden(
                    numericInput("celdaLselect", "Select Number of Feature Modules:", min = 1, max = 100, value = 25),
                    actionButton("celdaLbtn", "Select # of Modules")
                )
            ),
            mainPanel(
                plotlyOutput("modsplitplot", height = 300) %>% withSpinner(size = 3, color = "#0dc5c1", type = 8),
                plotlyOutput("modsplitplotdiff") %>% withSpinner(size = 3, color = "#0dc5c1", type = 8)
            )
        ),
        bsCollapsePanel("Cell Splitting",
            sidebarPanel(
                numericInput("celdaKinit", "Select Number of Initial Cell Modules:", min = 1, max = 10, value = 5),
                numericInput("celdaKmax", "Select Number of Maximum Cell Modules:", min = 15, max = 40, value = 25),
                actionButton("celdacellsplit", "Recursive Cell Split"),
                hidden(
                    numericInput("celdaKselect", "Select Number of Cell Clusters:", min = 2, max = 20, value = 10),
                    actionButton("celdaKbtn", "Select # of Clusters")
                )
            ),
            mainPanel(
                plotlyOutput("cellsplitplot", height = 300) %>% withSpinner(size = 3, color = "#0dc5c1", type = 8)
            )
        ),

        bsCollapsePanel("Umap/Tsne",
            sidebarPanel(
                actionButton("CeldaUmap", "Run UMAP"),
                actionButton("CeldaTsne", "Run tSNE")
            ),
            mainPanel()
        ),

        bsCollapsePanel("Cluster Data",
            fluidRow(
                    panel(
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
            style = "primary"),

        bsCollapsePanel("Factorize Data",
            fluidRow(
                column(4,
                    panel(
                        actionButton(inputId = "celda_factorize", "Factorize")
                    )
                )
            ),
            style = "primary")
    )

)