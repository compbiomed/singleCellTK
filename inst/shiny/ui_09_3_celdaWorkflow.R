# User Interface for Celda Workflow ---
shinyPanelCelda <- fluidPage(
    h1("Celda"),
    h5(tags$a(href = "https://www.sctk.science/articles/tab09_celda-workflow",
              "(help)", target = "_blank")),
    inlineCSS(list(".panel-danger>.panel-heading" = "background-color:#dcdcdc; color:#000000", ".panel-primary>.panel-heading" = "background-color:#f5f5f5; color:#000000; border-color:#dddddd", ".panel-primary" = "border-color:#dddddd;", ".panel-primary>.panel-heading+.panel-collapse>.panel-body" = "border-color:#dddddd;")),
    bsCollapse(id = "CeldaUI", open = "Data Input",
        bsCollapsePanel("Identify Number of Feature Modules",
           fluidRow(
               column(4,
                   panel(
                       selectInput("celdaassayselect", "Choose an Assay:",
                                   choices = c()),
                       selectInput("celdafeatureselect", "Choose Feature Selection Method:",
                                   choices = c("Simple Filter", "SeuratFindHVG", "Scran_modelGeneVar")),
                       conditionalPanel("input.celdafeatureselect == 'SeuratFindHVG'",
                                        selectInput("celdaseurathvgmethod", "Select HVG method:",
                                                    choices = c("vst", "dispersion", "mean.var.plot"))
                       ),
                       conditionalPanel("input.celdafeatureselect == 'Simple Filter'",
                                        numericInput("celdarowcountsmin",
                                                     "Keep features with this many counts:", value = 3),
                                        numericInput("celdacolcountsmin",
                                                     "In at least this many cells:", value = 3)
                       ),
                       conditionalPanel("input.celdafeatureselect != 'Simple Filter'",
                                        numericInput("celdafeaturenum",
                                                     "Select number of variable features:", min = 1, max = 5000, value = 2000)
                       ),
                       numericInput("celdaLinit", "Select Number of Initial Feature Modules:", min = 1, max = 25, value = 10),
                       numericInput("celdaLmax", "Select Number of Maximum Feature Modules:", min = 15, max = 200, value = 100),
                       actionButton("celdamodsplit", "Recursive Module Split"),
                       hidden(
                           numericInput("celdaLselect", "Select Number of Feature Modules:", min = 1, max = 100, value = 25),
                           actionButton("celdaLbtn", "Select Number of Modules")
                       )
                   )
               ),
               column(8,
                   fluidRow(
                       column(12,
                           hidden(
                               tags$div(class = "celda_modsplit_plots", tabsetPanel(id = "celdaModsplitTabset", type = "tabs"
                               ))
                           )
                       )
                   )
               )
           ),
            style = "primary"),

        bsCollapsePanel("Identify Number of Cell Clusters",
            fluidRow(
                column(4,
                    panel(
                        numericInput("celdaKinit", "Select Number of Initial Cell Modules:", min = 1, max = 10, value = 5),
                        numericInput("celdaKmax", "Select Number of Maximum Cell Modules:", min = 15, max = 40, value = 25),
                        actionButton("celdacellsplit", "Recursive Cell Split"),
                        hidden(
                            numericInput("celdaKselect", "Select Number of Cell Clusters:", min = 2, max = 20, value = 10),
                            actionButton("celdaKbtn", "Select Number of Clusters")
                        )
                    )
                ),
                column(8,
                    fluidRow(
                        column(12,
                            hidden(
                                tags$div(class = "celda_cellsplit_plots",
                                    tabsetPanel(
                                    tabPanel("Perplexity Plot",
                                             panel(
                                                 plotlyOutput(outputId = "plot_cellsplit_perp", height = "auto")
                                             )
                                    ),
                                    tabPanel("Perplexity Diff Plot",
                                             panel(
                                                 plotlyOutput(outputId = "plot_cellsplit_perpdiff", height = "auto")
                                             )
                                    ),
                                    tabPanel("Preliminary UMAP Plots",
                                            uiOutput("celdaKplots")
                                    )
                                )
                                )
                            )
                        )
                    )
                )
            ),
            style = "primary"),

        bsCollapsePanel("Visualization",
            tabsetPanel(
                tabPanel("UMAP",
                    fluidRow(
                        column(4,
                            panel(
                                numericInput("celdaUMAPmaxCells",
                                             label =
                                                 "Max.cells: Maximum number of cells to plot",
                                             value = 25000,
                                             min = 1,
                                             step = 1),
                                numericInput("celdaUMAPminClusterSize",
                                             label =
                                                 "Min.cluster.size: Do not subsample cell
                              clusters below this threshold",
                                             value = 100,
                                             min = 1,
                                             step = 1),
                                numericInput("celdaUMAPSeed",
                                             label =
                                                 "Seed: ",
                                             value = 12345),
                                numericInput("celdaUMAPmindist",
                                             label =
                                                 "Min.dist: Effective minimum distance
                                between embedded points",
                                             value = 0.75),
                                numericInput("celdaUMAPspread",
                                             label =
                                                 "Spread: ",
                                             value = 1),
                                numericInput("celdaUMAPnn",
                                             label =
                                                 "nNeighbors: ",
                                             value = 30),
                                actionButton("CeldaUmap", "Run UMAP")
                            )
                        ),
                        column(8,
                            panel(
                                plotlyOutput("celdaumapplot", height = "auto")
                            )
                        )
                    )
                ),
                tabPanel("TSNE",
                    fluidRow(
                        column(4,
                            panel(
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
                                actionButton("CeldaTsne", "Run tSNE")
                            )
                        ),
                        column(8,
                            panel(
                                plotlyOutput("celdatsneplot", height = "auto")
                            )
                        )
                    )
                ),
                tabPanel("Heatmap",
                         fluidRow(
                             panel(
                                 plotOutput("celdaheatmapplt")
                             )
                         )
                ),
                tabPanel("Module Heatmap",
                         fluidRow(
                             column(4,
                                    panel(
                                        numericInput(inputId = "celdamodheatmapnum",
                                                     label = "Select module to display on heatmap:", value = 10, step = 1),
                                        actionButton("celdamodheatmapbtn", "Plot Module Heatmap")
                                    )
                             ),
                             column(8,
                                    panel(
                                        plotOutput(outputId = "celdamodheatmapplt") %>% withSpinner(size = 3, color = "#0dc5c1", type = 8)
                                    )
                             )
                         )
                ),
                tabPanel("Probablity Map",
                        fluidRow(
                            panel(
                                plotOutput("celdaprobmapplt")
                            )
                        )
                )
            ),
            style = "primary")
    )
)
