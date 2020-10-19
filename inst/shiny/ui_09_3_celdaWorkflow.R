# User Interface for Celda Workflow ---
shinyPanelCelda <- fluidPage(
    inlineCSS(list(".panel-danger>.panel-heading" = "background-color:#dcdcdc; color:#000000", ".panel-primary>.panel-heading" = "background-color:#f5f5f5; color:#000000; border-color:#dddddd", ".panel-primary" = "border-color:#dddddd;", ".panel-primary>.panel-heading+.panel-collapse>.panel-body" = "border-color:#dddddd;")),
    bsCollapse(id = "CeldaUI", open = "Data Input",
        bsCollapsePanel("Identify # of Feature Modules",
           fluidRow(
               column(4,
                   panel(
                       numericInput("celdaLinit", "Select Number of Initial Feature Modules:", min = 1, max = 25, value = 10),
                       numericInput("celdaLmax", "Select Number of Maximum Feature Modules:", min = 15, max = 200, value = 100),
                       actionButton("celdamodsplit", "Recursive Module Split"),
                       hidden(
                           numericInput("celdaLselect", "Select Number of Feature Modules:", min = 1, max = 100, value = 25),
                           actionButton("celdaLbtn", "Select # of Modules")
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

        bsCollapsePanel("Identify # of Cell Clusters",
            fluidRow(
                column(4,
                    panel(
                        numericInput("celdaKinit", "Select Number of Initial Cell Modules:", min = 1, max = 10, value = 5),
                        numericInput("celdaKmax", "Select Number of Maximum Cell Modules:", min = 15, max = 40, value = 25),
                        actionButton("celdacellsplit", "Recursive Cell Split"),
                        hidden(
                            numericInput("celdaKselect", "Select Number of Cell Clusters:", min = 2, max = 20, value = 10),
                            actionButton("celdaKbtn", "Select # of Clusters")
                        )
                    )
                ),
                column(8,
                    fluidRow(
                        column(12,
                            hidden(
                                tags$div(class = "celda_cellsplit_plots", tabsetPanel(id = "celdaCellsplitTabset", type = "tabs"
                                ))
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
                                actionButton("CeldaUmap", "Run UMAP")
                            )
                        ),
                        column(8,
                            panel(
                                plotlyOutput("celdaumapplot", height = 300)
                            )
                        )
                    )
                ),
                tabPanel("TSNE",
                    fluidRow(
                        column(4,
                            panel(
                                actionButton("CeldaTsne", "Run tSNE")
                            )
                        ),
                        column(8,
                            panel(
                                plotlyOutput("celdatsneplot", height = 300)
                            )
                        )
                    )
                ),
                tabPanel("Heatmap",
                    fluidRow(
                        column(4,
                            panel(
                                actionButton("celdaheatmapbtn", "Plot Heatmap"),
                                materialSwitch(inputId = "heatmap_module", label = "Display module heatmap?", value = TRUE),
                                conditionalPanel(
                                    condition = 'input.heatmap_module == true',
                                    numericInput(inputId = "celdamodheatmapnum",
                                        label = "Select module to display on heatmap:", value = 10, step = 1),
                                ),
                            )
                        ),
                        column(8,
                            fluidRow(
                                column(12,
                                    hidden(
                                        tags$div(class = "celda_heatmap_plots", tabsetPanel(id = "celdaHeatmapTabset", type = "tabs"
                                        ))
                                    )
                                )
                            )
                        )
                    )
                ),
                tabPanel("Probablity Map",
                    column(4,
                        panel(
                            actionButton("celdaprobplotbtn", "Plot Probability Map")
                        )
                    ),
                    column(8,
                        panel(
                            plotOutput("celdaprobmapplt", height = 300)
                        )
                    )
                )
            ),
            style = "primary")
    )
)
