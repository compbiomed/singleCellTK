# User Interface for Celda Workflow ---
shinyPanelCelda <- fluidPage(
    inlineCSS(list(".panel-danger>.panel-heading" = "background-color:#dcdcdc; color:#000000", ".panel-primary>.panel-heading" = "background-color:#f5f5f5; color:#000000; border-color:#dddddd", ".panel-primary" = "border-color:#dddddd;", ".panel-primary>.panel-heading+.panel-collapse>.panel-body" = "border-color:#dddddd;")),
    bsCollapse(id = "CeldaUI", open = "Data Input",
        bsCollapsePanel("Module Splitting",
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
                   panel(
                       plotlyOutput("modsplitplot", height = 300) %>% withSpinner(size = 3, color = "#0dc5c1", type = 8),
                       plotlyOutput("modsplitplotdiff") %>% withSpinner(size = 3, color = "#0dc5c1", type = 8)
                   )
               )
           ),
            style = "primary"),

        bsCollapsePanel("Cell Splitting",
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
                    panel(
                        plotlyOutput("cellsplitplot", height = 300) %>% withSpinner(size = 3, color = "#0dc5c1", type = 8)
                    )
                )
            ),
            style = "primary"),

        bsCollapsePanel("Umap/Tsne",
            fluidRow(
                column(4,
                    panel(
                        actionButton("CeldaUmap", "Run UMAP"),
                        actionButton("CeldaTsne", "Run tSNE")
                    )
                ),
                column(8,
                    panel(
                        plotlyOutput("celdaumapplot", height = 300) %>% withSpinner(size = 3, color = "#0dc5c1", type = 8),
                        plotlyOutput("celdatsneplot", height = 300) %>% withSpinner(size = 3, color = "#0dc5c1", type = 8)
                    )
                )
            ),
            style = "primary"),

        bsCollapsePanel("Plot Heatmap",
            fluidRow(
                column(4,
                    panel(
                        numericInput("celdamodheatmapselect", "Select Module to Plot:", min = 1, max = 25, value = 10),
                        actionButton("celdamodheatmapbtn", "Select Module Number")
                        )
                    ),
                column(8,
                    panel(
                        plotlyOutput("celdamodheatmapplot", height = 300) %>% withSpinner(size = 3, color = "#0dc5c1", type = 8),
                        plotlyOutput("celdaprobplot", height = 300) %>% withSpinner(size = 3, color = "#0dc5c1", type = 8),
                        )
                    )
            ),
            style = "primary")
    )

)