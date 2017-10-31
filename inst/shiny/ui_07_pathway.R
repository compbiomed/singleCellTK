shiny_panel_pathway <- fluidPage(
  tags$div(
    class = "container",
    h1("Pathway Activity Analysis"),
    sidebarLayout(
      sidebarPanel(
        selectInput("pathwayAssay", "Select Assay:", currassays),
        selectInput("pathwayMethod", "Select Method:", c("GSVA")),
        radioButtons("genelistSource", "Gene list source:",
                     c("Manual Input", "MSigDB c2 (Human, Entrez ID only)")),
        uiOutput("selectPathwayGeneLists"),
        uiOutput("selectNumTopPaths"),
        selectInput("pathwayPlotVar",
                    "Select Condition(s) of interest for plot:", clusterChoice,
                    multiple = TRUE),
        radioButtons("pathwayOutPlot", "Plot Type:", c("Violin", "Heatmap")),
        withBusyIndicatorUI(actionButton("pathwayRun", "Run")),
        tags$hr(),
        downloadButton("downloadPathway", "Download Pathway Results")
      ),
      mainPanel(
        plotOutput("pathwayPlot", height = "600px")
      )
    )
  )
)
