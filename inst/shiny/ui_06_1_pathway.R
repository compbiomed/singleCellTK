shinyPanelPathway <- fluidPage(
  tags$div(
    class = "container",
    h1("Pathway Activity Analysis"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v08-tab06_Pathway-Activity-Analysis.html",
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        selectInput("pathwayAssay", "Select Assay:", currassays),
        selectInput("pathwayMethod", "Select Method:", "GSVA"),
        radioButtons("genelistSource", "Gene list source:",
                     c("Manual Input", "MSigDB c2 (Human, Entrez ID only)")),
        uiOutput("selectPathwayGeneLists"),
        uiOutput("selectNumTopPaths"),
        selectInput("pathwayPlotVar",
                    "Select Condition(s) of interest for plot:", clusterChoice,
                    multiple = TRUE),
        radioButtons("pathwayOutPlot", "Plot Type:", c("Heatmap", "Violin")),
        withBusyIndicatorUI(actionButton("pathwayRun", "Run")),
        tags$hr(),
        h3("Save pathway results:"),
        actionButton("savePathway", "Save Pathways"),
        downloadButton("downloadPathway", "Download Pathway Results")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Plot",
            plotOutput("pathwayPlot", height = "600px")
          ),
          tabPanel(
            "Results Table",
            DT::dataTableOutput("pathwaytable")
          )
        )
      )
    )
  )
)

