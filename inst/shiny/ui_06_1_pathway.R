shinyPanelPathway <- fluidPage(
  tags$div(
    class = "container",
    h1("Pathway Activity Analysis"),
    h5(tags$a(href = paste0(docs.artPath, "pathwayAnalysis.html"),
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        selectizeInput(
          inputId = "pathwayAssay", 
          label = "Select input matrix:", 
          choices = NULL, 
          selected = NULL, 
          multiple = FALSE,
          options = NULL),
        #uiOutput("pathwayAssay"),
        #selectInput("pathwayAssay", "Select Assay:", currassays),
        selectInput("pathwayMethod", "Select Method:", "GSVA"),
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
            shinyjqui::jqui_resizable(plotOutput("pathwayPlot"))
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

