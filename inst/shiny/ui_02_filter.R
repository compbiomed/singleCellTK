shiny_panel_filter <- fluidPage(
  tags$div(
    class = "container",
    h1("Data Summary"),
    fluidPage(
      fluidRow(
        column(8,
          tableOutput("summarycontents"),
          plotlyOutput("countshist"),
          plotlyOutput("geneshist"),
          dataTableOutput("contents")),
        column(
          4,
          wellPanel(
            checkboxInput("removeNoexpress", "Remove genes with 0 expression across all samples (Recommended)", value = TRUE),
            numericInput("minDetectGenect", label = "Minimum Detected Genes per Sample.", value = 1700, min = 1, max = 100000),
            numericInput("LowExpression", "% Low Gene Expression to Filter", value = 40, min = 0, max = 100),
            h2("Delete Outliers"),
            selectInput("deletesamplelist", "Select Samples:",
                        sampleChoice,
                        multiple = TRUE),
            actionButton("filterData", "Filter Data"),
            actionButton("resetData", "Reset"),
            tags$hr(),
            h3("Filter samples by annotation"),
            uiOutput("FilterSamples"),
            uiOutput("filterSampleOptions"),
            tags$hr(),
            h3("Filter genes by feature annotation"),
            uiOutput("filterFeatures"),
            uiOutput("filterFeatureOptions"),
            tags$hr(),
            h3("Delete an annotation column:"),
            selectInput("deletepdatacolumn", "Annotation Column:", clusterChoice),
            actionButton("deletepDatabutton", "Delete Column"),
            tags$hr(),
            downloadButton("downloadSCESet", "Download SCEset")
          )
        )
      )
    )
  ),
  includeHTML("www/footer.html")
)
