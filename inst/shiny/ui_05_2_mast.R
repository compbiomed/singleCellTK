shiny_panel_mast <- fluidPage(
  tags$div(
    class = "container",
    h1("MAST"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v07-tab05_Differential-Expression.html#mast",
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        selectInput("mastAssay", "Select Assay:", currassays),
        h3("Adaptive Thresholding:"),
        withBusyIndicatorUI(actionButton("runThreshPlot",
                                         "Run Thresholding")),
        h3("Hurdle Model:"),
        checkboxInput("useAdaptThresh", "Use Adaptive Thresholds",
                      value = TRUE),
        sliderInput("FCthreshold", "Select fold change threshold", log2(1),
                    log2(2), 0.6),
        sliderInput("hurdlethresh", "Select expression threshold", 0, 1,
                    0.1),
        selectInput("hurdlecondition", "Select Condition for Hurdle Model",
                    clusterChoice),
        uiOutput("hurdleconditionofinterestUI"),
        sliderInput("hurdlepvalue", "p-value (FDR) cutoff:", 0.01, 0.2, 0.05),
        withBusyIndicatorUI(actionButton("runDEhurdle",
                                         "Run DE Using Hurdle")),
        h5("Download Results"),
        downloadButton("downloadHurdleResult", "Download Results")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Adaptive thresholding", plotOutput("threshplot")),
          tabPanel("Results Table", DT::dataTableOutput("mastresults")),
          tabPanel("Violin Plot", plotOutput("hurdleviolin")),
          tabPanel("Linear Model", plotOutput("hurdlelm")),
          tabPanel("Heatmap", plotOutput("hurdleHeatmap"))
        )
      )
    )
  )
)
