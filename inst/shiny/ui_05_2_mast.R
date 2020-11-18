shinyPanelMAST <- fluidPage(
  tags$div(
    class = "container",
    h1("MAST"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v07-tab05_Differential-Expression.html#mast",
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        tags$h4("Select Assay:"),
        selectInput("mastAssay", label = NULL, currassays),
        tags$hr(),
        tags$h4("Adaptive Thresholding:"),
        withBusyIndicatorUI(actionButton("runThreshPlot", "Run Thresholding")),
        tags$hr(),
        tags$h4("Hurdle Model:"),
        checkboxInput("useAdaptThresh", "Use Adaptive Thresholds",
                      value = TRUE),
        sliderInput("FCthreshold", "Select fold change threshold:", log2(1),
                    log2(2), 0.6),
        sliderInput("hurdlethresh", "Select expression threshold:", 0, 1, 0.1),
        selectInput("hurdlecondition", "Select Condition for Hurdle Model:",
                    clusterChoice),
        uiOutput("hurdleconditionofinterestUI"),
        sliderInput("hurdlepvalue", "P-Value (FDR) Cutoff:", 0.01, 0.2, 0.05),
        withBusyIndicatorUI(actionButton("runDEhurdle", "Run DE Using Hurdle")),
        tags$hr(),
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
