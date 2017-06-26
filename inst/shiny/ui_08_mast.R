shiny_panel_mast <- fluidPage(
  tags$div(
    class="container",
    h1("MAST"),
    fluidPage(
      fluidRow(
        column(
          4,
          wellPanel(
            h3("Adaptive Thresholding:"),
            withBusyIndicatorUI(actionButton("runThreshPlot", "Run Thresholding")),
            h3("Hurdle Model:"),
            checkboxInput("useAdaptThresh", "Use Adaptive Thresholds", value=TRUE),
            sliderInput("FCthreshold", "Select fold change threshold", log2(1), log2(2), 0.6),
            sliderInput("hurdlethresh", "Select expression threshold", 0, 1, 0.1),
            sliderInput("samplesize", "Select Sample Size", 0, 500, 100),
            selectInput("hurdlecondition","Select Condition for Hurdle Model",clusterChoice),
            uiOutput("hurdleconditionofinterestUI"),
            sliderInput("hurdlepvalue", "p-value cutoff:", 0.01, 0.2, 0.05),
            selectInput("selectCorrection","Correction Method:",c("FDR", "Bonferroni")),
            withBusyIndicatorUI(actionButton("runDEhurdle", "Run DE Using Hurdle")),
            h5("Download Results"),
            downloadButton("downloadHurdlerResult","Download Results")
          )
        ),
        column(
          8,
          tabsetPanel(
            tabPanel("Adaptive thresholding", plotOutput("threshplot")),
            tabPanel("Results Table", dataTableOutput("mastresults")),
            tabPanel("Violin Plot", plotOutput("hurdleviolin")),
            tabPanel("Linear Model", plotOutput("hurdlelm")),
            tabPanel("Heatmap", plotOutput("hurdleHeatmap"))
          )
        )
      )
    )
  )
)
