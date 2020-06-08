shinyPanelMASTMarker <- fluidPage(
  tags$div(
    class = "container",
    h1("MAST Find Marker"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v07-tab05_Differential-Expression.html#mast",
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        selectInput('mastFMAssay', "Select Assay", currassays),
        selectInput("mastFMCluster", "Cluster Annotation", clusterChoice),
        numericInput("mastFMLogFC", "Log2FC greater than",
                     value = 0.25, min = 0, step = 0.05),
        numericInput("mastFMFDR", "FDR less than", value = 0.05,
                     max = 1, min = 0.01, step = 0.01),
        numericInput("mastFMFreq", "Use gene expressed in more than",
                     value = 0.1, min = 0, max = 1, step = 0.05),
        checkboxInput("mastFMUseThresh", "Use Adaptive Thresholds",
                      value = FALSE),
        withBusyIndicatorUI(actionButton("runMASTFM", "Find Marker"))
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Result Table",
            uiOutput("mastFMResClusterUI"),
            DT::dataTableOutput("mastFMResTable"),
            downloadButton("mastFMDownload", "Download Results"),
          ),
          tabPanel(
            "Heatmap",
            sidebarLayout(
              position = 'right',
              sidebarPanel(
                uiOutput('mastFMHMAssayUI'),
                radioButtons('mastFMHMOrder', "Order blocks by",
                             c("size", "name")),
                checkboxInput("mastFMHMdec", "Decreasing", TRUE),
                numericInput("mastFMHMFC", "Plot Log2FC greater than",
                             value = 1, min = 0, step = 0.05),
                numericInput("mastFMHMFDR", "Plot FDR less than",
                             value = 0.05, max = 1, step = 0.01),
                selectInput("mastFMHMrowData", "Additional feature annotation",
                            featureChoice, multiple = TRUE),
                selectInput("mastFMHMcolData", "Additional cell annotation",
                            clusterChoice, multiple = TRUE)
              ),
              mainPanel(
                plotOutput('mastFMHeatmap')
              )
            )
          )
        )
      )
    )
  )
)
