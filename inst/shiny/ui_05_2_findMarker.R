shinyPanelfindMarker <- fluidPage(
  tags$div(
    class = "container",
    h1("Find Marker"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v07-tab05_Differential-Expression.html#mast",
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        selectInput('fmAssay', "Select Assay", currassays),
        selectInput('fmMethod', "Select Differential Expression Method",
                    c("MAST", "DESeq2", "Limma")),
        selectInput("fmCluster", "Cluster Annotation", clusterChoice),
        numericInput("fmLogFC", "Log2FC greater than",
                     value = 0.25, min = 0, step = 0.05),
        numericInput("fmFDR", "FDR less than", value = 0.05,
                     max = 1, min = 0.01, step = 0.01),
        withBusyIndicatorUI(actionButton("runFM", "Find Marker"))
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Result Table",
            DT::dataTableOutput("fmResTable"),
            downloadButton("fmDownload", "Download Results"),
          ),
          tabPanel(
            "Heatmap",
            sidebarLayout(
              position = 'right',
              sidebarPanel(
                shinyjs::useShinyjs(),
                uiOutput('fmHMAssayUI'),
                radioButtons('fmHMOrder', "Order blocks by",
                             c("size", "name")),
                checkboxInput("fmHMdec", "Decreasing", TRUE),
                checkboxInput("fmUseTopN", "Plot Top N markers of each cluster",
                              TRUE),
                numericInput("fmTopN", NULL, 10, min = 1, step = 1),
                numericInput("fmHMFC", "Plot Log2FC greater than",
                             value = 0.25, min = 0, step = 0.05),
                numericInput("fmHMFDR", "Plot FDR less than",
                             value = 0.05, max = 1, step = 0.01),
                selectInput("fmHMrowData", "Additional feature annotation",
                            featureChoice, multiple = TRUE),
                selectInput("fmHMcolData", "Additional cell annotation",
                            clusterChoice, multiple = TRUE),
                withBusyIndicatorUI(actionButton('plotFM', "Plot"))
              ),
              mainPanel(
                plotOutput('fmHeatmap')
              )
            )
          )
        )
      )
    )
  )
)

