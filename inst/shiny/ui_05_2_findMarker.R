shinyPanelfindMarker <- fluidPage(
  tags$div(
    class = "container",
    h1("Find Marker"),
    tags$a(href = "https://www.sctk.science/articles/tab05_find-marker",
           "(help)", target = "_blank"),
    sidebarLayout(
      sidebarPanel(
        selectInput('fmMethod', "Select Differential Expression Method",
                    c("wilcox", "MAST", "DESeq2", "Limma", "ANOVA")),
        uiOutput('fmAssay'),
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
                #uiOutput('fmHMAssay'),
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

