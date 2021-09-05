shinyPanelfindMarker <- fluidPage(
  tags$div(
    class = "container",
    h1("Find Marker"),
    h5(tags$a(href = paste0(docs.artPath, "ui_find_marker.html"),
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        selectInput('fmMethod', "Select Differential Expression Method",
                    c("wilcox", "MAST", "DESeq2", "Limma", "ANOVA")),
        uiOutput('fmAssay'),
        selectInput("fmCluster", "Cluster Annotation", clusterChoice),
        selectInput("fmCovar", "Covariate(s)", clusterChoice, multiple = TRUE),
        numericInput("fmLogFC", "Log2FC greater than",
                     value = 0, min = 0, step = 0.05),
        numericInput("fmFDR", "FDR less than", value = 0.05,
                     max = 1, min = 0.01, step = 0.01),
        numericInput("fmMinClustExprPerc",
                     "Minimum Expressed Percentage in Cluster",
                     value = 0, min = 0, max = 1),
        numericInput("fmMaxCtrlExprPerc",
                     "Maximum Expressed Percentage in Control",
                     value = 1, min = 0, max = 1),
        numericInput("fmMinMeanExpr",
                     "Minimum Mean Expression Level in Cluster",
                     value = 0, min = 0),
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
            shinyjs::useShinyjs(),
            wellPanel(
              actionButton(inputId = "fmShowHMSetting", "Show Settings"),
              shinyjs::hidden(
                tags$div(
                  id = "fmHMsettings",
                  fluidRow(
                    column(
                      width = 6,
                      radioButtons('fmHMOrder', "Order blocks by",
                                   c("size", "name"), inline = TRUE),
                    ),
                    column(
                      width = 6,
                      checkboxInput("fmHMdec", "Decreasing", TRUE),
                    )
                  ),
                  fluidRow(
                    column(
                      width = 6,
                      checkboxInput("fmUseTopN", "Plot Top N markers of each cluster",
                                    TRUE)
                    ),
                    column(
                      width = 6,
                      numericInput("fmTopN", NULL, 5, min = 1, step = 1)
                    )
                  ),
                  fluidRow(
                    column(
                      width = 6,
                      numericInput("fmHMFC", "Plot Log2FC greater than",
                                   value = 0, min = 0, step = 0.05),
                    ),
                    column(
                      width = 6,
                      numericInput("fmHMFDR", "Plot FDR less than",
                                   value = 0.05, max = 1, step = 0.01),
                    )
                  ),
                  fluidRow(
                    column(
                      width = 4,
                      numericInput("fmHMMinClustExprPerc",
                                   "Plot Cluster Expression Percentage More Than",
                                   value = 0.5, min = 0, max = 1)
                    ),
                    column(
                      width = 4,
                      numericInput("fmHMMaxCtrlExprPerc",
                                   "Plot Control Expression Percentage Less Than",
                                   value = 0.4, min = 0, max = 1)
                    ),
                    column(
                      width = 4,
                      numericInput("fmHMMinMeanExpr",
                                   "Plot Mean Cluster Expression Level More Than",
                                   value = 0, min = 0)
                    )
                  ),
                  fluidRow(
                    column(
                      width = 6,
                      selectInput("fmHMrowData", "Additional feature annotation",
                                  featureChoice, multiple = TRUE),
                      withBusyIndicatorUI(actionButton('plotFM', "Plot"))
                    ),
                    column(
                      width = 6,
                      selectInput("fmHMcolData", "Additional cell annotation",
                                  clusterChoice, multiple = TRUE)
                    )
                  )
                )
              )
            ),
            shinyjqui::jqui_resizable(plotOutput('fmHeatmap'))
          )
        )
      )
    )
  )
)

