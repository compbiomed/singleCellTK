shiny_panel_filter <- fluidPage(
  useShinyalert(),
  tags$div(
    class = "container",
    h1("Data Summary and Filtering"),
    h5(tags$a(href = "https://www.sctk.science/articles/v04-tab02_Data-Summary-and-Filtering.html",
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        selectInput("filterAssaySelect", "Select Assay:", currassays),
        checkboxInput("removeNoexpress", "Remove genes with 0 expression across all samples (Recommended)", value = TRUE),
        numericInput("minDetectGene", label = "Minimum Detected Genes per Sample.", value = 1700, min = 1, max = 100000),
        numericInput("LowExpression", "% Low Gene Expression to Filter", value = 40, min = 0, max = 100),
        h2("Delete Outliers"),
        selectInput("deletesamplelist", "Select Samples:",
                    sampleChoice,
                    multiple = TRUE),
        withBusyIndicatorUI(actionButton("filterData", "Filter Data")),
        actionButton("resetData", "Reset"),
        tags$hr(),
        h3("Filter samples by annotation"),
        selectInput("filteredSample", "Select Annotation:", c("none", clusterChoice)),
        uiOutput("filterSampleOptions"),
        tags$hr(),
        h3("Filter genes by feature annotation"),
        selectInput("filteredFeature", "Select Feature:", c("none", featureChoice)),
        uiOutput("filterFeatureOptions"),
        tags$hr(),
        h3("Convert gene annotations"),
        selectInput("orgOrganism", "Select Organism:", as.character(grep("^org\\.", installed.packages()[, "Package"], value = TRUE))),
        uiOutput("orgConvertColumns"),
        withBusyIndicatorUI(actionButton("convertGenes", "Convert")),
        tags$hr(),
        h3("Delete an annotation column:"),
        selectInput("deleterowdatacolumn", "Annotation Column:", clusterChoice),
        actionButton("deleterowDatabutton", "Delete Column"),
        tags$hr(),
        h3("Randomly Subset"),
        numericInput("downsampleNum", "Number of samples to keep:", min = 2,
                     max = numSamples, value = numSamples, step = 1),
        withBusyIndicatorUI(actionButton("downsampleGo", "Subset Data")),
        tags$hr(),
        downloadButton("downloadSCE", "Download SCtkExperiment")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Data Summary",
            tableOutput("summarycontents"),
            plotlyOutput("countshist"),
            plotlyOutput("geneshist"),
            DT::dataTableOutput("contents")
          ),
          tabPanel(
            "Assay Details",
            br(),
            fluidRow(
              sidebarLayout(
                sidebarPanel(
                  h4("Assay Options:"),
                  selectInput("addAssayType", "Add Assay Type:", c("logcounts",
                              "cpm", "logcpm")),
                  withBusyIndicatorUI(actionButton("addAssay", "Add Assay")),
                  selectInput("delAssayType", "Delete Assay Type:", currassays),
                  withBusyIndicatorUI(actionButton("delAssay", "Delete Assay")),
                  tags$hr(),
                  h4("reducedDim Options:"),
                  selectInput("delRedDimType", "Delete reducedDim:", currreddim),
                  withBusyIndicatorUI(actionButton("delRedDim", "Delete reducedDim"))
                ),
                mainPanel(
                  fluidRow(
                    column(6,
                           h4("Available Assays:"),
                           tableOutput("assayList")
                    ),
                    column(6,
                           h4("Available Reduced Dims:"),
                           tableOutput("reducedDimsList")
                    )
                  )
                )
              )
            )
          ),
          tabPanel(
            "Annotation Data",
            tabsetPanel(
              tabPanel(
                "Data",
                DT::dataTableOutput("colDataDataFrame")
              ),
              tabPanel(
                "Options",
                wellPanel(
                  h4("Modify Annotation Data:"),
                  selectInput("annotModifyChoice", "Select Annotation:", c("none", clusterChoice)),
                  uiOutput("annotModifyUI"),
                  tags$hr(),
                  downloadButton("downloadcolData", "Download Annotation Data"),
                  tags$hr(),
                  fileInput(
                    "newAnnotFile", "Upload and replace annotation data:",
                    accept = c(
                      "text/csv",
                      "text/comma-separated-values",
                      ".csv"
                    )
                  )
                )
              )
            )
          )
        )
      ),
      position = "right"
    )
  )
)
