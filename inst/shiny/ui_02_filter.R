shinyPanelFilter <- fluidPage(
  useShinyalert(),
  tags$div(
    class = "container",
    h1("Data Summary & Filtering"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v04-tab02_Data-Summary-and-Filtering.html",
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        h3("Settings:"),
          HTML('<div class="accordion" id="filterAccordion">
            <div class="panel" style="background-color:transparent">'),
            # section format - accordionSelection(collapseId, accordionID, sectionTitle) from helpers.R
        
            HTML(accordionSection("collapse-AssaySettings", "Assay Settings", "filterAccordion")),
              selectInput("filterAssaySelect", "Select Assay:", currassays),
              checkboxInput("removeNoexpress", "Remove genes with 0 expression across all samples (Recommended)", value = TRUE),
              numericInput("minDetectGene", label = "Minimum Detected Genes per Sample.", value = 1700, min = 1, max = 100000),
              numericInput("LowExpression", "% Low Gene Expression to Filter", value = 40, min = 0, max = 100),
            HTML('</div>'),
            
            HTML(accordionSection("collapse-DeleteOutliers", "Delete Outliers", "filterAccordion")),
              selectInput("deletesamplelist", "Select Samples:",
                sampleChoice,
                multiple = TRUE),
              withBusyIndicatorUI(actionButton("filterData", "Filter Data")),
              actionButton("resetData", "Reset"),
            HTML('</div>'),
            
            HTML(accordionSection("collapse-FilterSamples", "Filter samples by annotation", "filterAccordion")),
              selectInput("filteredSample", "Select Annotation:", c("none", clusterChoice)),
              uiOutput("filterSampleOptions"),
            HTML('</div>'),
            
            HTML(accordionSection("collapse-FilterGenes", "Filter genes by feature annotation", "filterAccordion")),
              selectInput("filteredFeature", "Select Feature:", c("none", featureChoice)),
              uiOutput("filterFeatureOptions"),
            HTML('</div>'),
            
            HTML(accordionSection("collapse-ConvertAnnotations", "Convert gene annotations", "filterAccordion")),
              selectInput("orgOrganism", "Select Organism:", as.character(grep("^org\\.", 
                installed.packages()[, "Package"], value = TRUE))
              ),
              uiOutput("orgConvertColumns"),
              withBusyIndicatorUI(actionButton("convertGenes", "Convert")),
            HTML('</div>'),
            
            HTML(accordionSection("collapse-DeleteColumn", "Delete an annotation column", "filterAccordion")),
              selectInput("deleterowdatacolumn", "Annotation Column:", clusterChoice),
              actionButton("deleterowDatabutton", "Delete Column"),
            HTML('</div>'),
            
            HTML(accordionSection("collapse-RandomlySubset", "Randomly Subset", "filterAccordion")),
              numericInput("downsampleNum", "Number of samples to keep:", min = 2,
                max = numSamples, value = numSamples, step = 1),
              withBusyIndicatorUI(actionButton("downsampleGo", "Subset Data")),
            HTML('</div>'),
        
          HTML('</div>
        </div>'), 
        tags$hr(),
        downloadButton("downloadSCE", "Download SCtkExperiment")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Data Summary",
            wellPanel( 
              style="background-color:transparent",
              h4("Summary Contents:"),
              tableOutput("summarycontents"),
              tags$hr(),
              h4("Counts Histogram:"),
              plotlyOutput("countshist"),
              tags$hr(),
              h4("Genes Histogram:"),
              plotlyOutput("geneshist"),
              tags$hr(),
              h4("Data Table:"),
              DT::dataTableOutput("contents")
            )
          ),
          tabPanel(
            "Assay Details",
            wellPanel(  
              style="background-color:transparent",
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
            )
          ),
          tabPanel(
            "Annotation Data",
            wellPanel( 
              style="background-color:transparent",
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
        )
      )
    )
  )
)
