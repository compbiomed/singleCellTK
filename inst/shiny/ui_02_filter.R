shinyPanelFilter <- fluidPage(
  useShinyalert(),
  tags$div(
    class = "container",
    h1("Data Summary & Filtering"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v04-tab02_Data-Summary-and-Filtering.html",
              "(help)", target = "_blank")),

    tabsetPanel(
      tabPanel(
        "Data Summary",
        wellPanel(
          sidebarLayout(
            sidebarPanel(
              fluidRow(
                column(5,
                  h3("Settings:")
                ),
                column(3,
                  br(),
                  actionButton("f_hideAllSections", "Hide All")
                ),
                column(3,
                  br(),
                  actionButton("f_showAllSections", "Show All")
                )
              ),
              br(),
              
              # SHINYJS COLLAPSE --------------------------
              # Section 1 - Assay Settings
              actionButton("f_button1", "Assay Settings"),
              # open by default
              tags$div( id="f_collapse1",
                wellPanel(
                  selectInput("filterAssaySelect", "Select Assay:", currassays)
                )
              ),
              # Section 2 - Delete Outliers
              actionButton("f_button2", "Delete Outliers"),
              shinyjs::hidden(
                tags$div( id="f_collapse2",
                  wellPanel(
                    checkboxInput("removeNoexpress", "Remove genes with 0 expression across all samples (Recommended)", value = TRUE),
                    numericInput("minDetectGene", label = "Minimum Detected Genes per Sample.", value = 1700, min = 1, max = 100000),
                    numericInput("LowExpression", "% Low Gene Expression to Filter", value = 40, min = 0, max = 100),
                    selectInput("deletesamplelist", "Select Samples:", sampleChoice, multiple = TRUE),
                    fluidRow(
                      column(6, withBusyIndicatorUI(actionButton("filterData", "Filter Data"))),
                      column(6, actionButton("resetData", "Reset"))
                    )
                  )
                )
              ),
              # Section 3 - Filter Samples by Annotation
              actionButton("f_button3", "Filter Samples by Annotation"),
              shinyjs::hidden(
                tags$div( id="f_collapse3",
                  wellPanel(
                    selectInput("filteredSample", "Select Annotation:", c("none", clusterChoice)),
                    uiOutput("filterSampleOptions")
                  )
                )
              ),
              # Section 4 - Filter Genes by Feature Annotation
              actionButton("f_button4", "Filter Genes by Feature Annotation"),
              shinyjs::hidden(
                tags$div( id="f_collapse4",
                  wellPanel(
                    selectInput("filteredFeature", "Select Feature:", c("none", featureChoice)),
                    uiOutput("filterFeatureOptions")
                  )
                )
              ),
              # Section 5 - Convert Gene Annotations
              actionButton("f_button5", "Convert Gene Annotations"),
              shinyjs::hidden(
                tags$div( id="f_collapse5",
                  wellPanel(  
                    selectInput("orgOrganism", "Select Organism:", 
                      as.character(grep("^org\\.", installed.packages()[, "Package"], value = TRUE))
                    ),
                    uiOutput("orgConvertColumns"),
                    withBusyIndicatorUI(actionButton("convertGenes", "Convert"))
                  )
                )
              ),
              # Section 6 - Delete an Annotation Column
              actionButton("f_button6", "Delete an Annotation Column"),
              shinyjs::hidden(
                tags$div( id="f_collapse6",
                  wellPanel(
                    selectInput("deleterowdatacolumn", "Annotation Column:", clusterChoice),
                    actionButton("deleterowDatabutton", "Delete Column")
                  )
                )
              ),
              # Section 7 - Randomly Subset
              actionButton("f_button7", "Randomly Subset"),
              shinyjs::hidden(
                tags$div( id="f_collapse7",
                  wellPanel(
                    numericInput("downsampleNum", "Number of samples to keep:", min = 2,
                                 max = numSamples, value = numSamples, step = 1),
                    withBusyIndicatorUI(actionButton("downsampleGo", "Subset Data"))
                  )
                )
              ),
              tags$hr(),
              downloadButton("downloadSCE", "Download SCtkExperiment")
            ),
            mainPanel(
              wellPanel(
                style = "background-color:transparent",
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
            )
          )
        )
      ),
      tabPanel(
        "Assay Details",
        wellPanel(
          sidebarLayout(
            sidebarPanel(
              h3("Assay Options:"),
              selectInput("addAssayType", "Add Assay Type:", c("logcounts",
                                                               "cpm", "logcpm")),
              withBusyIndicatorUI(actionButton("addAssay", "Add Assay")),
              selectInput("delAssayType", "Delete Assay Type:", currassays),
              withBusyIndicatorUI(actionButton("delAssay", "Delete Assay"))
            ),
            mainPanel(
              h4("Available Assays:"),
              tableOutput("assayList")
            )
          )
        )
      ),
      tabPanel(
        "Annotation Data",
        wellPanel(
          sidebarLayout(
            sidebarPanel(
              h3("Modify Annotation Data:"),
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
            ),
            mainPanel(
              tags$h4("Data:"),
              tags$br(),
              DT::dataTableOutput("colDataDataFrame")
            )
          )
        )
      )
    )
  )
)
