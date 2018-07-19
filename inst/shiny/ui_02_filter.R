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
          br(),
          fluidRow(
          sidebarPanel(
            h3("Settings:"),
            
            # SHINYJS ACCORDION --------------------------
            # Section 1 - Assay Settings
            actionButton("button1", "Assay Settings"),
            # collapse open by default
            tags$div( id="collapse1",
              wellPanel(
                selectInput("filterAssaySelect", "Select Assay:", currassays),
                checkboxInput("removeNoexpress", "Remove genes with 0 expression across all samples (Recommended)", value = TRUE),
                numericInput("minDetectGene", label = "Minimum Detected Genes per Sample.", value = 1700, min = 1, max = 100000),
                numericInput("LowExpression", "% Low Gene Expression to Filter", value = 40, min = 0, max = 100)
              )
            ),
            # Section 2 - Delete Outliers
            actionButton("button2", "Delete Outliers"),
            shinyjs::hidden(
              tags$div( id="collapse2",
                wellPanel(
                  selectInput("deletesamplelist", "Select Samples:", sampleChoice, multiple = TRUE),
                  fluidRow(
                    column(6, withBusyIndicatorUI(actionButton("filterData", "Filter Data"))),
                    column(6, actionButton("resetData", "Reset"))
                  )
                )
              )
            ),
            # Section 3 - Filter Samples by Annotation
            actionButton("button3", "Filter Samples by Annotation"),
            shinyjs::hidden(
              tags$div( id="collapse3",
                wellPanel(
                  selectInput("filteredSample", "Select Annotation:", c("none", clusterChoice)),
                  uiOutput("filterSampleOptions")
                )
              )
            ),
            # Section 4 - Filter Genes by Feature Annotation
            actionButton("button4", "Filter Genes by Feature Annotation"),
            shinyjs::hidden(
              tags$div( id="collapse4",
                wellPanel(
                  selectInput("filteredFeature", "Select Feature:", c("none", featureChoice)),
                  uiOutput("filterFeatureOptions")
                )
              )
            ),
            # Section 5 - Convert Gene Annotations
            actionButton("button5", "Convert Gene Annotations"),
            shinyjs::hidden(
              tags$div( id="collapse5",
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
            actionButton("button6", "Delete an Annotation Column"),
            shinyjs::hidden(
              tags$div( id="collapse6",
                wellPanel(
                  selectInput("deleterowdatacolumn", "Annotation Column:", clusterChoice),
                  actionButton("deleterowDatabutton", "Delete Column")
                )
              )
            ),
            # Section 7 - Randomly Subset
            actionButton("button7", "Randomly Subset"),
            shinyjs::hidden(
              tags$div( id="collapse7",
                wellPanel(
                  numericInput("downsampleNum", "Number of samples to keep:", min = 2,
                               max = numSamples, value = numSamples, step = 1),
                  withBusyIndicatorUI(actionButton("downsampleGo", "Subset Data"))
                )
              )
            ),
            
            # # HTML ACCORDION -----------------------------
            # HTML('<div class="accordion" id="filterAccordion">
            #   <div class="panel" style="background-color:transparent">'),
            #     # section format - accordionSelection(collapseId, accordionID, sectionTitle) from helpers.R
            #     HTML(accordionSection("collapse-AssaySettings", "Assay Settings", "filterAccordion")),
            #       wellPanel(  
            #         selectInput("filterAssaySelect", "Select Assay:", currassays),
            #         checkboxInput("removeNoexpress", "Remove genes with 0 expression across all samples (Recommended)", value = TRUE),
            #         numericInput("minDetectGene", label = "Minimum Detected Genes per Sample.", value = 1700, min = 1, max = 100000),
            #         numericInput("LowExpression", "% Low Gene Expression to Filter", value = 40, min = 0, max = 100)
            #       ),
            #     HTML('</div>'),
            # 
            #     HTML(accordionSection("collapse-DeleteOutliers", "Delete Outliers", "filterAccordion")),
            #       wellPanel(  
            #         selectInput("deletesamplelist", "Select Samples:",
            #           sampleChoice,
            #           multiple = TRUE),
            #         fluidRow(
            #           column(6, withBusyIndicatorUI(actionButton("filterData", "Filter Data"))),
            #           column(6, actionButton("resetData", "Reset"))
            #         )
            #       ),
            #     HTML('</div>'),
            # 
            #     HTML(accordionSection("collapse-FilterSamples", "Filter samples by annotation", "filterAccordion")),
            #       wellPanel(
            #         selectInput("filteredSample", "Select Annotation:", c("none", clusterChoice)),
            #         uiOutput("filterSampleOptions")
            #       ),
            #     HTML('</div>'),
            # 
            #     HTML(accordionSection("collapse-FilterGenes", "Filter genes by feature annotation", "filterAccordion")),
            #       wellPanel(
            #         selectInput("filteredFeature", "Select Feature:", c("none", featureChoice)),
            #         uiOutput("filterFeatureOptions")
            #       ),
            #     HTML('</div>'),
            # 
            #     HTML(accordionSection("collapse-ConvertAnnotations", "Convert gene annotations", "filterAccordion")),
            #       wellPanel(  
            #         selectInput("orgOrganism", "Select Organism:", as.character(grep("^org\\.",
            #           installed.packages()[, "Package"], value = TRUE))),
            #         uiOutput("orgConvertColumns"),
            #         withBusyIndicatorUI(actionButton("convertGenes", "Convert"))
            #       ),
            #     HTML('</div>'),
            # 
            #     HTML(accordionSection("collapse-DeleteColumn", "Delete an annotation column", "filterAccordion")),
            #       wellPanel(
            #         selectInput("deleterowdatacolumn", "Annotation Column:", clusterChoice),
            #         actionButton("deleterowDatabutton", "Delete Column")
            #       ),
            #     HTML('</div>'),
            # 
            #     HTML(accordionSection("collapse-RandomlySubset", "Randomly Subset", "filterAccordion")),
            #       wellPanel(
            #         numericInput("downsampleNum", "Number of samples to keep:", min = 2,
            #           max = numSamples, value = numSamples, step = 1),
            #         withBusyIndicatorUI(actionButton("downsampleGo", "Subset Data"))
            #       ),
            #     HTML('</div>'),
            # 
            #   HTML('</div>
            # </div>'),
            
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
          br(),
          fluidRow(
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
      tabPanel(
        "Visualize",
        wellPanel(
          br(),
          fluidRow(
            sidebarLayout(
              sidebarPanel(
                h4("Visualization Options:"),
                selectInput("visAssaySelect", "Select Assay:", currassays),
                selectInput("visPlotMethod", "Visualization Method:", c("boxplot", "scatterplot", "barplot", "heatmap")),
                selectInput("visCondn", "Condition:", c(clusterChoice)),
                selectizeInput("selectvisGenes", label = "Select Gene(s):", NULL, multiple = TRUE),
                withBusyIndicatorUI(actionButton("plotvis", "Plot"))
              ),
              mainPanel(
                fluidRow(
                  plotOutput("visPlot")
                )
              )
            )
          )
        )
      )
    )
  )
)
