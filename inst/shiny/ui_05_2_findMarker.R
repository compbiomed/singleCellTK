shinyPanelfindMarker <- fluidPage(
  tags$script("Shiny.addCustomMessageHandler('close_dropDownFM', function(x){
                  $('html').click();
                });"),
    h1("Find Marker"),
    h5(tags$a(href = paste0(docs.artPath, "find_marker.html"),
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        selectInput('fmMethod', "Select Differential Expression Method",
                    c("wilcox", "MAST", "DESeq2", "Limma", "ANOVA")),
        selectizeInput(
          inputId = "fmAssay", 
          label = "Select input matrix:", 
          choices = NULL, 
          selected = NULL, 
          multiple = FALSE,
          options = NULL),
        selectInput("fmCluster", "Cluster Annotation", clusterChoice),
        selectInput("fmCovar", "Covariate(s)", clusterChoice, multiple = TRUE),
        # numericInput("fmLogFC", "Log2FC greater than",
        #              value = 0, min = 0, step = 0.05),
        numericInput("fmFDR", "Keep FDR less than", value = 0.05,
                     max = 1, min = 0.01, step = 0.01),
        # numericInput("fmMinClustExprPerc",
        #              "Minimum Expressed Percentage in Cluster",
        #              value = 0, min = 0, max = 1),
        # numericInput("fmMaxCtrlExprPerc",
        #              "Maximum Expressed Percentage in Control",
        #              value = 1, min = 0, max = 1),
        # numericInput("fmMinMeanExpr",
        #              "Minimum Mean Expression Level in Cluster",
        #              value = 0, min = 0),
        actionButton("runFM", "Run")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Heatmap",
            shinyjs::useShinyjs(),
            fluidRow(
              column(
                width = 4,
                dropdown(
                  fluidRow(
                    column(12,
                           fluidRow(actionBttn(inputId = "closeDropDownFM", label = NULL, style = "simple", color = "danger", icon = icon("times"), size = "xs"), align = "right"),
                    )
                  ),
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
                      withBusyIndicatorUI(
                        actionBttn(
                        inputId = "plotFM",
                        label = "Update",
                        style = "bordered",
                        color = "primary",
                        size = "sm"
                      )
                      )
                    ),
                    column(
                      width = 6,
                      selectInput("fmHMcolData", "Additional cell annotation",
                                  clusterChoice, multiple = TRUE)
                    )
                  ),
                  inputId = "dropDownFM",
                  icon = icon("cog"),
                  status = "primary",
                  circle = FALSE,
                  inline = FALSE,
                  width = "500px"
                )
              ),
              column(
                width = 7,
                fluidRow(
                  h6(
                    "The heatmap plots the expression level of top markers found for each cluster"
                  ),
                  align="center"
                )
              )
            ),
            hr(),
            br(),
            shinyjqui::jqui_resizable(plotOutput('fmHeatmap'))
          ),
          tabPanel(
            "Result Table",
            filterTableUI(id = "filterfmResTable"),
            #DT::dataTableOutput("fmResTable"),
            downloadButton("fmDownload", "Download Results"),
          )
        )
      )
    )
)

