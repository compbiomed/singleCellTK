shinyPanelBubbleplot <- fluidPage(
  tags$div(
    class = "container",
    h1("Bubbleplot"),
    p("Generic bubbleplot plotting panel for customized figure.",
      style = "color:grey;"),
    h5(tags$a(href = paste0(docs.artPath, "bubbleplot.html"),
              "(help)", target = "_blank")),
    panel(
      h3("Assay to Plot"),
      fluidRow(
        column(width = 4,
               selectizeInput(
                 inputId = "hmAssay", 
                 label = "Select input matrix:", 
                 choices = NULL, 
                 selected = NULL, 
                 multiple = FALSE,
                 options = NULL)
               
               
        )
      ),
      fluidRow(column(width = 4,
                      actionButton("hmImportRun", "Import"))
               ),
      
      hr(),
      # Subset ####
      h3("Cell/Feature Subsetting"),
      p("Only to plot cells/features of interests", style = "color:grey;"),
      tabsetPanel(
        id = 'hmSubsetTSP',
        tabPanel(
          title = "Cell", value = 'hmSubsetCellTP',
          tagList(
            uiOutput('hmCellColUI'),
            DT::dataTableOutput("hmCellColTable"),
            actionButton('hmCellColTable_addAll', "Add all filtered"),
            actionButton('hmCellColTable_clear', "Clear selection"),
          )
        ),
        tabPanel(
          title = "Feature", value = 'hmSubsetGeneTP',
          tagList(
            uiOutput('hmGeneColUI'),
            DT::dataTableOutput("hmGeneColTable"),
            actionButton('hmGeneColTable_addAll', "Add all filtered"),
            actionButton('hmGeneColTable_clear', "Clear selection"),
          )
        ),
      ),
      uiOutput("hmCellSumUI"),
      uiOutput("hmGeneSumUI"),
      hr(),
      
      # Others ####
      h3("Bubbleplot Setting"),
      p("Settings for title, label, color scheme and etc.",
        style = "color:grey;"),
          panel(
            fluidRow(
              column(
                width = 4,
                dropdown(
                  panel(
                    fluidRow(
                      column(width = 6,
                             textInput("deG1Name", "Title", NULL)
                      )
                    ),
                    fluidRow(
                      column(width = 6,
                             textInput("deG2Name", "X-axis Label", NULL)
                      ),
                      column(width = 6,
                             textInput("deG2Name", "Y-axis Label", NULL)
                      )
                    ),
                    fluidRow(
                      column(
                        width = 4,
                        colourpicker::colourInput('hmCSLow', 'Low color',value = 'white')
                      ),
                      column(
                        width = 4,
                        colourpicker::colourInput('hmCSHigh', 'High color',value = 'blue')
                      )
                    )
                  ),
                  inputId = "dropDownDeHM",
                  icon = icon("cog"),
                  status = "primary",
                  circle = FALSE,
                  inline = TRUE,
                  width = "500px"
                )
              )
            )
          ),
      panel(
        fluidRow(
          column(width = 4,
                 textInput("deG1Name", "Title", NULL)
          )
        ),
        fluidRow(
          column(width = 4,
                 textInput("deG2Name", "X-axis Label", NULL)
          ),
          column(width = 4,
                 textInput("deG2Name", "Y-axis Label", NULL)
          )
        ),
        fluidRow(
          column(
            width = 4,
            colourpicker::colourInput('hmCSLow', 'Low color',value = 'white')
          ),
          column(
            width = 4,
            colourpicker::colourInput('hmCSHigh', 'High color',value = 'blue')
          )
        )
      ),
      hr(),
      withBusyIndicatorUI(actionButton("plotHeatmap", "Plot Bubbleplot")),
      div(
        style = 'height:800px;',
        plotOutput("Heatmap")
      )
    )
  )
)

