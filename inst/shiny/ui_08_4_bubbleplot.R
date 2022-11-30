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
                      actionButton("bpImportRun", "Import"))
      ),
      
      hr(),
      # Subset ####
      h3("Cell/Feature Subsetting"),
      p("Only to plot cells/features of interests", style = "color:grey;"),
      tabsetPanel(
        id = 'bpSubsetTSP',
        tabPanel(
          title = "Cell", value = 'bpSubsetCellTP',
          tagList(
            uiOutput('bpCellColUI'),
            DT::dataTableOutput("bpCellColTable"),
            actionButton('bpCellColTable_addAll', "Add all filtered"),
            actionButton('bpCellColTable_clear', "Clear selection"),
          )
        ),
        tabPanel(
          title = "Feature", value = 'hmSubsetGeneTP',
          tagList(
            uiOutput('bpGeneColUI'),
            DT::dataTableOutput("bpGeneColTable"),
            actionButton('bpGeneColTable_addAll', "Add all filtered"),
            actionButton('bpGeneColTable_clear', "Clear selection"),
          )
        ),
      ),
      uiOutput("bpCellSumUI"),
      uiOutput("bpGeneSumUI"),
      hr(),
      
      # Others ####
      h3("Bubbleplot Setting"),
      p("Settings for title, label, color scheme and etc.",
        style = "color:grey;"),
      panel(
        fluidRow(
          column(width = 4,
                 textInput("bpTitle", "Title", NULL)
          )
        ),
        fluidRow(
          column(width = 4,
                 textInput("bpX", "X-axis Label", NULL)
          ),
          column(width = 4,
                 textInput("bpY", "Y-axis Label", NULL)
          )
        ),
        fluidRow(
          column(
            width = 4,
            colourpicker::colourInput('bpLow', 'Low color',value = 'white')
          ),
          column(
            width = 4,
            colourpicker::colourInput('bpHigh', 'High color',value = 'blue')
          )
        )
      ),
      hr(),
      withBusyIndicatorUI(actionButton("plotBubbleplot", "Plot Bubbleplot")),
      div(
        style = 'height:800px;',
        plotOutput("Bubbleplot")
      )
    )
  )
)