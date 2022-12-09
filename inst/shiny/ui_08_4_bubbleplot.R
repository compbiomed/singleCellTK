shinyPanelBubbleplot <- fluidPage(
  tags$div(
    class = "container",
    h1("Bubbleplot"),
    p("Generic bubbleplot plotting panel for customized figure.",
      style = "color:grey;"),
    panel(
      h3("Assay to Plot"),
      fluidRow(
        column(width = 4,
               selectizeInput(
                 inputId = "bpAssay", 
                 label = "Select input matrix:", 
                 choices = NULL, 
                 selected = NULL, 
                 multiple = FALSE,
                 options = NULL)
        )
      ),
      
      hr(),
      # Subset ####
      h3("Cluster/Feature Subsetting"),
      p("Select cluster and features of interests", style = "color:grey;"),
      tabsetPanel(
        tabPanel(
          title = "Cluster",
          uiOutput("bpClusterUI")),
        tabPanel(
          title = "Feature",
          uiOutput('bpRowUI'),
          selectizeInput(
            'bpFeatures',
            "Select Features",
            choices = NULL, multiple = TRUE, width = '550px')
        ),
      ),
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
            colourpicker::colourInput('bpLow', 'Low color', value = 'white')
          ),
          column(
            width = 4,
            colourpicker::colourInput('bpHigh', 'High color', value = 'blue')
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
