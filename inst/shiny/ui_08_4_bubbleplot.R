shinyPanelBubbleplot <- fluidPage(
  tags$div(
    class = "container",
    h1("Bubbleplot"),
    p("Generic bubbleplot plotting panel for customized figure.",
      style = "color:grey;"),
    panel(
      
      fluidRow(
        column(width = 4,
               hr(),
               # Subset ####
               h3("Cluster"),
               uiOutput("bpClusterUI"),
               h3("Feature"),
               fluidRow(
                 column(width = 12,
                        selectizeInput(
                          inputId = "bpAssay", 
                          label = "Select input matrix:", 
                          choices = NULL, 
                          selected = NULL, 
                          multiple = FALSE,
                          options = NULL)
                 )
               ),
               uiOutput('bpRowUI'),
               selectizeInput(
                 'bpFeatures',
                 "Select Features",
                 choices = NULL, multiple = TRUE, width = '550px'),
               h5(style="display: inline-block; margin-top: 0px; margin-bottom: 20px","Scale Data"),
               switchInput(
                 inputId = "scaleBubble",
                 onLabel = "Yes",
                 offLabel = "No",
                 value=FALSE,
                 size="mini",
                 inline = TRUE
               ),
               withBusyIndicatorUI(actionButton("plotBubbleplot", "Plot Bubbleplot")),
               hr(),
        ),
        
        column(
          width = 3,
          dropdown(
            fluidRow(
              column(
                width = 12,
                fluidRow(actionBttn(inputId = "closeDropDownBubble", label = NULL, style = "simple", color = "danger", icon = icon("times"), size = "xs"), align = "right"),
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
                withBusyIndicatorUI(
                  actionBttn(
                    inputId = "updateBubbleplot",
                    label = "Update",
                    style = "bordered",
                    color = "primary",
                    size = "sm"
                  )
                ),
              )
            ),
            inputId = "dropDownBubble",
            icon = icon("cog"),
            status = "primary",
            circle = FALSE,
            inline = TRUE
          )
        ),
        
        column(
          8,
          wellPanel(
            h5(strong("Plotting Region")),
            div(
              style = 'height:800px;',
              plotOutput("Bubbleplot")
            ),
            hr(),
          )
        )
      ),
      
    )
  )
)
