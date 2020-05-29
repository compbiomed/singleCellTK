shinyPanelMASTDE <- fluidPage(
  tags$div(
    class = "container",
    h1("MAST Differential Expression"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v07-tab05_Differential-Expression.html#mast",
              "(help)", target = "_blank")),
    fluidRow(
      panel(style = "margin:2px;",
        selectInput("mastAssay", "Select Assay:", currassays),
        radioButtons('mastCondMethod', "Condition Selection:",
                     choiceNames = c("By annotations",
                                     "Manually select individual cells",
                                     "Manually enter cell IDs"),
                     choiceValues = c(1, 2, 3), inline = TRUE),
        fluidRow(
          column(width = 6,
            textInput("mastG1Name", "Name of Condition1", "Condition1",
                      placeholder = "Required")
          ),
          column(width = 6,
            textInput("mastG2Name", "Name of Condition2", "Condition2",
              placeholder = "Required")
          )
        ),
        conditionalPanel(
          condition = "input.mastCondMethod == 1",
          # Select by class UI ####
          panel(
            selectInput("mastC1Class", "Choose Annotation Class:",
              clusterChoice),
            fluidRow(
              column(width = 6,
                uiOutput("mastC1G1UI"),
                uiOutput("mastC1G1CellCheckUI"),
                uiOutput("mastC1G1NCell")
              ),
              column(width = 6,
                uiOutput("mastC1G2UI"),
                uiOutput("mastC1G2CellCheckUI"),
                uiOutput("mastC1G2NCell")
              )
            )
          )
        ),
        conditionalPanel(
          condition = "input.mastCondMethod == 2",
          panel(
            tabsetPanel(
              tabPanel("Condition 1",
                tagList(
                  selectInput('mastC2G1Col',
                    "Columns to display (Changes here clear the selection)",
                    clusterChoice, multiple = TRUE, width = '550px'),
                  DT::dataTableOutput("mastC2G1Table"),
                  actionButton('mastC2G1Table_addAll', "Add all filtered"),
                  actionButton('mastC2G1Table_clear', "Clear selection"),
                )
              ),
              tabPanel("Condition 2",
                tagList(
                  selectInput('mastC2G2Col',
                    "Columns to display (Changes here clear the selection)",
                    clusterChoice, multiple = TRUE, width = '550px'),
                  DT::dataTableOutput("mastC2G2Table"),
                  actionButton('mastC2G2Table_addAll', "Add all filtered"),
                  actionButton('mastC2G2Table_clear', "Clear selection"),
                )
              )
            ),
            h4("Summary", style = 'margin-top:10px'),
            uiOutput("mastC2G1info"),
            uiOutput("mastC2G2info")
          )
        ),
        conditionalPanel(
          condition = "input.mastCondMethod == 3",
          # Direct enter name UI ####
          panel(
            fluidRow(
              column(width = 6,
                h4("The subject group:"),
                textAreaInput("mastC3G1Cell", "Cell IDs:", height = '150px',
                  placeholder = "Enter cell IDs here, \none per line with no symbol separator."),
                uiOutput("mastC3G1NCell")
              ),
              column(width = 6,
                h4("The object group:"),
                textAreaInput("mastC3G2Cell", "Cell IDs:", height = '150px',
                  placeholder = "Leave this blank for the all the rest cells."),
                uiOutput("mastC3G2NCell")
              )
            )
          )
        ),
        h4("Parameters:"),
        fluidRow(style = 'margin:4px;',
          div(style="display:inline-block;vertical-align:center;width:230px;margin-right:8px;",
            numericInput("mastFreq", "Use gene expressed in more than:",
              min = 0, max = 1, step = 0.05, value = 0.1),),
          div(style="display: inline-block;vertical-align:center; width: 160px;margin-left:8px;margin-right:8px;",
            numericInput("mastFDRThresh", "Output FDR less than:",
              min = 0.01, max = 1, step = 0.01, value = 0.5)),
          div(style="display: inline-block;vertical-align:center; width: 300px;margin-left:8px;margin-right:8px;",
            numericInput("mastFCThresh", "Output Log2FC Absolute value greater than:",
              min = 0, step = 0.05, value = 0)),
          div(style="display: inline-block;vertical-align:center; width: 230px;margin-left:8px;margin-right:8px;",
            checkboxInput("mastPosOnly", "Only output up-regulated genes",
              value = FALSE)),
          div(style="display: inline-block;vertical-align:top; width: 200px;margin-left:8px;",
            checkboxInput("useAdaptThresh", "Use Adaptive Thresholds",
              value = FALSE))
        ),
        fluidRow(style = 'margin:4px;',
          div(style="display: inline-block;vertical-align:center; width: 150px;margin-right:8px;",
            uiOutput("mastCompNameUI")),
          div(style="display: inline-block;vertical-align:center; width: 150px;margin-left:8px;margin-right:8px",
            withBusyIndicatorUI(actionButton("runMAST", "Run MAST"))),
          div(style="display: inline-block;vertical-align:center; width: 360px;margin-left:8px;",
            uiOutput("mastNameWarn"))
        )
      )
    ),
    fluidRow(
      uiOutput("mastResSelUI"),
      tabsetPanel(
        tabPanel("Adaptive thresholding", plotOutput("threshplot")),
        tabPanel("Results Table",
          DT::dataTableOutput("mastresults"),
          downloadButton("mastDownload", "Download Results", height = "800px")),
        tabPanel("Violin Plot", plotOutput("hurdleviolin", height = "800px")),
        tabPanel("Linear Model", plotOutput("hurdlelm", height = "800px")),
        tabPanel("Heatmap", plotOutput("hurdleHeatmap", height = "800px"))
      )
    )
  )
)
