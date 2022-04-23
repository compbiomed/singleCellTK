shinyPanelDiffex <- fluidPage(
  tags$script("Shiny.addCustomMessageHandler('close_dropDownDeHM', function(x){
                  $('html').click();
                });"),
  tags$script("Shiny.addCustomMessageHandler('close_dropDownDeViolin', function(x){
                  $('html').click();
                });"),
  tags$script("Shiny.addCustomMessageHandler('close_dropDownDeReg', function(x){
                  $('html').click();
                });"),
  tags$div(
    class = "container",
    h1("Differential Expression"),
    h5(tags$a(href = paste0(docs.artPath, "differential_expression.html"),
              "(help)", target = "_blank")),
    fluidRow(
      panel(
        style = "margin:2px;",
        h3("Method and Matrix"),
        fluidRow(
          column(
            4,
            selectInput('deMethod', "Choose analysis method",
                        c('wilcox', 'MAST', 'DESeq2', 'Limma', 'ANOVA'))
          ),
          column(
            4,
            selectizeInput(
              inputId = "deAssay", 
              label = "Select input matrix:", 
              choices = NULL, 
              selected = NULL, 
              multiple = FALSE,
              options = NULL)
            #uiOutput("deAssay")
          )
        ),
        useShinyjs(),
        actionButton("deViewThresh", label = "View Thresholding"),
        shinyjs::hidden(
          wellPanel(
            id = "deThreshpanel",
            textOutput("deSanityWarnThresh"),
            uiOutput("deThreshPlotDiv"),
            actionButton("deHideThresh", label = "Hide")
            )
          ),
        h3("Condition Setting"),
        p("Three approaches of setting provided for flexibility. ",
          style = "color:grey;"),
        radioButtons('deCondMethod', "Condition Selection:",
                     choiceNames = c("By annotations",
                                     "Manually select individual cells",
                                     "Manually enter cell IDs"),
                     choiceValues = c(1, 2, 3), inline = TRUE),
        fluidRow(
          column(width = 6,
                 textInput("deG1Name", "Name of Condition1", NULL,
                           placeholder = "Required")
          ),
          column(width = 6,
                 textInput("deG2Name", "Name of Condition2", NULL,
                           placeholder = "Required")
          )
        ),
        conditionalPanel(
          condition = "input.deCondMethod == 1",
          # Select by class UI ####
          panel(
            selectInput("deC1Class", "Choose Annotation Class:",
                        clusterChoice),
            fluidRow(
              column(width = 6,
                     uiOutput("deC1G1UI"),
                     uiOutput("deC1G1CellCheckUI"),
                     uiOutput("deC1G1NCell")
              ),
              column(width = 6,
                     uiOutput("deC1G2UI"),
                     uiOutput("deC1G2CellCheckUI"),
                     uiOutput("deC1G2NCell")
              )
            )
          )
        ),
        conditionalPanel(
          condition = "input.deCondMethod == 2",
          panel(
            tabsetPanel(
              tabPanel(
                "Condition 1",
                tagList(
                  selectInput(
                    'deC2G1Col',
                    "Columns to display",
                    clusterChoice, multiple = TRUE, width = '550px'
                  ),
                  DT::dataTableOutput("deC2G1Table"),
                  actionButton('deC2G1Table_addAll', "Add all filtered"),
                  actionButton('deC2G1Table_clear', "Clear selection"),
                )
              ),
              tabPanel(
                "Condition 2",
                tagList(
                  selectInput(
                    'deC2G2Col',
                    "Columns to display",
                    clusterChoice, multiple = TRUE, width = '550px'
                  ),
                  p("Leave unselected for all the others.",
                    style = 'color:grey;'),
                  DT::dataTableOutput("deC2G2Table"),
                  actionButton('deC2G2Table_addAll', "Add all filtered"),
                  actionButton('deC2G2Table_clear', "Clear selection"),
                )
              )
            ),
            h4("Summary", style = 'margin-top:10px'),
            uiOutput("deC2G1info"),
            uiOutput("deC2G2info")
          )
        ),
        conditionalPanel(
          condition = "input.deCondMethod == 3",
          # Direct enter name UI ####
          panel(
            fluidRow(
              column(width = 6,
                     h4("The Condition of Interests:"),
                     textAreaInput(
                       "deC3G1Cell", "Cell IDs:", height = '150px',
                       placeholder = "Enter cell IDs here, \none per line with no symbol separator."),
                     uiOutput("deC3G1NCell")
              ),
              column(width = 6,
                     h4("The Control Condition:"),
                     textAreaInput(
                       "deC3G2Cell", "Cell IDs:", height = '150px',
                       placeholder = "Leave this blank for all the others."),
                     uiOutput("deC3G2NCell")
              )
            )
          )
        ),
        h3("Parameters"),
        fluidRow(
          column(
            width = 3,
            selectInput("deCovar", "Select Covariates",
                        clusterChoice, multiple = TRUE)
          ),
          column(
            width = 3,
            numericInput("deFDRThresh", "Output FDR less than:",
                         min = 0, max = 1, step = 0.01, value = 0.05)
          ),
          column(
            width = 3,
            numericInput("deFCThresh",
                         "Output Log2FC Absolute value greater than:",
                         min = 0, step = 0.05, value = NULL)
          ),
          column(
            width = 3,
            style = 'margin-top: 18px;',
            checkboxInput("dePosOnly", "Only up-regulated genes",
                          value = FALSE)
          )
        ),
        fluidRow(
          column(
            width = 3,
            numericInput("deMinExp1", 
                         "Output Group1 mean expression greater than:",
                         min = 0, step = 0.1, value = NULL)
          ),
          column(
            width = 3,
            numericInput("deMaxExp2", 
                         "Output Group2 mean expression less than:",
                         min = 0, step = 0.1, value = NULL)
          ),
          column(
            width = 3,
            numericInput("deMinExpPerc1",
                         "Output Group1 expression percentage greater than:",
                         min = 0, max = 1, step = 0.05, value = NULL)
          ),
          column(
            width = 3,
            numericInput("deMaxExpPerc2",
                         "Output Group2 expression percentage less than:",
                         min = 0, max = 1, step = 0.05, value = NULL)
          )
        ),
        fluidRow(
          column(
            width = 3,
            textInput("deAnalysisName",
                      "Name of Differential Expression Analysis:",
                      placeholder = 'Required.')
          ),
          column(
            width = 2,
            style = 'margin-top: 25px;',
            withBusyIndicatorUI(actionButton("runDE", "Run"))
          )
        )
      )
    ),
    h3("Visualization"),
    p("For preview and result presentation.", style = "color:grey;"),
    fluidRow(
      selectInput("deResSel", "Select Differential Expression Analysis", choices = NULL),
      tabsetPanel(
        tabPanel(
          "Heatmap",
          panel(

            fluidRow(
              column(
                width = 4,
                dropdown(
                  fluidRow(
                    column(12,
                           fluidRow(actionBttn(inputId = "closeDropDownDeHM", label = NULL, style = "simple", color = "danger", icon = icon("times"), size = "xs"), align = "right"),
                    )
                  ),
                  fluidRow(
                    column(
                      width = 6,
                      checkboxInput('deHMDoLog', "Do log transformation", FALSE),
                    ),
                    column(
                      width = 6,
                      checkboxInput('deHMPosOnly', "Only up-regulated",
                                    value = FALSE)
                    )
                  ),
                  fluidRow(
                    column(
                      width = 6,
                      numericInput("deHMFC", "Aboslute log2FC value greater than:",
                                   value = 0.5, min = 0, step = 0.05),
                    ),
                    column(
                      width = 6,
                      numericInput("deHMFDR", "FDR value less than", value = 0.05,
                                   max = 1, step = 0.01),
                    )
                  ),
                  
                  fluidRow(
                    column(
                      width = 6,
                      numericInput("deHMMinExp1", 
                                   "Group1 mean expression greater than:",
                                   value = NULL, min = 0, step = 0.1),
                    ),
                    column(
                      width = 6,
                      numericInput("deHMMaxExp2", 
                                   "Group2 mean expression less than:", 
                                   value = NULL, min = 0, step = 0.1),
                    )
                  ),
                  fluidRow(
                    column(
                      width = 6,
                      numericInput("deHMMinExpPerc1", 
                                   "Group1 expression percentage greater than:",
                                   value = NULL, min = 0, max = 1, step = 0.05),
                    ),
                    column(
                      width = 6,
                      numericInput("deHMMaxExpPerc2", 
                                   "Group2 expression percentage less than:", 
                                   value = NULL, min = 0, max = 1, step = 0.05),
                    )
                  ),
                  
                  fluidRow(
                    column(
                      width = 6,
                      selectInput("deHMcolData", "Additional cell annotation",
                                  choices = clusterChoice, multiple = TRUE),
                    ),
                    column(
                      width = 6,
                      selectInput("deHMrowData", "Additional feature annotation",
                                  choices = featureChoice, multiple = TRUE),
                    )
                  ),
                  fluidRow(
                    column(
                      width = 6,
                      uiOutput('deHMSplitColUI'),
                    ),
                    column(
                      width = 6,
                      uiOutput('deHMSplitRowUI'),
                    )
                  ),
                  fluidRow(
                    column(
                      width = 6,
                      checkboxInput(
                        inputId = "deHMrowLabel",
                        label = "Add row labels",
                        value = FALSE  
                      )
                    ),
                    column(
                      width = 4,
                      withBusyIndicatorUI(
                        actionBttn(
                          inputId = "dePlotHM",
                          label = "Update",
                          style = "bordered",
                          color = "primary",
                          size = "sm"
                        )
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
              ),
              column(
                width = 7,
                fluidRow(
                  h6(
                    "A heatmap of the expression of DEGs in the selected cells, columns splitted by condition setting, and rows splitted by up-/down-regulation (whether log2FC is positive or negative, respectively)"),
                  align="center"
                )
              )
            ),
            hr(),
            br(),

            shinyjqui::jqui_resizable(plotOutput("deHeatmap"))
          )
        ),
        tabPanel("Results Table",
                 DT::dataTableOutput("deResult"),
                 downloadButton("deDownload", "Download Result Table")),
        tabPanel(
          "Violin Plot",
          panel(


            fluidRow(
              column(
                width = 4,
                dropdown(
                  fluidRow(
                    column(12,
                           fluidRow(actionBttn(inputId = "closeDropDownDeViolin", label = NULL, style = "simple", color = "danger", icon = icon("times"), size = "xs"), align = "right"),
                    )
                  ),
                  fluidRow(
                    div(style="display: inline-block;vertical-align:center; width: 100px;margin-left:10px",
                        p('Plot the top')),
                    div(style="display: inline-block;vertical-align:center; width: 60px;",
                        numericInput('deVioNRow', label = NULL, value = 4, min = 1)),
                    div(style="display: inline-block;vertical-align:center; width: 12px;",
                        p('x')),
                    div(style="display: inline-block;vertical-align:center; width: 60px;",
                        numericInput('deVioNCol', label = NULL, value = 4, min = 1)),
                    div(style="display: inline-block;vertical-align:center; width: 10px;",
                        p('=')),
                    div(style="display: inline-block;vertical-align:center; width: 30px;",
                        uiOutput('deVioTotalUI')),
                    div(style="display: inline-block;vertical-align:center; width: 50px;",
                        p('genes'))
                  ),
                  fluidRow(
                    column(
                      width = 6,
                      selectInput('deVioLabel', "Label features by",
                                  c("Default ID", featureChoice))
                    ),
                    column(
                      width = 4,
                      style = 'margin-top: 23px;',
                      withBusyIndicatorUI(
                        actionBttn(
                          inputId = "dePlotVio",
                          label = "Update",
                          style = "bordered",
                          color = "primary",
                          size = "sm"
                        )
                      )
                    )
                  ),
                  inputId = "dropDownDeViolin",
                  icon = icon("cog"),
                  status = "primary",
                  circle = FALSE,
                  inline = TRUE,
                  width = "500px"
                )
              ),
              column(
                width = 7,
                fluidRow(
                  h6(
                    "Violin plots of the expression of top DEGs in the selected analysis. The violin plot for each DEG will be grouped by the condition setting."),
                  align="center"
                )
              )
            ),
            hr(),
            br(),
            textOutput("deSanityWarnViolin"),
            shinyjqui::jqui_resizable(plotOutput("deViolinPlot"))
          ),

        ),
        tabPanel(
          "Linear Model",
          panel(

            fluidRow(
              column(
                width = 4,
                dropdown(
                  fluidRow(
                    column(12,
                           fluidRow(actionBttn(inputId = "closeDropDownDeReg", label = NULL, style = "simple", color = "danger", icon = icon("times"), size = "xs"), align = "right"),
                    )
                  ),
                  fluidRow(
                    div(style="display: inline-block;vertical-align:center; width: 100px;margin-left:10px",
                        p('Plot the top')),
                    div(style="display: inline-block;vertical-align:center; width: 60px;",
                        numericInput('deRegNRow', label = NULL, value = 4, min = 1)),
                    div(style="display: inline-block;vertical-align:center; width: 12px;",
                        p('x')),
                    div(style="display: inline-block;vertical-align:center; width: 60px;",
                        numericInput('deRegNCol', label = NULL, value = 4, min = 1)),
                    div(style="display: inline-block;vertical-align:center; width: 10px;",
                        p('=')),
                    div(style="display: inline-block;vertical-align:center; width: 30px;",
                        uiOutput('deRegTotalUI')),
                    div(style="display: inline-block;vertical-align:center; width: 50px;",
                        p('genes'))
                  ),
                  fluidRow(
                    column(
                      width = 6,
                      selectInput('deRegLabel', "Label features by",
                                  c("Default ID", featureChoice))
                    ),
                    column(
                      width = 4,
                      style = 'margin-top: 23px;',
                      withBusyIndicatorUI(
                        actionBttn(
                          inputId = "dePlotReg",
                          label = "Update",
                          style = "bordered",
                          color = "primary",
                          size = "sm"
                        )
                      )
                    )
                  ),
                  inputId = "dropDownDeReg",
                  icon = icon("cog"),
                  status = "primary",
                  circle = FALSE,
                  inline = TRUE,
                  width = "500px"
                )
              ),
              column(
                width = 7,
                fluidRow(
                  h6(
                    "Linear regression plots of the expression of top DEGs in the selected analysis. The regression plot for each DEG will be grouped by the condition setting."),
                  align="center"
                )
              )
            ),
            hr(),
            br(),
            textOutput("deSanityWarnReg"),
            shinyjqui::jqui_resizable(plotOutput("deRegPlot"))
          )
        )
      )
    )
  ),
  nonLinearWorkflowUI(id = "nlw-de")
)

