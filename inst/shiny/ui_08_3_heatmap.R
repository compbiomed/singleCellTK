shinyPanelHeatmap <- fluidPage(
  h1("Heatmap"),
  p("Generic heatmap plotting panel for customized figure.",
    style = "color:grey;"),
  useShinyjs(),
  panel(
    fluidRow(
      div(style="display:inline-block;vertical-align:top;width:120px;margin-left:20px;",
          h4("Assay to Plot")),
      div(style="display:inline-block;vertical-align:bottom;width:120px;margin-top:5px;",
          selectInput("hmAssay", NULL, currassays))
    ),
    fluidRow(
      div(style="display:inline-block;vertical-align:top;width:230px;margin-left:20px;",
          h4("Import from precalculation")),
      div(style="display:inline-block;vertical-align:bottom;width:150px;margin-top:5px;",
          selectInput("hmImport", NULL,
                      c("None", "MAST DEG", "MAST Marker"), selected = "None")),
      div(style="display:inline-block;vertical-align:bottom;width:50px;margin-left:8px;margin-bottom:15px;",
          actionButton("hmImportRun", "Import"))
    ),
    # Subset ####
    h3("Cell/Feature Subsetting"),
    p("Only to plot cells/features of interests", style = "color:grey;"),
    actionLink('hmHideDiv1', "Show/Hide"),
    hidden(
      div(
        id='hmDiv1',
        tabsetPanel(
          tabPanel(
            "Cell",
            tagList(
              selectInput(
                'hmCellCol',
                "Columns to display (Changes here clear the selection)",
                clusterChoice, multiple = TRUE, width = '550px'),
              DT::dataTableOutput("hmCellColTable"),
              actionButton('hmCellColTable_addAll', "Add all filtered"),
              actionButton('hmCellColTable_clear', "Clear selection"),
            )
          ),
          tabPanel(
            "Feature",
            tagList(
              selectInput(
                'hmGeneCol',
                "Columns to display (Changes here clear the selection)",
                featureChoice, multiple = TRUE, width = '550px'),
              DT::dataTableOutput("hmGeneColTable"),
              actionButton('hmGeneColTable_addAll', "Add all filtered"),
              actionButton('hmGeneColTable_clear', "Clear selection"),
            )
          ),
          tabPanel(
            "Manual Enter",
            fluidRow(
              column(
                width = 6,
                selectInput("hmCellTextBy", "Match input cell identifiers by:",
                            c("Row Names", clusterChoice), selected = "Row Names"),
                checkboxInput('hmCellTextEM', "Exact Match", value = TRUE),
                checkboxInput('hmCellTextFM', "First Match", value = TRUE),
                textAreaInput(
                  inputId = 'hmCellText',
                  label = "Enter Cell Identifiers",
                  placeholder = "One cell per line, with no symbol separator."
                ),
                uiOutput('hmCellNEnteredUI'),
                actionButton("hmCellAddFromText", "Add")
              ),
              column(
                width = 6,
                selectInput("hmGeneTextBy", "Match input feature identifiers by:",
                            c("Row Names", featureChoice), selected = "Row Names"),
                checkboxInput('hmGeneTextEM', "Exact Match", value = TRUE),
                checkboxInput('hmGeneTextFM', "First Match", value = TRUE),
                textAreaInput(
                  inputId = 'hmGeneText',
                  label = "Enter Feature Identifiers",
                  placeholder = "One Feature per line, with no symbol separator."
                ),
                uiOutput('hmGeneNEnteredUI'),
                actionButton("hmGeneAddFromText", "Add")
              )
            )
          )
        ),
      )
    ),
    uiOutput("hmCellSumUI"),
    uiOutput("hmGeneSumUI"),
    hr(),
    # Annotaion ####
    h3("Annotation Setting"),
    p("Stick additional information at sides of the plot",
      style = "color:grey;"),
    actionLink('hmHideDiv2', "Show/Hide"),
    hidden(
      div(
        id = 'hmDiv2',
        panel(
          "Add annotations"
        ),
      )
    ),
    hr(),
    # Others ####
    h3("Heatmap Setting"),
    p("Settings for split, label, dendrogram, colormap, legend and etc.",
      style = "color:grey;"),
    actionLink("hmHideDiv3", "Show/Hide"),
    hidden(
      div(
        id = 'hmDiv3',
        panel(
          p("Set split, label, dendrogram, colormap, legend and etc.")
        ),
      )
    ),
    hr(),
    withBusyIndicatorUI(actionButton("plotHeatmap", "Plot Heatmap")),
    panel(
      plotOutput("Heatmap", width="600px", height="600px")
    )
  )
)
