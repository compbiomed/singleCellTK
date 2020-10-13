shinyPanelHeatmap <- fluidPage(
  h1("Heatmap"),
  p("Generic heatmap plotting panel for customized figure.",
    style = "color:grey;"),
  panel(
    fluidRow(
      div(style="display:inline-block;vertical-align:top;width:120px;margin-left:20px;",
          h4("Assay to Plot")),
      div(style="display:inline-block;vertical-align:bottom;width:120px;margin-top:5px;",
          selectInput("hmAssay", NULL, currassays))
    ),
    fluidRow(
      div(style="display:inline-block;vertical-align:top;width:180px;margin-left:20px;",
          h4("Import from analysis")),
      div(style="display:inline-block;vertical-align:bottom;width:150px;margin-top:5px;",
          selectInput("hmImport", NULL,
                      c("None",
                        "Differential Expression",
                        "Find Marker"),
                      selected = "None")),
      div(style="display:inline-block;vertical-align:bottom;width:50px;margin-left:8px;margin-bottom:15px;",
          actionButton("hmImportRun", "Import"))
    ),
    conditionalPanel(
      condition = "input.hmImport == 'Differential Expression'",
      uiOutput("hmImpDEGUI")
    ),
    hr(),
    # Subset ####
    h3("Cell/Feature Subsetting"),
    p("Only to plot cells/features of interests", style = "color:grey;"),
    tabsetPanel(id = 'hmSubsetTSP',
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
      tabPanel(
        title = "Manual Enter", value = 'hmSubsetManualTP',
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
    uiOutput("hmCellSumUI"),
    uiOutput("hmGeneSumUI"),
    hr(),
    # Annotation ####
    h3("Annotation Setting"),
    p("Stick additional information at sides of the plot",
      style = "color:grey;"),
    tabsetPanel(id = 'hmAnnTSP',
      tabPanel(
        title = "Cell", value = 'hmAnnCellTP',
        uiOutput('hmCellAnnUI'),
        uiOutput('hmCellAnnAssUI')
      ),
      tabPanel(
        title = "Feature", value = 'hmAnnGeneTP',
        uiOutput('hmGeneAnnUI'),
        uiOutput('hmGeneAnnAssUI')
      )
    ),
    hr(),
    # Others ####
    h3("Heatmap Setting"),
    p("Settings for split, label, dendrogram, color scheme and etc.",
      style = "color:grey;"),
    actionLink("hmHideDiv3", "Show/Hide"),
    panel(
      fluidRow(
        column(
          width = 6,
          uiOutput('hmColSplitUI'),
          checkboxGroupInput('hmAddLabel', "Add cell/feature labels",
                             choiceNames = c('cells', 'features'),
                             inline = TRUE, choiceValues = c(1, 2)),
          conditionalPanel(
            condition = "input.hmAddLabel.includes('1')",
            selectInput('hmAddCellLabel', "Add cell labels from",
                        c("Default cell IDs", clusterChoice),
                        selected = "Default cell IDs")
          ),
          conditionalPanel(
            condition = "input.hmAddLabel.includes('2')",
            selectInput('hmAddGeneLabel', "Add feature labels from",
                        c("Default feature IDs", featureChoice),
                        selected = "Default feature IDs")
          ),
          uiOutput('hmTrimUI')
        ),
        column(
          width = 6,
          uiOutput('hmRowSplitUI'),
          checkboxGroupInput('hmShowDendro', "Show dendrograms for",
                             choiceNames = c('cells', 'features'),
                             inline = TRUE, choiceValues = c(1, 2),
                             selected = c(1, 2)),
          checkboxInput("hmScale", "Z-Score Scale", value = TRUE)
        )
      ),
      h4("Color Scheme"),
      selectizeInput(
        inputId = 'hmCSPalette',
        label = "Choose from preset:",
        choices = NULL,
      ),
      fluidRow(
        column(
          width = 4,
          colourpicker::colourInput('hmCSLow', 'Low color',value = 'blue')
        ),
        column(
          width = 4,
          colourpicker::colourInput('hmCSMedium', 'Medium color',value = 'white')
        ),
        column(
          width = 4,
          colourpicker::colourInput('hmCSHigh', 'High color',value = 'red')
        )
      )
      # TODO: Do we add save preset button?
    ),
    hr(),
    withBusyIndicatorUI(actionButton("plotHeatmap", "Plot Heatmap")),
    div(
      style = 'height:800px;',
      plotOutput("Heatmap")
    )
    #plotOutput("Heatmap")
  )
)

