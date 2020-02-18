shinyPanelBatchcorrect <- fluidPage(
  includeCSS('styles.css'),
  tabsetPanel(
  tabPanel(
    "Normalization", fluid = TRUE,
    sidebarLayout(
      sidebarPanel(
        h4("Assay Options:"),
        selectInput(
          "assayModifyAction",
          "Assay Actions:",
          c(
            "Log Transform" = "log",
            "Create CPM" = "cpm",
            "Rename" = "rename",
            "Delete" = "delete"
          )
        ),
        conditionalPanel(condition = "input.assayModifyAction == 'cpm'",
                         h5(
                           "Select a count assay to use for CPM calculation:"
                         )),
        selectInput("modifyAssaySelect", "Select Assay:", currassays),
        conditionalPanel(
          condition = "input.assayModifyAction != 'delete'",
          textInput("modifyAssayOutname", "Assay Name", "",
                    placeholder = "What should the assay be called?")
        ),
        withBusyIndicatorUI(actionButton("modifyAssay", "Run"))
      ),
      mainPanel(fluidRow(column(
        12,
        h4("Available Assays:"),
        tableOutput("assayList")
      )))
    )),
  tabPanel(
    "Batch Correction", 
    sidebarLayout(
      sidebarPanel(
        selectInput("combatAssay", "Select Assay:", currassays),
        tags$hr(),
        h4("Plot Batch Effect:"),
        selectInput("batchVarPlot", "Select Batch Annotation:", c("none", clusterChoice)),
        selectInput("conditionVarPlot", "Select Condition Annotation:", c("none", clusterChoice)),
        tags$hr(),
        h4("Run Batch Correction:"),
        selectInput("batchMethod", "Select Method:", "ComBat"),
        selectInput("combatBatchVar", "Select Batch Condition:", clusterChoice),
        selectInput("combatConditionVar", "Select Additional Covariates:",
                    clusterChoice, multiple = TRUE),
        radioButtons("combatParametric", "Adjustments:", c("Parametric",
                                                           "Non-parametric"),
                     selected = "Parametric"),
        checkboxInput("combatMeanOnly", "Correct mean of the batch effect only",
                      value = FALSE),
        checkboxInput("combatRef", "Run reference batch combat:",
                      value = FALSE),
        uiOutput("selectCombatRefBatchUI"),
        textInput("combatSaveAssay", "Assay Name to Use:", value = "combat"),
        withBusyIndicatorUI(actionButton("combatRun", "Run"))
      ),
      mainPanel(
        uiOutput("combatStatus"),
        plotOutput("combatBoxplot", height = "600px")
      )
    )
  )
))