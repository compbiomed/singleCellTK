shiny_panel_batchcorrect <- fluidPage(
  tags$div(
    class = "container",
    h1("Batch Correction"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v06-tab04_Batch-Correction.html",
              "(help)", target = "_blank")),
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
        uiOutput("selectCombat_refbatchUI"),
        textInput("combatSaveAssay", "Assay Name to Use:", value = "combat"),
        withBusyIndicatorUI(actionButton("combatRun", "Run"))
      ),
      mainPanel(
        uiOutput("combatStatus"),
        plotOutput("combatBoxplot", height = "600px")
      )
    )
  )
)
