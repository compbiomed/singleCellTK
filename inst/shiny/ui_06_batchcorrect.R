shiny_panel_batchcorrect <- fluidPage(
  tags$div(
    class = "container",
    h1("Batch Correction"),
    sidebarLayout(
      sidebarPanel(
        selectInput("combatAssay", "Select Assay:", ""),
        selectInput("combatBatchVar", "Select Batch Condition:", ""),
        selectInput("combatConditionVar", "Select Additional Covariates:", "", multiple = TRUE),
        radioButtons("combatParametric", "Adjustments:", c("Parametric", "Non-parametric"), selected = "Parametric"),
        checkboxInput("combatMeanOnly", "Correct mean of the batch effect only", value = FALSE),
        checkboxInput("combatRef", "Run reference batch combat:", value = FALSE),
        conditionalPanel(
          condition = "input.combatRef == true",
          selectInput("combatRefBatch", "Choose Reference Batch:", "")
        ),
        textInput("combatSaveAssay", "Assay Name to Use:", value = "combat"),
        withBusyIndicatorUI(actionButton("combatRun", "Run"))
      ),
      mainPanel(
        "Main panel"
      )
    )
  )
)
