shinyPanelRemove <- fluidPage(
  includeCSS('styles.CSS'),
  
  fluidRow(
    column(12,
           fluidRow(
             column(12,
                    panel(
                      heading = "Remove Data",
                      h6("Select data to remove:"),
                      uiOutput(
                        outputId = "assaysList"
                      ),
                      uiOutput(
                        outputId = "reducedDimsList"
                      ),
                      hr(),
                      uiOutput(outputId = "removeDataWarningUI"),
                      withBusyIndicatorUI(
                        actionButton(
                          inputId = "delRedDim", 
                          label = "Delete")
                        )
                      )
                    )
           )
           )
  )
)