shinyPanelRemove <- fluidPage(
  includeCSS('styles.CSS'),
  h1("Remove Data"),
  hr(),
  h5("Select data to remove:"),
  fluidRow(
    column(12,
                      uiOutput(
                        outputId = "assaysList"
                      ),
                      uiOutput(
                        outputId = "rowDataList"
                      ),
                      uiOutput(
                        outputId = "colDataList"
                      ),
                      uiOutput(
                        outputId = "reducedDimsList"
                      ),
                      uiOutput(
                        outputId = "altExpList"
                      ),
                      hr(),
                      uiOutput(outputId = "removeDataWarningUI"),
                      withBusyIndicatorUI(
                        actionButton(
                          inputId = "delRedDim", 
                          label = "Delete")
                        )
                    )
           ),
  br()
  )