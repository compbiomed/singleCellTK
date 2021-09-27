shinyPanelRemove <- fluidPage(
  includeCSS('styles.CSS'),
  h1("Remove Data"),
  h5(tags$a(href = paste0(docs.artPath, "delete_data.html"),
            "(help)", target = "_blank")),
  hr(),
  h5("Select data to remove:"),
  fluidRow(
    column(
      12,
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
      h6("Warning: This action is inreversible. "),
      withBusyIndicatorUI(
        actionButton(
          inputId = "delRedDim",
          label = "Delete")
      )
    )
  ),
  br()
)
