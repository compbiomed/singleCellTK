source("ui_01_import.R", local = TRUE) #creates shinyPanelImport variable

shinyPanelData <- fluidPage(
  includeCSS('styles.css'),
  tabsetPanel(
    tabPanel("Import", shinyPanelImport)
  )
)