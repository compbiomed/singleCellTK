source("ui_01_import.R", local = TRUE) #creates shinyPanelImport variable
source("ui_export.R", local = TRUE) #creates shinyPanelExport variable

shinyPanelData <- fluidPage(
  includeCSS('styles.css'),
  tabsetPanel(
    tabPanel("Import", shinyPanelImport),
    tabPanel("Export", shinyPanelExport)
  )
)