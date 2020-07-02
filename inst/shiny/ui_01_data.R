source("ui_01_import.R", local = TRUE) #creates shinyPanelImport variable
source("ui_01_columnAnnotation.R", local = TRUE) #creates shinyPanelColumnAnnotation variable
source("ui_01_rowAnnotation.R", local = TRUE) #creates shinyPanelRowAnnotation variable
source("ui_export.R", local = TRUE) #creates shinyPanelExport variable

shinyPanelData <- fluidPage(
  includeCSS('styles.css'),
  tabsetPanel(
    tabPanel("Import", shinyPanelImport),
    tabPanel("Column Annotation", shinyPanelColumnAnnotation),
    tabPanel("Row Annotation", shinyPanelRowAnnotation),
    tabPanel("Export", shinyPanelExport)
  )
)