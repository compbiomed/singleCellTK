source("ui_02_qc.R", local = TRUE) #creates shinyPanelQC variable
source("ui_02_filter.R", local = TRUE) #creates shinyPanelFilter variable

shinyPanelQCFilter <- fluidPage(
  includeCSS('styles.css'),
  tabsetPanel(
    tabPanel("QC", shinyPanelQC),
    tabPanel("Filtering", shinyPanelFilter)
  )
)