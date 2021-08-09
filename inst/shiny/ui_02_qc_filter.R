source("ui_02_qc.R", local = TRUE) #creates shinyPanelQC variable
source("ui_02_filter.R", local = TRUE) #creates shinyPanelFilter variable

shinyPanelQCFilter <- fluidPage(
  includeCSS('styles.CSS'),
  h1("Data QC & Filtering"),
  tabsetPanel(id = "QCFilterTabsetPanel",
    tabPanel("QC", shinyPanelQC),
    tabPanel("Filtering", shinyPanelFilter)
  ),
  nonLinearWorkflowUI(id = "nlw-qcf")
)

