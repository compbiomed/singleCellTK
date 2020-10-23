source("ui_02_qc.R", local = TRUE) #creates shinyPanelQC variable
source("ui_02_filter.R", local = TRUE) #creates shinyPanelFilter variable

shinyPanelQCFilter <- fluidPage(
  includeCSS('styles.CSS'),
  h1("Data QC & Filtering"),
  h5(tags$a(href = "https://www.sctk.science/articles/tab02_data-summary-and-filtering",
            "(help)", target = "_blank")),
  tabsetPanel(
    tabPanel("QC", shinyPanelQC),
    tabPanel("Filtering", shinyPanelFilter)
  )
)

