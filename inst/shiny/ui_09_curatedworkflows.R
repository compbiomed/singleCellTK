source("ui_celda.R", local = TRUE)
shinyPanelCuratedWorkflows <- fluidPage(
  tabsetPanel(
    tabPanel("CELDA", shinyPanelCelda),
    tabPanel("Seurat"),
    tabPanel("Bioconductor/OSCA")
  )
)

