shinyPanelExport <- fluidPage(
  tags$div(
  class = "container",
  style = "margin-bottom: 10px",
  h1("Export Data"),
  h5(tags$a(href = paste0(docs.artPath, "export_data.html"),
            "(help)", target = "_blank")),
  tags$hr(),
  fluidRow(
    column(
      6,
      tags$h5(style = "font-weight: bold; margin-bottom: 15px", "Choose export type"),
      radioButtons(
        "exportChoice",
        label = NULL,
        c(
          "Download as RDS file" = "rds",
          "Python annData object" = "annData"
        )
      ),
      conditionalPanel(
        condition = "input.exportChoice == 'annData' || input.exportChoice == 'rds'",
        tags$label(id="exportFileNameLabel", "File Name"),
      ),
      
      uiOutput("exportFileName"),
      downloadButton("exportData", "Download")
    ),
    column(
      6,
      conditionalPanel(
        condition = "input.exportChoice == 'annData'",
        tags$h5(style = "font-weight: bold; margin-bottom: 15px", "Set export specifications"),
        selectizeInput(
          inputId = "exportAssay", 
          label = "Select input matrix:", 
          choices = NULL, 
          selected = NULL, 
          multiple = FALSE,
          options = NULL),
        
      ),
    )
    )
))

