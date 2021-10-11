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
      shinyDirButton("outputDirectory", label = "Select directory", title = "Download"),
      # A UI to display what users select
      verbatimTextOutput("outputDirectoryPath", placeholder = TRUE),

      tags$h5(style = "font-weight: bold; margin-bottom: 15px", "Choose export type"),
      radioButtons(
        "exportChoice",
        label = NULL,
        c(
          "Download as RDS file" = "rds",
          "Python annData object" = "annData",
          "Flat text files" = "textfile"
        )
      ),
      conditionalPanel(
        condition = "input.exportChoice == 'annData' || input.exportChoice == 'rds'",
        tags$label(id="exportFileNameLabel", "File Name"),
      ),
      conditionalPanel(
        condition = "input.exportChoice == 'textfile'",
        tags$label(id="exportFileNameLabel", "File prefix"),
      ),
      uiOutput("exportFileName"),
      actionButton("exportData", "Download")
    ),
    column(
      6,
      conditionalPanel(
        condition = "input.exportChoice == 'textfile'",
        tags$h5(style = "font-weight: bold; margin-bottom: 15px", "Set export specifications"),
        checkboxInput("exportFlatGzip", "Gzip Compress", value = TRUE)
      ),
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
        #uiOutput("exportAssay"),
        # tags$label(id="compressionLabel", "Compression"),
        # selectInput("compression", label = NULL, c("None", "lzf", "gzip"), width='140px'),
        tags$label(id="compressionOptsLabel", "Compression Opts"),
        numericInput("compressionOpts", label = NULL, 1, min = 1, max = 100, width='140px'),
        checkboxInput("forceDense", "Force Dense", value = FALSE)
      ),
      conditionalPanel(
        condition = "input.exportChoice == 'textfile' || input.exportChoice == 'annData'",
        checkboxInput("exportOverwrite", "Overwrite", value = TRUE)
      )
    )
    )
))

