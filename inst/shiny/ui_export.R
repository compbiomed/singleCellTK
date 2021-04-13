shinyPanelExport <- fluidPage(
  tags$div(
  class = "container",
  style = "margin-bottom: 10px",
  h1("Export Data"),
  tags$hr(),
  fluidRow(
    column(
      6,
      #shinyDirectoryInput::directoryInput('outputDirectory', label = 'Select directory', value = '~'),

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
      tags$label(id="exportFileNameLabel", "File Name"),
      uiOutput("exportFileName"),
      actionButton("exportData", "Download")
    ),
    column(
      6,
      conditionalPanel(
        condition = "input.exportChoice === 'textfile'",
        tags$h5(style = "font-weight: bold; margin-bottom: 15px", "Set export specifications"),
        checkboxInput("exportFlatGzip", "Gzip Compress", value = TRUE)
      ),
      conditionalPanel(
        condition = "input.exportChoice === 'annData'",
        tags$h5(style = "font-weight: bold; margin-bottom: 15px", "Set export specifications"),
        tags$label(id="exportAssayLabel", "Assay"),
        uiOutput("exportAssay"),
        tags$label(id="compressionLabel", "Compression"),
        # selectInput("compression", label = NULL, c("None", "lzf", "gzip"), width='140px'),
        tags$label(id="compressionOptsLabel", "Compression Opts"),
        numericInput("compressionOpts", label = NULL, 1, min = 1, max = 100, width='140px'),
        checkboxInput("forceDense", "Force Dense", value = FALSE)
      ),
      conditionalPanel(
        condition = "input.exportChoice === 'textfile' || input.exportChoice === 'annData'",
        checkboxInput("exportOverwrite", "Overwrite", value = TRUE)
      )
    )
    )
))

