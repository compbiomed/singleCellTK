shinyPanelExport <- fluidPage(
  tags$div(
  class = "container",
  style = "margin-bottom: 10px",
  h1("Export Data"),
  tags$hr(),
  fluidRow(
    column(
      6,
      shinyDirectoryInput::directoryInput('outputDirectory', label = 'Select directory', value = '~'),
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
      actionButton("exportData", "Download")
    ),
    column(
      6,
      conditionalPanel(
        condition = "input.exportChoice === 'textfile'",
        tags$h5(style = "font-weight: bold; margin-bottom: 15px", "Set export specifications"),
        tags$label(id="gzipLabel", "Gzip"),
        selectInput("gzip", label=NULL, c("True", "False"), width = '140px')
      ),
      conditionalPanel(
        condition = "input.exportChoice === 'annData'",
        tags$h5(style = "font-weight: bold; margin-bottom: 15px", "Set export specifications"),
        tags$label(id="exportAssayLabel", "Assay"),
        selectInput("exportAssay", label=NULL, c(""), width='140px'),
        tags$label(id="compressionLabel", "Compression"),
        selectInput("compression", label = NULL, c("None", "lzf", "gzip"), width='140px'),
        tags$label(id="compressionOptsLabel", "Compression Opts"),
        numericInput("compressionOpts", label = NULL, 0, min = 1, max = 100, width='140px'),
        tags$label(id="forceDenseLabel", "Force Dense"),
        selectInput("forceDense", label = NULL, c("False", "True"), width='140px'),
      ),
      conditionalPanel(
        condition = "input.exportChoice === 'textfile' || input.exportChoice === 'annData'",
        tags$label(id="overwriteLabel", "Overwrite"),
        selectInput("overwrite", label = NULL,c("True", "False"), width = '140px'),
      )
    )
    )
))

