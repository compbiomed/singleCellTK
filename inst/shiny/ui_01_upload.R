shiny_panel_upload <- fluidPage(
  useShinyjs(),
  tags$style(appCSS),
  tags$div(
    class="jumbotron",
    tags$div(
      class="container",
      h1("Single Cell Toolkit"),
      p("Filter, cluster, and analyze single cell RNA-Seq data")
    )
  ),
  tags$div(
    class="container",
    tags$div(id="uploadAlert", alertText),
    fileInput('countsfile', 'Upload a matrix of counts here',
              accept = c(
                'text/csv',
                'text/comma-separated-values',
                'text/tab-separated-values',
                'text/plain',
                '.csv',
                '.tsv'
              )
    ),
    fileInput('annotfile', 'Optional: Upload a matrix of annotations here',
              accept = c(
                'text/csv',
                'text/comma-separated-values',
                'text/tab-separated-values',
                'text/plain',
                '.csv',
                '.tsv'
              )
    ),
    fileInput('featurefile', 'Optional: Upload a matrix of feature annotations here',
              accept = c(
                'text/csv',
                'text/comma-separated-values',
                'text/tab-separated-values',
                'text/plain',
                '.csv',
                '.tsv'
              )
    ),
    withBusyIndicatorUI(
      actionButton("uploadData", "Upload")
    ),
    tags$div(
      class="container",
      p("")
    )
  ),
  includeHTML('www/footer.html')
)
