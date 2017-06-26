exampleDatasets <- c("GSE36552", "GSE60361_subset", "GSE66507", "GSE73121",
                     "maits_SCESet")

shiny_panel_upload <- fluidPage(
  useShinyjs(),
  tags$style(appCSS),
  tags$div(
    class = "jumbotron",
    tags$div(
      class = "container",
      h1("Single Cell Toolkit"),
      p("Filter, cluster, and analyze single cell RNA-Seq data")
    )
  ),
  tags$div(
    class = "container",
    tags$div(id = "uploadAlert", alertText),
    radioButtons("uploadChoice", "Upload:",
                 c("Files" = "files",
                   "Example data" = "example")),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'files'", "uploadChoice"),
      fileInput("countsfile", "Upload a matrix of counts here",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values",
                  "text/tab-separated-values",
                  "text/plain",
                  ".csv",
                  ".tsv"
                )
      ),
      fileInput("annotfile", "Optional: Upload a matrix of annotations here",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values",
                  "text/tab-separated-values",
                  "text/plain",
                  ".csv",
                  ".tsv"
                )
      ),
      fileInput("featurefile", "Optional: Upload a matrix of feature annotations here",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values",
                  "text/tab-separated-values",
                  "text/plain",
                  ".csv",
                  ".tsv"
                )
      )
    ),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'example'", "uploadChoice"),
      selectInput("selectExampleData", "Or, choose example data:",
                  exampleDatasets)
    ),
    withBusyIndicatorUI(
      actionButton("uploadData", "Upload")
    ),
    tags$div(
      class = "container",
      p("")
    )
  ),
  includeHTML("www/footer.html")
)
