exampleDatasets <- c("mouse_brain_subset", "maits")

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
    conditionalPanel(condition = sprintf("input['%s'] == 'files'", "uploadChoice"),
      h3("Upload data in tab separated text format:"),
      fluidRow(
        column(width = 4,
          fileInput("countsfile", "Input assay (eg. counts, required):",
                   accept = c(
                     "text/csv",
                     "text/comma-separated-values",
                     "text/tab-separated-values",
                     "text/plain",
                     ".csv",
                     ".tsv"
                   )
          ),
          selectInput("inputAssayType", "Input Assay Type:", c("counts",
                                                               "normcounts",
                                                               "logcounts",
                                                               "cpm", "logcpm",
                                                               "tpm", "logtpm")),
          checkboxInput("createLogcounts", "Also create log2 input assay on upload", value = TRUE)
        ),
        column(width = 4,
          fileInput("annotfile", "Sample annotations (optional):",
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
        column(width = 4,
          fileInput("featurefile", "Feature annotations (optional):",
                   accept = c(
                     "text/csv",
                     "text/comma-separated-values",
                     "text/tab-separated-values",
                     "text/plain",
                     ".csv",
                     ".tsv"
                   )
          )
        )
      )
    ),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'example'", "uploadChoice"),
      selectInput("selectExampleData", "Or, choose example data:",
                  exampleDatasets),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'mouse_brain_subset'", "selectExampleData"),
        h3(tags$a(href="https://doi.org/10.1126/science.aaa1934", "Mouse Brain Subset: GSE60361")),
        "A subset of 30 samples from a single cell RNA-Seq experiment from Zeisel, et al. Science 2015. The data was produced from cells from the mouse somatosensory cortex (S1) and hippocampus (CA1). 15 of the cells were identified as oligodendrocytes and 15 of the cell were identified as microglia.",
        tags$br(),
        tags$br()
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'maits'", "selectExampleData"),
        h3(tags$a(href="https://doi.org/10.1186/s13059-015-0844-5", "MAITs data from MAST package")),
        "96 Single-cell transcriptome profiling from Mucosal Associated Invariant T cells (MAITs), measured on the Fluidigm C1.",
        tags$br(),
        tags$br()
      )
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
