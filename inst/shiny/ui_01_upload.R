exampleDatasets <- c("mouse_brain_subset", "maits")

shiny_panel_upload <- fluidPage(
  useShinyjs(),
  tags$style(appCSS),
  tags$div(
    class = "jumbotron",
    tags$div(
      class = "container",
      h1("Single Cell Toolkit"),
      p("Filter, cluster, and analyze single cell RNA-Seq data"),
      p(
        "Need help?",
        tags$a(href = "https://compbiomed.github.io/sctk_docs/",
               "Read the docs.", target = "_blank")
      )
    )
  ),
  tags$div(
    class = "container",
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v03-tab01_Upload.html",
              "(Upload tab help)", target = "_blank")),
    tags$div(id = "uploadAlert", alertText),
    radioButtons("uploadChoice", "Upload:",
                 c("Files" = "files",
                   "Example data" = "example")),
    conditionalPanel(condition = sprintf("input['%s'] == 'files'", "uploadChoice"),
      h3("Upload data in tab separated text format:"),
      fluidRow(
        column(width = 4,
          h4("Example count file:"),
          HTML('<table class="table"><thead><tr class="header"><th>Gene</th>
               <th>Cell1</th><th>Cell2</th><th>&#x2026;</th><th>CellN</th>
               </tr></thead><tbody><tr class="odd"><td>Gene1</td><td>0</td>
               <td>0</td><td>&#x2026;</td><td>0</td></tr><tr class="even">
               <td>Gene2</td><td>5</td><td>6</td><td>&#x2026;</td><td>0</td>
               </tr><tr class="odd"><td>Gene3</td><td>4</td><td>3</td>
               <td>&#x2026;</td><td>8</td></tr><tr class="even">
               <td>&#x2026;</td><td>&#x2026;</td><td>&#x2026;</td>
               <td>&#x2026;</td><td>&#x2026;</td></tr><tr class="odd">
               <td>GeneM</td><td>10</td><td>10</td><td>&#x2026;</td><td>10</td>
               </tr></tbody></table>'),
          tags$a(href = "https://drive.google.com/open?id=1n0CtM6phfkWX0O6xRtgPPg6QuPFP6pY8",
                 "Download an example count file here.", target = "_blank"),
          tags$br(),
          tags$br(),
          fileInput(
            "countsfile", "Input assay (eg. counts, required):",
            accept = c(
              "text/csv", "text/comma-separated-values",
              "text/tab-separated-values", "text/plain", ".csv", ".tsv"
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
          h4("Example sample annotation file:"),
          HTML('<table class="table"><thead><tr class="header"><th>Cell</th>
               <th>Annot1</th><th>…</th></tr></thead><tbody><tr class="odd">
               <td>Cell1</td><td>a</td><td>…</td></tr><tr class="even">
               <td>Cell2</td><td>a</td><td>…</td></tr><tr class="odd">
               <td>Cell3</td><td>b</td><td>…</td></tr><tr class="even">
               <td>…</td><td>…</td><td>…</td></tr><tr class="odd"><td>CellN</td>
               <td>b</td><td>…</td></tr></tbody></table>'),
          tags$a(href = "https://drive.google.com/open?id=10IDmZQUiASN4wnzO4-WRJQopKvxCNu6J",
                 "Download an example annotation file here.", target = "_blank"),
          tags$br(),
          tags$br(),
          fileInput(
            "annotfile", "Sample annotations (optional):",
            accept = c(
              "text/csv", "text/comma-separated-values",
              "text/tab-separated-values", "text/plain", ".csv", ".tsv"
            )
          )
        ),
        column(width = 4,
          h4("Example feature file:"),
          HTML('<table class="table"><thead><tr class="header"><th>Gene</th>
               <th>Annot2</th><th>…</th></tr></thead><tbody><tr class="odd">
               <td>Gene1</td><td>a</td><td>…</td></tr><tr class="even">
               <td>Gene2</td><td>a</td><td>…</td></tr><tr class="odd">
               <td>Gene3</td><td>b</td><td>…</td></tr><tr class="even">
               <td>…</td><td>…</td><td>…</td></tr><tr class="odd"><td>GeneM</td>
               <td>b</td><td>…</td></tr></tbody></table>'),
          tags$a(href = "https://drive.google.com/open?id=1gxXaZPq5Wrn2lNHacEVaCN2a_FHNvs4O",
                "Download an example feature file here.", target = "_blank"),
          tags$br(),
          tags$br(),
          fileInput(
            "featurefile", "Feature annotations (optional):",
            accept = c(
              "text/csv", "text/comma-separated-values",
              "text/tab-separated-values", "text/plain", ".csv", ".tsv"
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
        h3(tags$a(href = "https://doi.org/10.1126/science.aaa1934", "Mouse Brain Subset: GSE60361", target = "_blank")),
        "A subset of 30 samples from a single cell RNA-Seq experiment from Zeisel, et al. Science 2015. The data was produced from cells from the mouse somatosensory cortex (S1) and hippocampus (CA1). 15 of the cells were identified as oligodendrocytes and 15 of the cell were identified as microglia.",
        tags$br(),
        tags$br()
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'maits'", "selectExampleData"),
        h3(tags$a(href = "https://doi.org/10.1186/s13059-015-0844-5", "MAITs data from MAST package", target = "_blank")),
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
