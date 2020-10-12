exampleDatasets <- c("mouseBrainSubset", "maits", "campBrainSubset")
if ("scRNAseq" %in% rownames(installed.packages())){
  exampleDatasets <- c(exampleDatasets, "fluidigm_pollen_et_al",
                       "th2_mahata_et_al", "allen_tasic_et_al")
}

shinyPanelUpload <- fluidPage(
  useShinyjs(),
  tags$style(appCSS),
  tags$div(
    class = "jumbotron", style = "background-color:#ededed",
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
  tags$br(),
  tags$div(
    class = "container",
    h1("Upload"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v03-tab01_Upload.html",
      "(help)", target = "_blank")),
    tags$hr(),
    tags$div(id = "uploadAlert", alertText),
    h3("Choose data source:"),
    radioButtons("uploadChoice", label = NULL, c("Upload files" = "files",
                                                 "Upload SummarizedExperiment/SCExperiment RDS File" = "rds",
                                                 "Use example data" = "example")
    ),
    tags$hr(),
    conditionalPanel(condition = sprintf("input['%s'] == 'files'", "uploadChoice"),
      h3("Upload data in tab separated text format:"),
      fluidRow(
        column(width = 4,
          wellPanel(
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
              "countsfile",
              HTML(
                paste("Input assay (eg. counts, required):",
                tags$span(style = "color:red", "*", sep = ""))
              ),
              accept = c(
                "text/csv", "text/comma-separated-values",
                "text/tab-separated-values", "text/plain", ".csv", ".tsv"
              )
            )
          ),
          h4("Input Assay Type:"),
          selectInput("inputAssayType", label = NULL,
                      c("counts", "normcounts", "logcounts", "cpm",
                        "logcpm", "tpm", "logtpm")
          ),
          checkboxInput("createLogcounts",
                        "Also create log2 input assay on upload", value = TRUE)
        ),
        column(width = 4,
          wellPanel(
            h4("Example sample annotation file:"),
            HTML('<table class="table"><thead><tr class="header"><th>Cell</th>
                 <th>Annot1</th><th>&#x2026;</th></tr></thead><tbody><tr class="odd">
                 <td>Cell1</td><td>a</td><td>&#x2026;</td></tr><tr class="even">
                 <td>Cell2</td><td>a</td><td>&#x2026;</td></tr><tr class="odd">
                 <td>Cell3</td><td>b</td><td>&#x2026;</td></tr><tr class="even">
                 <td>&#x2026;</td><td>&#x2026;</td><td>&#x2026;</td></tr><tr class="odd"><td>CellN</td>
                 <td>b</td><td>&#x2026;</td></tr></tbody></table>'),
            tags$a(href = "https://drive.google.com/open?id=10IDmZQUiASN4wnzO4-WRJQopKvxCNu6J",
                   "Download an example annotation file here.", target = "_blank"),
            tags$br(),
            tags$br(),
            fileInput(
              "annotFile", "Sample annotations (optional):",
              accept = c(
                "text/csv", "text/comma-separated-values",
                "text/tab-separated-values", "text/plain", ".csv", ".tsv"
              )
            )
          )
        ),
        column(width = 4,
          wellPanel(
            h4("Example feature file:"),
            HTML('<table class="table"><thead><tr class="header"><th>Gene</th>
               <th>Annot2</th><th>&#x2026;</th></tr></thead><tbody><tr class="odd">
                 <td>Gene1</td><td>a</td><td>&#x2026;</td></tr><tr class="even">
                 <td>Gene2</td><td>a</td><td>&#x2026;</td></tr><tr class="odd">
                 <td>Gene3</td><td>b</td><td>&#x2026;</td></tr><tr class="even">
                 <td>&#x2026;</td><td>&#x2026;</td><td>&#x2026;</td></tr><tr class="odd"><td>GeneM</td>
                 <td>b</td><td>&#x2026;</td></tr></tbody></table>'),
            tags$a(href = "https://drive.google.com/open?id=1gxXaZPq5Wrn2lNHacEVaCN2a_FHNvs4O",
                  "Download an example feature file here.", target = "_blank"),
            tags$br(),
            tags$br(),
            fileInput(
              "featureFile", "Feature annotations (optional):",
              accept = c(
                "text/csv", "text/comma-separated-values",
                "text/tab-separated-values", "text/plain", ".csv", ".tsv"
              )
            )
          )
        )
      )
    ),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'example'", "uploadChoice"),
      h3("Choose Example Dataset:"),
      selectInput("selectExampleData", label = NULL, exampleDatasets),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'mouseBrainSubset'", "selectExampleData"),
        h3(tags$a(href = "https://doi.org/10.1126/science.aaa1934", "Mouse Brain Subset: GSE60361", target = "_blank")),
        "A subset of 30 samples from a single cell RNA-Seq experiment from Zeisel, et al. Science 2015. The data was produced from cells from the mouse somatosensory cortex (S1) and hippocampus (CA1). 15 of the cells were identified as oligodendrocytes and 15 of the cell were identified as microglia.",
        tags$br(),
        tags$br()
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'campBrainSubset'", "selectExampleData"),
        h3(tags$a(href = "https://doi.org/10.1038/nn.4495", "500 cells from Campbell et. al, 2017, Mouse Brain Subset", target = "_blank")),
        "A subset of 500 cells from a single cell RNA-Seq experiment from Campbell, et al. Nature Neuroscience 2017 using droplet-based sequencing technology. This study was perfomed to identify various hypothalamic arcuate–median eminence complex (Arc-ME) cell types. This contains information such as the diet of the mice, sex and proposed cell type for each cell. ",
        tags$br(),
        tags$br()
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'maits'", "selectExampleData"),
        h3(tags$a(href = "https://doi.org/10.1186/s13059-015-0844-5", "MAITs data from MAST package", target = "_blank")),
        "96 Single-cell transcriptome profiling from Mucosal Associated Invariant T cells (MAITs), measured on the Fluidigm C1.",
        tags$br(),
        tags$br()
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'fluidigm_pollen_et_al'", "selectExampleData"),
        h3(tags$a(href = "http://dx.doi.org/10.1038/nbt.2967", "130 cells from (Pollen et al. 2014), 65 at high coverage and 65 at low coverage", target = "_blank")),
        "Transcriptomes of cell populations in both of low-coverage (~0.27 million reads per cell) and high-coverage (~5 million reads per cell) to identify cell-type-specific biomarkers, and to compare gene expression across samples specifically for cells of a given type as well as to reconstruct developmental lineages of related cell types. (data loaded from scRNASeq package)",
        tags$br(),
        tags$br()
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'th2_mahata_et_al'", "selectExampleData"),
        h3(tags$a(href = "http://dx.doi.org/10.1016/j.celrep.2014.04.011", "96 T helper cells from (Mahata et al. 2014)", target = "_blank")),
        "96 T helper cells from 6-week-old mouse, day 4.5 in vitro Th2 differentiation. (data loaded from scRNASeq package)",
        tags$br(),
        tags$br()
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'allen_tasic_et_al'", "selectExampleData"),
        h3(tags$a(href = "http://dx.doi.org/10.1038/nn.4216", "Mouse visual cortex cells from (Tasic et al. 2016)", target = "_blank")),
        "Subset of 379 cells from the mouse visual cortex. (data loaded from scRNASeq package)",
        tags$br(),
        tags$br()
      )
    ),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'rds'", "uploadChoice"),
      h3("Choose an RDS file that contains a SummarizedExperiment/SCExperiment Object:"),
      fileInput(
        "rdsFile", "SCExperiment RDS file:", accept = c(".rds", ".RDS")
      )
    ),
    withBusyIndicatorUI(
      actionButton("uploadData", "Upload")
    ),
    tags$div(
      class = "container",
      p("")
    )
  )
  #includeHTML("www/footer.html")
)

