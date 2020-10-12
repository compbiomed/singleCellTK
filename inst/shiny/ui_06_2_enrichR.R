shinyPanelEnrichR <- fluidPage(
  tags$div(
    class = "container",
    h1("Gene Set Enrichment Analysis using enrichR"),
    h5(tags$a(href = "", "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        h3("Choose data source:"),
        radioButtons(
          "geneListChoice", label = NULL, c("Select Gene(s)" = "selectGenes",
                                            "Upload file" = "geneFile",
                                            "Saved top genes" = "biomarker")
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'selectGenes'", "geneListChoice"),
          selectizeInput("enrichGenes", label = "Select Gene(s):", NULL, multiple = TRUE)
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'geneFile'", "geneListChoice"),
          fileInput(
            "enrFile",
            tags$b(tags$i("Please upload a file with only gene names or Entrez Gene Symbols")),
            accept = c("text/csv", "text/comma-separated-values",
                       "text/tab-separated-values", "text/plain",
                       ".csv", ".tsv")
          ),
          tags$b(tags$i("Sample files:")),
          br(),
          tags$a(href = "https://drive.google.com/open?id=1iJZ6H_G2brbeww9B0dA5seMyYUZYyFrU",
                 "Gene Names", target = "_blank"
          ),
          HTML('&nbsp;'),
          HTML('&nbsp;'),
          HTML('&nbsp;'),
          tags$a(href = "https://drive.google.com/open?id=1BLrwW0uMi2pxsX0m1zJrlOTnBLkIiOhk",
                 "Entrez ids", target = "_blank"
          ),
          br(),
          tags$div(
            tags$b(tags$i("Options:")),
            # Input: Checkbox if file has header ----
            fluidRow(
              column(width = 1,
                     checkboxInput("header", "Header", value = TRUE)
              ),
              # Input: Select separator ----
              column(width = 1,
                     offset = 2,
                     radioButtons("sep", "Separator",
                                  choices = c(Comma = ",",
                                              Semicolon = ";",
                                              Tab = "\t"),
                                  selected = ",")
              ),
              # Input: Select quotes ----
              column(width = 1,
                     offset = 3,
                     radioButtons("quote", "Quote",
                                  choices = c(None = "",
                                              "Double Quote" = '"',
                                              "Single Quote" = "'"),
                                  selected = '""')
              )
            )
          )
        ),
        conditionalPanel(
          helpText("To use this, first run Differential expression and save top genes."),
          condition = sprintf("input['%s'] == 'biomarker'", "geneListChoice"),
          uiOutput("enrBioGenes")
        ),
        selectizeInput("enrichDb", label = "Select DB:", c("ALL", enrichedDB),
                       multiple = TRUE),
        helpText("Selecting 'ALL' or leaving it blank will run enrichR on all available enrichR databases (N = 130) which will take significant amount of time."),
        withBusyIndicatorUI(actionButton("enrichRun", "Run")),
        br(),
        downloadButton("downloadEnrichR", "Download results")
      ),
      mainPanel(
        uiOutput("enrTabs")
      )
    )
  )
)

