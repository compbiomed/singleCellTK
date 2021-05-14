shinyPanelEnrichR <- fluidPage(
  tags$div(
    class = "container",
    h1("Gene Set Enrichment Analysis using enrichR"),
    h5(tags$a(href = paste0(docs.artPath, "ui_enrichR.html"),
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        h3("Choose data source:"),
        radioButtons(
          "geneListChoice", label = NULL, c("Select Gene(s)" = "selectGenes",
                                            "Upload file" = "geneFile")
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
          hr(),
          h4("Sample files:"),
          tags$a(href = "https://drive.google.com/open?id=1iJZ6H_G2brbeww9B0dA5seMyYUZYyFrU",
                 "Gene Names", target = "_blank"
          ),
          tags$a(href = "https://drive.google.com/open?id=1BLrwW0uMi2pxsX0m1zJrlOTnBLkIiOhk",
                 "Entrez ids", target = "_blank"
          ),
          hr(),
          h4("Options:"),
          # Input: Checkbox if file has header ----
          checkboxInput("header", "Header", value = TRUE),
          # Input: Select separator ----
          radioButtons("sep", "Separator",
                       choices = c(Comma = ",",
                                   Semicolon = ";",
                                   Tab = "\t"),
                       selected = ",",
                       inline = TRUE),
          # Input: Select quotes ----
          radioButtons("quote", "Quote",
                       choices = c(None = "",
                                   "Double Quote" = '"',
                                   "Single Quote" = "'"),
                       selected = '""',
                       inline = TRUE)
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

