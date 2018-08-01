library(shinycssloaders)
shinyPanelEnrichR <- fluidPage(
  tags$div(
    class = "container",
    h1("Gene Enrichment Analysis using enrichR"),
    h5(tags$a(href = "", "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        selectInput("enrichAssay", "Select Assay", currassays),
        h3("Choose data source:"),
        radioButtons("geneListChoice", label = NULL, c("Select Gene(s)" = "selectGenes",
                                                       "Upload file" = "geneFile",
                                                       "Biomarker" = "biomarker")
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'selectGenes'", "geneListChoice"),
          selectizeInput("enrichGenes", label = "Select Gene(s):", NULL, multiple = TRUE)
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'geneFile'", "geneListChoice"),
          fileInput("enrFile",
                    tags$b(tags$i("Please upload a file with only gene names(Entrez Gene Symbols)")),
                    accept = c("text/csv", "text/comma-separated-values",
                               "text/tab-separated-values", "text/plain", ".csv", ".tsv")
          ),
          tags$b(tags$i("Options:")),
          # Input: Checkbox if file has header ----
          checkboxInput("header", "Header", TRUE),
          # Input: Select separator ----
          radioButtons("sep", "Separator",
                       choices = c(Comma = ",",
                                   Semicolon = ";",
                                   Tab = "\t"),
                       selected = ","),
          # Input: Select quotes ----
          radioButtons("quote", "Quote",
                       choices = c(None = "",
                                   "Double Quote" = '"',
                                   "Single Quote" = "'"),
                       selected = '""')
        ),
        # conditionalPanel(
        #   condition = sprintf("input['%s'] == 'biomarker'", "geneListChoice")
        # ),
        selectizeInput("enrichDb", label = "Select DB:", c("ALL", enrichedDB), 
                       multiple = TRUE),
        withBusyIndicatorUI(actionButton("enrichRun", "Run")),
        br(),
        downloadButton("downloadEnrichR", "Download results")
      ),
      mainPanel(
        uiOutput("enrTabs") %>% shinycssloaders::withSpinner()
      )
    )
  )
)
