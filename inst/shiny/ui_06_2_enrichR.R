shinyPanelEnrichR <- fluidPage(
  tags$div(
    class = "container",
    h1("Gene Set Enrichment Analysis using enrichR"),
    h5(tags$a(href = paste0(docs.artPath, "cnsl_enrichR.html"),
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        h4("Choose data source:"),
        radioButtons(
          "geneListChoice", label = NULL, c("Import DEG" = "deg",
                                            "Select Gene(s)" = "selectGenes",
                                            "Upload file" = "geneFile")
        ),
        hr(),
        conditionalPanel(
          condition = "input.geneListChoice == 'deg'",
          selectInput("enrDEGSelect", "Select DE Analysis:", NULL),
          numericInput("enrDEGlog2fc", "Use Log2FC greater than", 0.5),
          numericInput("enrDEGFDR", "Use FDR less than", 0.05),
          checkboxInput("enrDEGUpOnly", "Only use upregulated genes", TRUE),
          numericInput("enrDEGminMean1", "Use mean exp in group1 greater than", 
                       value = 0, min = 0, step = 0.5),
          numericInput("enrDEGmaxMean2", "Use mean exp in group2 less than", 
                       value = 10, min = 0, step = 0.5),
          numericInput("enrDEGminPerc1", "Use exp% in group1 greater than", 
                       value = 0, min = 0, max = 1, step = 0.05),
          numericInput("enrDEGmaxPerc2", "Use exp% in group2 less than", 
                       value = 1, min = 0, max = 1, step = 0.05),
          uiOutput("enrDEGText"),
          verbatimTextOutput(outputId = "enrDEGRes", placeholder = TRUE)
        ),
        conditionalPanel(
          condition = "input.geneListChoice == 'selectGenes'",
          selectizeInput("enrichGenes", label = "Select Gene(s):", NULL, multiple = TRUE)
        ),
        conditionalPanel(
          condition = "input.geneListChoice == 'geneFile'",
          fileInput(
            "enrFile",
            tags$b(tags$i("Please upload a file with only gene names or Entrez Gene Symbols")),
            accept = c("text/csv", "text/comma-separated-values",
                       "text/tab-separated-values", "text/plain",
                       ".csv", ".tsv")
          ),
          h5("Sample files:"),
          tags$a(href = "https://drive.google.com/open?id=1iJZ6H_G2brbeww9B0dA5seMyYUZYyFrU",
                 "Gene Names", target = "_blank"
          ),
          tags$a(href = "https://drive.google.com/open?id=1BLrwW0uMi2pxsX0m1zJrlOTnBLkIiOhk",
                 "Entrez ids", target = "_blank"
          ),
          h5("Options:"),
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
                       inline = TRUE),
          selectInput("enrFileBy", "Match feature type:",
                      choices = c("rownames", featureChoice))
        ),
        hr(),
        selectInput("enrFeatureName", 
                    label = "Select symbol annotation:",
                    choices = c("rownames", featureChoice)),
        selectizeInput("enrichDb", label = "Select DB:", enrichedDB,
                       multiple = TRUE, 
                       options = list(placeholder = "Use all (default)")),
        textInput('enrAnalysisNameSet', "Set analysis name:", 
                  placeholder = "Required"),
        withBusyIndicatorUI(actionButton("enrichRun", "Run")),
      ),
      mainPanel(
        dropdown(
          fluidRow(
            selectInput('enrAnalysisNameSel', label = "Select analysis name:", 
                        choices = NULL),
            selectizeInput('enrDbShow', 
                           label = "Show result of specific database:", 
                           choices = NULL, multiple = TRUE,
                           options = list(placeholder = "Show all (default)"))
          ),
          inputId = "enrDropDown",
          icon = icon("cog"),
          status = "primary",
          circle = FALSE,
          inline = TRUE
        ),
        hr(),
        DT::dataTableOutput("enrDataTable"),
        br(),
        downloadButton("downloadEnrichR", "Download results")
      )
    )
  )
)

