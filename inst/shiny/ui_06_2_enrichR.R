shinyPanelEnrichR <- fluidPage(
  tags$div(
    class = "container",
    h1("Gene Enrichment Analysis using enrichR"),
    h5(tags$a(href = "", "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        selectInput("enrichAssay", "Select Assay", currassays),
        selectizeInput("enrichGenes", label = "Select Gene(s):", NULL, multiple = TRUE),
        selectizeInput("enrichDb", label = "Select DB:", c("ALL", enrichedDB), multiple = TRUE),
        withBusyIndicatorUI(actionButton("enrichRun", "Run")),
        br(),
        downloadButton("downloadEnrichR","Download results")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Results", DT::dataTableOutput("enrichTable"))
        )
      )
    )
  )
)
