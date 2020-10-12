ShinyPanelHypeR <- fluidPage(
  tags$div(
    class = "container",
    h1("Enrichment Analysis using hypeR"),
    h5(tags$a(href = "https://github.com/montilab/hypeR", "(help)", target = "_blank")),
    sidebarLayout(
    sidebarPanel(
      # Put your geneset selector module anywhere
      hypeR::genesets_UI("genesets"),
      tags$br(),
      HTML(paste(tags$h6(tags$i("Note: Please wait until a checkmark appears after clicking.")))),
      selectizeInput("hypgenes", label = "Signature", NULL, multiple=TRUE),
      actionButton("enrichment", "Enrichment")
    ),
    mainPanel(
      wellPanel(
        tags$h4(tags$b("Enrichment Results"), align="center"),
        reactable::reactableOutput("hyptab", height = "400px"),
        shiny::plotOutput("plot", height = "400px")
      )
    )
  )
)
)

