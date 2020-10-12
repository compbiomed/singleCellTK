shinyPanelVis <- fluidPage(
  tags$div(
    class = "container",
    h1("Gene Visualization"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v05-tab03_Dimensionality-Reduction-and-Clustering.html",
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        h3("Options:"),
        selectInput("visAssaySelect", "Select Assay:", currassays),
        selectInput("visPlotMethod", "Visualization Method:", c("boxplot", "scatterplot", "barplot", "heatmap")),
        selectInput("visCondn", "Condition:", c("none", clusterChoice)),
        helpText("To convert the condition to a factor or a numeric value, Go to Data Summary tab -> Annotation data -> Select condition -> select Field type as 'factor' or 'numeric' accordingly"),
        h3("Choose data source:"),
        radioButtons(
          "visGeneList", label = NULL, c("Select Gene(s)" = "selVisRadioGenes",
                                         "Saved top genes" = "visBiomarker")
          ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == 'selVisRadioGenes'", "visGeneList"),
          selectizeInput("selectvisGenes", label = "Select Gene(s):", NULL, multiple = TRUE)
          ),
        shiny::conditionalPanel(
          helpText("To use this, first run Differential expression and save top genes."),
          helpText("Note: currently selects first 'n' genes from the list"),
          condition = sprintf("input['%s'] == 'visBiomarker'", "visGeneList"),
          uiOutput("visBioGenes")
          ),
        uiOutput("visOptions"),
        withBusyIndicatorUI(actionButton("plotvis", "Plot"))
        ),
      mainPanel(
        fluidRow(
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] != 'heatmap'", "visPlotMethod"),
            plotlyOutput("visPlot1", height = "600px")
          ),
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'heatmap'", "visPlotMethod"),
            plotOutput("visPlot2", height = '600px')
          )
          )
        )
      )
    )
  )

