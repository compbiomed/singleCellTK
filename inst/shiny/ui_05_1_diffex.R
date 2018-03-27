shiny_panel_diffex <- fluidPage(
  tags$div(
    class = "container",
    h1("Differential Expression"),
    h5(tags$a(href = "https://compbiomed.github.io/singleCellTK/articles/v07-tab05_Differential-Expression.html",
              "(help)", target="_blank")),
    sidebarLayout(
      sidebarPanel(
        #TODO: Remove DESeq, add edgeR, add more custom options?
        selectInput("diffexAssay", "Select Assay:", currassays),
        selectInput("selectDiffex", "Select Method:", c("limma", "DESeq2",
                                                        "ANOVA")),
        uiOutput("selectDiffex_conditionUI"),
        uiOutput("selectDiffex_conditionlevelUI"),
        sliderInput("selectNGenes", "Display Top N Genes:", 5, 500, 500, 5),
        checkboxInput("applyCutoff", "Apply p-value Cutoff"),
        conditionalPanel(
          condition = "input.applyCutoff == true",
          sliderInput("selectPval", "p-value cutoff:", 0.01, 0.2, 0.05),
          selectInput("selectCorrection", "Correction Method:", c("FDR"))
        ),
        withBusyIndicatorUI(actionButton("runDiffex",
                                         "Run Differential Expression")),
        downloadButton("downloadGeneList", "Download Results"),
        h3("Save gene list as biomarker:"),
        textInput("biomarkerName", "Biomarker Name: ", value = ""),
        withBusyIndicatorUI(actionButton("saveBiomarker", "Save Biomarker"))
      ),
      mainPanel(
        tabsetPanel(
          id = "dataset",
          tabPanel(
            "Heatmap",
            br(),
            tabsetPanel(
              tabPanel("Heatmap", plotOutput("diffPlot")),
              tabPanel(
                "Options",
                wellPanel(
                  h3("General Options"),
                  checkboxInput("displayHeatmapRowLabels", "Display Row Labels",
                                value = TRUE),
                  checkboxInput("displayHeatmapColumnLabels",
                                "Display Column Labels", value = TRUE),
                  checkboxInput("displayHeatmapColumnDendrograms",
                                "Display Column Dendrograms", value = TRUE),
                  checkboxInput("displayHeatmapRowDendrograms",
                                "Display Row Dendrograms", value = TRUE),
                  checkboxInput("clusterRows", "Cluster Heatmap Rows", value = TRUE),
                  checkboxInput("clusterColumns", "Cluster Heatmap Columns",
                                value = TRUE),
                  textInput("heatmapColumnsTitle", "Columns Title",
                            value = "Differential Expression"),
                  tags$hr(),
                  h3("Colorbar Options"),
                  checkboxInput("displayHeatmapColorBar", "Color Bar",
                                value = TRUE),
                  uiOutput("colorBarCondition"),
                  uiOutput("HeatmapSampleAnnotations")
                )
              )
            )
          ),
          tabPanel(
            "Results Table",
            dataTableOutput("diffextable")
          )
        )
      )
    )
  )
)
