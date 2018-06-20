shinyPanelDiffex <- fluidPage(
  tags$div(
    class = "container",
    h1("Differential Expression"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v07-tab05_Differential-Expression.html",
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        #TODO: Remove DESeq, add edgeR, add more custom options?
        selectInput("diffexAssay", "Select Assay:", currassays),
        selectInput("selectDiffex", "Select Method:", c("limma (use log values)" = "limma",
                                                        "DESeq2 (use counts)" = "DESeq2",
                                                        "ANOVA (use log values)" = "ANOVA")),
        uiOutput("selectDiffexConditionUI"),
        uiOutput("selectDiffexConditionLevelUI"),
        withBusyIndicatorUI(actionButton("runDiffex",
                                         "Run Differential Expression")),
        tags$hr(),
        h3("Plot Options:"),
        sliderInput("selectNGenes", "Display Top N Genes:", 5, 500, 500, 5),
        checkboxInput("applyCutoff", "Apply p-value Cutoff"),
        conditionalPanel(
          condition = "input.applyCutoff == true",
          sliderInput("selectPval", "p-value (adjusted) cutoff:", 0.01, 0.2, 0.05),
          selectInput("selectCorrection", "Correction Method:", c("fdr", "holm",
                                                                  "hochberg",
                                                                  "hommel",
                                                                  "bonferroni",
                                                                  "BH", "BY",
                                                                  "none"))
        ),
        tags$hr(),
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
                  uiOutput("colorBarConditionUI"),
                  uiOutput("HeatmapSampleAnnotations")
                )
              )
            )
          ),
          tabPanel(
            "Results Table",
            DT::dataTableOutput("diffextable")
          )
        )
      )
    )
  )
)
