shiny_panel_diffex <- fluidPage(
  tags$div(
    class="container",
    h1("Differential Expression"),
    fluidPage(
      fluidRow(
        column(
          4,
          wellPanel(
            selectInput("selectDiffex","Differential Expression",c("limma", "DESeq", "DESeq2")),
            selectInput("selectDiffex_condition","Select Condition",clusterChoice),
            uiOutput("selectDiffex_conditionofinterestUI"),
            sliderInput("selectNGenes", "Display Top N Genes:", 5, 500, 500, 5),
            checkboxInput("applyCutoff", "Apply p-value Cutoff"),
            checkboxInput("clusterRows", "Cluster Heatmap Rows", value=TRUE),
            checkboxInput("clusterColumns", "Cluster Heatmap Columns", value=TRUE),
            sliderInput("selectPval", "p-value cutoff:", 0.01, 0.2, 0.05),
            selectInput("selectCorrection","Correction Type",c("FDR")),
            withBusyIndicatorUI(actionButton("runDiffex", "Run Differential Expression")),
            downloadButton("downloadGeneList","Download Results"),
            h3("Save gene list as biomarker:"),
            textInput("biomarkerName", "Biomarker Name: ", value = ""),
            actionButton("saveBiomarker","Save Biomarker")
          )
        ),
        column(
          8,
          tabsetPanel(
            id = 'dataset',
            tabPanel('Heatmap', 
                     fluidPage(
                       fluidRow(
                         column(
                           4,
                           wellPanel(
                             "General Options",
                             checkboxInput("displayHeatmapRowLabels", "Display Row Labels", value=TRUE),
                             checkboxInput("displayHeatmapColumnLabels", "Display Column Labels", value=TRUE),
                             checkboxInput("displayHeatmapColumnDendrograms", "Display Column Dendrograms", value=TRUE),
                             checkboxInput("displayHeatmapRowDendrograms", "Display Row Dendrograms", value=TRUE),
                             textInput("heatmapColumnsTitle", "Columns Title", value = "Differential Expression")
                           ),
                           wellPanel(
                             "Colorbar Options",
                             checkboxInput("displayHeatmapColorBar", "Color Bar", value=TRUE),
                             uiOutput("colorBarCondition"),
                             uiOutput("HeatmapSampleAnnotations")
                           )
                         ),
                         column(
                           8,
                           plotOutput("diffPlot")
                         )
                       )
                     )),
            tabPanel('Results Table', dataTableOutput('diffextable')),
            tabPanel('Interactive Heatmap', d3heatmapOutput("interactivediffPlot"))
          )
        )
      )
    )
  ),
  includeHTML('www/footer.html')
)
