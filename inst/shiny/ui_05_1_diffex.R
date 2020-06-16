shinyPanelDiffex <- fluidPage(
  tags$div(
    class = "container",
    h1("Differential Expression"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v07-tab05_Differential-Expression.html",
              "(help)", target = "_blank")),
    sidebarLayout(
      sidebarPanel(
        fluidRow(
          column(5, h3("Settings:")),
          column(3, br(), actionButton("Diffex_hideAllSections", "Hide All")),
          column(3, br(), actionButton("Diffex_showAllSections", "Show All"))
        ),
        br(),
        #TODO: Remove DESeq, add edgeR, add more custom options?
        actionButton("diffex1", "Options"),
        tags$div(
          id = "de1",
          wellPanel(
            actionButton("diffex6", "Saved Results"),
            shinyjs::hidden(
              tags$div(
                id = "de6",
                wellPanel(
                  uiOutput("savedRes"),
                  withBusyIndicatorUI(
                    actionButton("loadResults", "Load Saved Results")
                  )
                )
              )
            ),
            selectInput("diffexAssay", "Select Assay:", currassays),
            selectInput("selectDiffex", "Select Method:",
                        c("limma (use log values)" = "limma",
                          "DESeq2 (use counts)" = "DESeq2",
                          "ANOVA (use log values)" = "ANOVA")),
            uiOutput("selectDiffexConditionUI"),
            uiOutput("selectDiffexConditionLevelUI"),
            selectInput("selectCorrection", "Correction Method:",
                        c("fdr", "holm", "hochberg", "hommel", "bonferroni",
                          "BH", "BY", "none")),
            withBusyIndicatorUI(actionButton("runDiffex",
                                             "Run Differential Expression"))
          )
        ),
        actionButton("diffex2", "Heatmap"),
        shinyjs::hidden(
          tags$div(
            id = "de2",
            wellPanel(
              numericInput("selectNGenes", "Enter Top 'N' Genes value:", value = 500),
              uiOutput("diffexNgenes"),
              checkboxInput("applyCutoff", "Apply p-value Cutoff"),
              conditionalPanel(
                condition = "input.applyCutoff == true",
                sliderInput("selectPval", "p-value (adjusted) cutoff:", 0.01, 0.2, 0.05)
              ),
              checkboxInput("applylogFCCutoff", "Apply logFC Cutoff"),
              conditionalPanel(
                condition = "input.applylogFCCutoff == true",
                numericInput("selectlogFCDiffex", "Select logFC cutoff", value = 2, step = 0.5),
                uiOutput("logFCDiffexRange"),
                checkboxInput("applyAbslogFCDiffex", "absolute logFC value")
              ),
              checkboxInput("applyScaleDiffex", "Scale Expression values?", value = TRUE),
              br(),
              actionButton("diffex3", "Options"),
              tags$div(
                id = "de3",
                wellPanel(
                  h3("General Options"),
                  checkboxInput("displayHeatmapRowLabels",
                                "Row Labels", value = TRUE),
                  checkboxInput("displayHeatmapColumnLabels",
                                "Column Labels", value = TRUE),
                  checkboxInput("displayHeatmapColumnDendrograms",
                                "Column Dendrograms", value = TRUE),
                  checkboxInput("displayHeatmapRowDendrograms",
                                "Row Dendrograms", value = TRUE),
                  checkboxInput("clusterRows", "Cluster Rows",
                                value = TRUE),
                  checkboxInput("clusterColumns", "Cluster Columns",
                                value = TRUE),
                  textInput("heatmapColumnsTitle", "Columns Title",
                            value = "Differential Expression"),
                  tags$hr(),
                  h3("Colorbar Options"),
                  checkboxInput("displayHeatmapColorBar", "Display Color Bar",
                                value = TRUE),
                  uiOutput("colorBarConditionUI"),
                  uiOutput("heatmapSampleAnnotations")
                )
              ),
              br(),
              withBusyIndicatorUI(actionButton("runPlotDiffex", "Plot heatmap"))
            )
          )
        ),
        actionButton("diffex5", "Save Results"),
        shinyjs::hidden(
          tags$div(
            id = "de5",
            wellPanel(
              helpText("Save results in rowData() of a sctk object to use it in a current session or download the sctk object for future use"),
              textInput("ResultsName", "Name : ", value = ""),
              withBusyIndicatorUI(actionButton("saveResults", "Save Results")),
              uiOutput("saveDiffResultsNote")
            )
          )
        ),
        #this section was previously called 'biomarker'
        actionButton("diffex7", "Save top 'N' significant genes"),
        shinyjs::hidden(
          tags$div(
            id = "de7",
            wellPanel(
             helpText("Save these to be used in enrichment analysis, visualise and filtering gene annotation"),
              br(),
              numericInput("selectBioNGenes", "Enter Top 'N' Genes value: ", value = 100),
              uiOutput("BioNgenes"),
              checkboxInput("applyBioCutoff1", "Apply p-value Cutoff"),
              conditionalPanel(
                condition = "input.applyBioCutoff1 == true",
                sliderInput("selectAdjPVal", "Select p-value(adjusted) cutoff", 0.01, 0.2, 0.05)
              ),
              checkboxInput("applyBioCutoff2", "Apply logFC Cutoff"),
              conditionalPanel(
                condition = "input.applyBioCutoff2 == true",
                numericInput("selectlogFC", "Select logFC cutoff", value = 2, step = 0.5),
                uiOutput("logFCBioRange"),
                checkboxInput("applyAbslogFC", "absolute logFC value")
              ),
              textInput("biomarkerName", "Name : ", value = ""),
              withBusyIndicatorUI(actionButton("saveBiomarker", "Save Biomarker")),
              uiOutput("bioMarkerNote")
            )
          )
        ),
        actionButton("diffex4", "Download Results table"),
        shinyjs::hidden(
          tags$div(
            id = "de4",
            wellPanel(downloadButton("downloadGeneList", "Download(.csv)"))
          )
        )
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Results Table",
            DT::dataTableOutput("diffextable")
          ),
          tabPanel(
            "Heatmap",
            plotOutput("diffPlot")
          )
        )
      )
    )
  )
)

