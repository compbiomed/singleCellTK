orgpkgs <- c("Anopheles" = "org.Ag.eg.db", "Arabidopsis" = "org.At.tair.db",
             "Bovine" = "org.Bt.eg.db", "Worm" = "org.Ce.eg.db",
             "Canine" = "org.Cf.eg.db", "Fly" = "org.Dm.eg.db",
             "Zebrafish" = "org.Dr.eg.db",
             "E coli strain K12" = "org.EcK12.eg.db",
             "E coli strain Sakai" = "org.EcSakai.eg.db",
             "Chicken" = "org.Gg.eg.db", "Human" = "org.Hs.eg.db",
             "Mouse" = "org.Mm.eg.db", "Rhesus" = "org.Mmu.eg.db",
             "Malaria" = "org.Pf.plasmo.db", "Chimp" = "org.Pt.eg.db",
             "Rat" = "org.Rn.eg.db", "Yeast" = "org.Sc.sgd.db",
             "Pig" = "org.Ss.eg.db", "Xenopus" = "org.Xl.eg.db")

shinyPanelFilter <- fluidPage(
  useShinyalert(),
  tags$div(
    class = "container",
    h4("QC & Filtering"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v04-tab02_Data-Summary-and-Filtering.html",
              "(help)", target = "_blank")),
    tabsetPanel(
      tabPanel(
        "Data Summary",
        wellPanel(
          sidebarLayout(
            sidebarPanel(
              fluidRow(
                column(5,
                       h3("Settings:")
                ),
                column(3,
                       br(),
                       actionButton("f_hideAllSections", "Hide All")
                ),
                column(3,
                       br(),
                       actionButton("f_showAllSections", "Show All")
                )
              ),
              br(),
              # SHINYJS COLLAPSE --------------------------
              # Section 1 - Assay Settings
              actionButton("f_button1", "Assay Settings"),
              # open by default
              tags$div(id = "f_collapse1",
                       br(),
                       wellPanel(
                         selectInput("filterAssaySelect", "Select Assay:", currassays)
                       )
              ),
              # Section 2 - Delete Outliers
              actionButton("f_button2", "Delete Outliers"),
              shinyjs::hidden(
                tags$div(id = "f_collapse2",
                         wellPanel(
                           checkboxInput("removeNoexpress", "Remove genes with 0 expression across all Cells (Recommended)", value = TRUE),
                           numericInput("minDetectGene", label = "Minimum Detected Genes per Cells", value = 1700, min = 1, max = 100000),
                           numericInput("LowExpression", "% Low Gene Expression to Filter", value = 40, min = 0, max = 100),
                           selectInput("deletesamplelist", "Select Cells:", sampleChoice, multiple = TRUE),
                           fluidRow(
                             column(6, withBusyIndicatorUI(actionButton("filterData", "Filter Data"))),
                             column(6, actionButton("resetData", "Reset All"))
                           )
                         )
                )
              ),
              # Section 3 - Filter Samples by Annotation
              actionButton("f_button3", "Filter Cells by Annotation"),
              shinyjs::hidden(
                tags$div(id = "f_collapse3",
                         wellPanel(
                           selectInput("filteredSample", "Select Annotation:", c("none", clusterChoice)),
                           uiOutput("filterSampleOptions")
                         )
                )
              ),
              # Section 4 - Filter Genes by Feature Annotation
              actionButton("f_button4", "Filter Genes by Feature Annotation"),
              shinyjs::hidden(
                tags$div(id = "f_collapse4",
                         wellPanel(
                           selectInput("filteredFeature", "Select Feature:", c("none", featureChoice)),
                           uiOutput("filterFeatureOptions")
                         )
                )
              ),
              # Section 5 - Convert Gene Annotations
              actionButton("f_button5", "Convert Gene Annotations"),
              shinyjs::hidden(
                tags$div(id = "f_collapse5",
                         wellPanel(
                           selectInput("orgOrganism", "Select Organism:",
                                       orgpkgs[orgpkgs %in% as.character(
                                         grep("^org\\.", installed.packages()[, "Package"],
                                              value = TRUE))]),
                           uiOutput("orgConvertColumns"),
                           withBusyIndicatorUI(actionButton("convertGenes", "Convert"))
                         )
                )
              ),
              # Section 6 - Delete an Annotation Column
              actionButton("f_button6", "Delete an Annotation Column"),
              shinyjs::hidden(
                tags$div(id = "f_collapse6",
                         wellPanel(
                           selectInput("deleterowdatacolumn", "Annotation Column:", clusterChoice),
                           actionButton("deleterowDatabutton", "Delete Column")
                         )
                )
              ),
              # Section 7 - Randomly Subset
              actionButton("f_button7", "Randomly Subset"),
              shinyjs::hidden(
                tags$div(id = "f_collapse7",
                         wellPanel(
                           numericInput("downsampleNum", "Number of cells to keep:", min = 2,
                                        max = numSamples, value = numSamples, step = 1),
                           withBusyIndicatorUI(actionButton("downsampleGo", "Subset Data"))
                         )
                )
              ),
              tags$hr(),
              downloadButton("downloadSCE", "Download SCtkExperiment(.rds)")
            ),
            mainPanel(
              wellPanel(
                style = "background-color:transparent",
                h4("Summary Contents:"),
                tableOutput("summarycontents"),
                tags$hr(),
                h4("Counts Histogram:"),
                plotlyOutput("countshist"),
                tags$hr(),
                h4("Genes Histogram:"),
                plotlyOutput("geneshist"),
                tags$hr(),
                h4("Data Table:"),
                helpText("Note: Shows table only if cells < 50"),
                DT::dataTableOutput("contents")
              )
            )
          )
        )
      ),
      tabPanel(
        "Assay Normalization",
        wellPanel(
          br(),
          #fluidRow(
            sidebarLayout(
              sidebarPanel(
                fluidRow(
                  column(5, h3("Settings:")),
                  column(3, br(), actionButton("Norm_hideAllSections", "Hide All")),
                  column(3, br(), actionButton("Norm_showAllSections", "Show All"))
                ),
                br(),
                actionButton("norm1", "Basic Options"),
                shinyjs::hidden(
                tags$div(
                  id = "nm1",
                  wellPanel(
                    selectInput("assayModifyAction", "Assay Actions:",
                                c("Log Transform" = "log", "Create CPM" = "cpm",
                                  "Rename" = "rename", "Delete" = "delete")),
                    conditionalPanel(
                      condition = "input.assayModifyAction == 'cpm'",
                      h5("Select a count assay to use for CPM calculation:")
                    ),
                    selectInput("modifyAssaySelect", "Select Assay:", currassays),
                    conditionalPanel(
                      condition = "input.assayModifyAction != 'delete'",
                      textInput("modifyAssayOutname", "Assay Name", "",
                                placeholder = "What should the assay be called?")
                    ),
                    withBusyIndicatorUI(actionButton("modifyAssay", "Run"))
                    )
                  )
                ),
                  actionButton("norm2", "Advanced Options"),
                  shinyjs::hidden(
                  tags$div(
                    id = "nm2",
                    wellPanel(
                      selectInput("normalizeAssay", "Normalization:",
                                  c("Library Size normalization (by scater)" = "libSizeNorm",
                                    "Normalization by deconvolution (by scran)" = "deconvoNorm")),
                      conditionalPanel(
                        condition = "input.normalizeAssay == 'libSizeNorm'",
                        selectInput("libSizeNormAssaySelect", "Assay:", currassays),
                        helpText("Note: select 'counts' assay"),
                        textInput("libSizeAssayName", "Name:")
                      ),
                      conditionalPanel(
                        condition = "input.normalizeAssay == 'deconvoNorm'",
                        selectInput("deconvoAssaySelect", "Assay:", currassays),
                        textInput("clusterMin", "Min no. of cells per cluster", "10"),
                        textInput("deconvoAssayName", "Name:")
                      ),
                      withBusyIndicatorUI(actionButton("modifyNorm", "Run"))
                    )
                  )
                )
              ),
              mainPanel(
                fluidRow(
                  column(6,
                         h4("Available Assays:"),
                         tableOutput("assayList")
                         )
                )
              )
            )
        )
      ),
      tabPanel(
        "Annotation Data",
        wellPanel(
          sidebarLayout(
            sidebarPanel(
          h4("Modify Annotation Data:"),
          selectInput("annotModifyChoice", "Select Annotation:", c("none", clusterChoice)),
          uiOutput("annotModifyUI"),
          uiOutput("annotModifyUIHelpText"),
          tags$hr(),
          downloadButton("downloadcolData", "Download Annotation Data"),
          tags$hr(),
          fileInput(
            "newAnnotFile", "Upload and replace annotation data:",
            accept = c(
              "text/csv",
              "text/comma-separated-values",
              ".csv"
              )
            )
          ),
          mainPanel(
            tags$h4("Column Annotation data(colData):"),
            tags$br(),
            DT::dataTableOutput("colDataDataFrame")
            )
          )
        )
      ),
      tabPanel(
        "Batch Correction",
        wellPanel(
          sidebarLayout(
            sidebarPanel(
              selectInput("combatAssay", "Select Assay:", currassays),
              tags$hr(),
              h4("Plot Batch Effect:"),
              selectInput("batchVarPlot", "Select Batch Annotation:", c("none", clusterChoice)),
              selectInput("conditionVarPlot", "Select Condition Annotation:", c("none", clusterChoice)),
              tags$hr(),
              h4("Run Batch Correction:"),
              selectInput("batchMethod", "Select Method:", "ComBat"),
              selectInput("combatBatchVar", "Select Batch Condition:", clusterChoice),
              selectInput("combatConditionVar", "Select Additional Covariates:",
                          clusterChoice, multiple = TRUE),
              radioButtons("combatParametric", "Adjustments:", c("Parametric",
                                                                 "Non-parametric"),
                           selected = "Parametric"),
              checkboxInput("combatMeanOnly", "Correct mean of the batch effect only",
                            value = FALSE),
              checkboxInput("combatRef", "Run reference batch combat:",
                            value = FALSE),
              uiOutput("selectCombatRefBatchUI"),
              textInput("combatSaveAssay", "Assay Name to Use:", value = "combat"),
              withBusyIndicatorUI(actionButton("combatRun", "Run"))
            ),
            mainPanel(
              uiOutput("combatStatus"),
              plotOutput("combatBoxplot", height = "600px")
            )
          )
        )
      )
    )
  )
)
