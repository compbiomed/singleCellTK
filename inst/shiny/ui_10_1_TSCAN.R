# User Interface for TSCAN  ---
shinyPanelTSCAN <- fluidPage(
  # Dropdown closing script ####
  tags$script("Shiny.addCustomMessageHandler('close_dropDownTSCAN', function(x){
                  $('html').click();
                });"),
  tags$script("Shiny.addCustomMessageHandler('close_dropDownTscanDE', function(x){
                  $('html').click();
                });"),
  tags$script("Shiny.addCustomMessageHandler('close_dropDownTscanClusterDEG', function(x){
                  $('html').click();
                });"),

  h1("Trajectory Analysis - TSCAN"),
  h5(tags$a(href = paste0(docs.artPath, "trajectoryAnalysis.html"),
            "(help)", target = "_blank")),
  inlineCSS(list(".panel-danger>.panel-heading" = "background-color:#dcdcdc; color:#000000",
                 ".panel-primary>.panel-heading" = "background-color:#f5f5f5; color:#000000; border-color:#dddddd",
                 ".panel-primary" = "border-color:#dddddd;",
                 ".panel-primary>.panel-heading+.panel-collapse>.panel-body" = "border-color:#dddddd;")),
  bsCollapse(
    id = "TSCANUI",
    open = "Calculate Pseudotime Values",
    bsCollapsePanel(
      # Collapse 1, Get MST ####
      "Calculate Pseudotime Values",
      fluidRow(
        column(
          4,
          panel(
            selectInput("TSCANReddim", "Select input dimension reduction:", currreddim),
            selectInput("TSCANclusterName", "Select clustering result: ",
                        "Auto generate clusters", selected = NULL),
            conditionalPanel(
              condition = 'input.TSCANclusterName == "Auto generate clusters"',
              numericInput(inputId = "seed_TSCAN",
                           label = "Seed value for reproducibility of result:",
                           value = 12345,
                           step = 1)
            ),
            actionButton("TSCANRun", "Run")
          )
        ),
        column(
          8,
          panel(
            fluidRow(
              column(
                width = 3,
                dropdown(
                  fluidRow(
                    column(
                      12,
                      fluidRow(actionBttn(inputId = "closeDropDownTSCAN",
                                          label = NULL, style = "simple",
                                          color = "danger", icon = icon("times"),
                                          size = "xs"),
                               align = "right"),
                      selectInput("TSCANVisRedDim", "Select 2D embedding for visualization:", currreddim),
                      actionBttn(
                        inputId = "TSCANPlot",
                        label = "Update",
                        style = "bordered",
                        color = "primary",
                        size = "sm"
                      )
                    )
                  ),
                  inputId = "dropDownTSCAN",
                  icon = icon("cog"),
                  status = "primary",
                  circle = FALSE,
                  inline = TRUE
                ),
              ),
              column(
                width = 9,
                fluidRow(
                  column(
                    width = 12,
                    h6(
                      "A scatter plot of the selected low-dimensionality representation of the dataset will be generated, with the calculated pseudotime colored on each dot (cell). The MST is also projected to the cells."
                    )
                  ),
                  align="center"
                )
              )
            ),
            hr(),
            shinyjqui::jqui_resizable(plotOutput("TSCANPlot"))
          )
        )
      ),
      style = "primary"
    ),

    bsCollapsePanel(
      # Collapse 2, DEG along selected path ####
      "Identify Genes Differentially Expressed For Path",
      fluidRow(
        column(
          4,
          panel(
            selectizeInput(
              inputId = "TSCANassayselect",
              label = "Select input matrix:",
              choices = NULL,
              selected = NULL,
              multiple = FALSE,
              options = NULL),
            pickerInput("pathIndexx", "Select path terminal node:",
                        choices = "", multiple = FALSE),
            pickerInput("discardCluster",
                        "Select cluster(s) to discard (OPTIONAL):",
                        choices = NULL,
                        selected = NULL,
                        multiple = TRUE,
                        options = list(
                          `none-selected-text` = "No cluster discarded"
                        )),
            actionButton("runTSCANDEG", "Run")
          )
        ),
        column(
          8,
          panel(
            fluidRow(
              column(
                width = 3,
                dropdown(
                  fluidRow(
                    column(
                      12,
                      fluidRow(
                        actionBttn(inputId = "closeDropDownTscanDE",
                                   label = NULL, style = "simple",
                                   color = "danger", icon = icon("times"),
                                   size = "xs"),
                        align = "right"
                      ),
                      selectInput("tscanDEexpPathIndex",
                                  "Select path terminal node:",
                                  choices = "", multiple = FALSE),
                      numericInput(inputId = "tscanDEHMTopGenes",
                                   label = "Number of top features for heatmap",
                                   value = 30,
                                   step = 1),
                      numericInput(inputId = "tscanDERegTopGenes",
                                   label = "Number of top features for regulation plots",
                                   value = 10,
                                   step = 1),
                      selectInput("tscanDEFeatureDisplay",
                                  "Display ID Type",
                                  c("Rownames (Default)",
                                    featureChoice)),
                      actionBttn(
                        inputId = "tscanDEPlot",
                        label = "Update",
                        style = "bordered",
                        color = "primary",
                        size = "sm"
                      )
                    )
                  ),
                  inputId = "dropDownTscanDE",
                  icon = icon("cog"),
                  status = "primary",
                  circle = FALSE,
                  inline = TRUE
                )
              ),
              column(
                width = 9,
                fluidRow(
                  column(
                    width = 12,
                    h6(
                      "Visualization on top genes that have significant expression changes along the pseudotime path of insterest."
                    )
                  ),
                  align = "center"
                )
              )
            ),
            hr(),
            tabsetPanel(
              tabPanel(
                # Tab 2.1, DEG Heatmap ####
                "Heatmap",
                panel(
                  fluidRow(
                    column(
                      width = 12,
                      h6(
                        "A heatmap of the expression of the top DE genes along the path in the cells on the path. "
                      )
                    ),

                  ),
                  hr(),
                  shinyjqui::jqui_resizable(
                    plotOutput(outputId = "heatmapPlot")
                  )
                )
              ),
              tabPanel(
                # Tab 2.2, Gene expression increasing along pseudotime ####
                "Up-regulated Genes",
                panel(
                  fluidRow(
                    column(
                      width = 12,
                      h6(
                        "A cell scatter plot showing the expression change along the pseudotime. Genes with top significance in increasing expression along the pseudotime are displayed."
                      )
                    ),
                    align="center"
                  ),
                  hr(),
                  shinyjqui::jqui_resizable(
                    plotOutput(outputId = "UpregGenesPlot")
                  )
                )
              ),
              tabPanel(
                # Tab 2.3, Gene expression decreasing along pseudotime ####
                "Down-regulated Genes",
                panel(
                  fluidRow(
                    column(
                      width = 12,
                      h6(
                        "A cell scatter plot showing the expression change along the pseudotime. Genes with top significance in decreasing expression along the pseudotime are displayed."
                      )
                    ),
                    align="center"
                  ),
                  hr(),
                  shinyjqui::jqui_resizable(
                    plotOutput(outputId = "DownregGenesPlot")
                  )
                )
              )
            )
          )
        )
      ),
      style = "primary"
    ),

    bsCollapsePanel(
      # Collapse 3, DEG between branches of selected cluster ####
      "Identify Genes Differentially Expressed For Branched Cluster",
      fluidRow(
        column(
          4,
          panel(
            selectInput("TSCANUseCluster",
                        "Select branched cluster of interest:",
                        choices = NULL),
            selectizeInput(
              inputId = "TSCANBranchAssaySelect",
              label = "Select input matrix:",
              choices = NULL,
              selected = NULL,
              multiple = FALSE,
              options = NULL),
            numericInput(inputId = "fdrThreshold_TSCAN",
                         label = "FDR less than:",
                         value = 0.05,
                         step = 0.01),
            actionButton("findDEGenes", "Run")
          )
        ),
        column(
          8,
          panel(
            fluidRow(
              column(
                width = 3,
                dropdown(
                  fluidRow(
                    actionBttn(inputId = "closeDropDownTscanClusterDEG",
                               label = NULL, style = "simple",
                               color = "danger", icon = icon("times"),
                               size = "xs"),
                    align = "right"
                  ),
                  selectInput("plotTSCANClusterDEG_useCluster",
                              "Select branched cluster of interest:",
                              choices = "", multiple = FALSE),
                  pickerInput("plotTSCANClusterDEG_pathIndex",
                              "Select Path Index:",
                              choices = NULL,
                              choicesOpt = NULL,
                              selected = NULL),
                  selectInput("plotTSCANClusterDEG_useReducedDim",
                              "Select 2D embedding for visualization:",
                              currreddim),
                  numericInput("plotTSCANClusterDEG_topN",
                               label = "Number of top features to plot:",
                               value = 4, min = 1, step = 1),
                  selectInput("plotTSCANClusterDEG_featureDisplay",
                              "Display ID Type",
                              c("Rownames (Default)",
                                featureChoice)),
                  actionBttn(
                    inputId = "plotTSCANClusterDEG",
                    label = "Update",
                    style = "bordered",
                    color = "primary",
                    size = "sm"
                  ),
                  inputId = "dropDownTscanClusterDEG",
                  icon = icon("cog"),
                  status = "primary",
                  circle = FALSE,
                  inline = TRUE
                )
              ),
              column(
                width = 9,
                fluidRow(
                  column(
                    width = 12,
                    h6(
                      "Visualization and tables of top DE genes on different branch path of the cluster of interest."
                    )
                  ),
                  align="center"
                )
              )
            ),
            hr(),
            tabsetPanel(
              tabPanel(
                # Tab 3.1, feature cluster expression scatter ####
                "Top Feature Plot",
                panel(
                  fluidRow(
                    column(
                      width = 12,
                      h6(
                        "Scatter plots on the selected low-dimension representation of cells in the selected cluster, colored by the expression of top features differentially expressed in the selected branch path. Local MST overlaid."
                      )
                    ),
                    align="center"
                  ),
                  hr(),
                  shinyjqui::jqui_resizable(
                    plotOutput(outputId = "tscanCLusterDEG")
                  )
                )
              ),
              tabPanel(
                # Tab 3.2, Datatable ####
                "Top Feature Table",
                panel(
                  fluidRow(
                    column(
                      width = 12,
                      h6(
                        "A table of top features differentially expressed in the selected branch path of the selected cluster, with statistical metrics displayed."
                      )
                    ),
                    align="center"
                  ),
                  hr(),
                  DT::dataTableOutput("tscanCLusterDEGTable")
                )
              ),
              tabPanel(
                # Tab 3.3, scatter plots of pseudotime ####
                "Pseudotime",
                panel(
                  fluidRow(
                    column(
                      width = 12,
                      h6(
                        "Scatter plots on the selected low-dimension representation of cells in the selected cluster, colored by the recomputed pseudotime value on each of the branching path. Local MST overlaid."
                      )
                    ),
                    align="center"
                  ),
                  hr(),
                  shinyjqui::jqui_resizable(
                    plotOutput(outputId = "tscanCLusterPeudo")
                  )
                )
              )
            )
          )
        )
      ),
      style = "primary"
    ),
    bsCollapsePanel(
      # Collanpse 4, free plot ####
      "Plot feature expression on trajectory",
      fluidRow(
        column(
          4,
          panel(
            selectizeInput(
              inputId = "plotTSCANDimReduceFeatures_useAssay",
              label = "Select expression matrix:",
              choices = NULL,
              selected = NULL,
              multiple = FALSE,
              options = NULL),
            selectizeInput("plotTSCANDimReduceFeatures_features",
                           label = "Select feature(s):", NULL,
                           multiple = TRUE),
            selectInput("plotTSCANDimReduceFeatures_useReducedDim",
                        "Select 2D embedding for visualization:", currreddim),
            pickerInput("plotTSCANDimReduceFeatures_useCluster",
                        "Show cluster(s) of interest:",
                        choices = NULL,
                        selected = NULL,
                        multiple = TRUE,
                        options = list(
                          `none-selected-text` = "Show all"
                        )),
            selectInput("plotTSCANDimReduceFeatures_featureDisplay",
                        "Display ID Type",
                        c("Rownames (Default)",
                          featureChoice)),
            actionButton("plotTSCANDimReduceFeatures", "Plot")
          )
        ),
        column(
          8,
          panel(
            fluidRow(
              column(
                width = 12,
                h6(
                  "Scatter plots on the selected low-dimension representation of cells in the selected cluster, colored by the expression of selected features, with the MST overlaid."
                )
              ),
              align="center"
            ),
            hr(),
            shinyjqui::jqui_resizable(
              plotOutput("TscanDimReduceFeatures")
            )
          )
        )
      ),
      style = "primary"
    )
  ),
  nonLinearWorkflowUI(id = "nlw-Traj")
)
