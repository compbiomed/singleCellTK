# User Interface for Scanpy Workflow ---
shinyPanelScanpy <- fluidPage(
    tags$script("Shiny.addCustomMessageHandler('close_dropDownScanpy', function(x){
                  $('html').click();
                });"), # close_dropDownSeuratHM
    h1("Scanpy"),
    h5(tags$a(href = paste0(docs.artPath, "cnsl_scanpy_curated_workflow.html"), #cnsl_seurat_curated_workflow.html
              "(help)", target = "_blank")),
    inlineCSS(list(".panel-danger>.panel-heading" = "background-color:#dcdcdc; color:#000000", ".panel-primary>.panel-heading" = "background-color:#f5f5f5; color:#000000; border-color:#dddddd", ".panel-primary" = "border-color:#dddddd;", ".panel-primary>.panel-heading+.panel-collapse>.panel-body" = "border-color:#dddddd;")),
    conditionalPanel(
        condition = "false",
        selectInput(
            "activePanelSelectScanpy", #activePanelSelectSeurat
            label = "Active Panel:",
            choices = c("",
                        "Normalize Data",
                        "Highly Variable Genes",
                        "Dimensionality Reduction",
                        "tSNE/UMAP",
                        "Clustering",
                        "Find Markers",
                        "Heatmap Plot"),
            selected = ""
        )
    ),
    bsCollapse(id = "ScanpyUI", open = "scanpy QC & Filtering", # ScanPy
            bsCollapsePanel("Normalize Data",
                                                                   fluidRow(
                                                                     column(4,
                                                                            panel(heading = "Options",
                                                                                  selectizeInput(
                                                                                    inputId = "scanpySelectNormalizationAssay", # seuratSelectNormalizationAssay
                                                                                    label = "Select input matrix:",
                                                                                    choices = NULL,
                                                                                    selected = NULL,
                                                                                    multiple = FALSE,
                                                                                    options = NULL),
                                                                                  numericInput(inputId = "scanpy_targetSum", label = "Specify targetSum:", value = 1),
                                                                                  numericInput(inputId = "scanpy_maxFraction", label = "Specify maxFraction:", value = 0.05),
                                                                                  actionButton(inputId = "scanpy_normalize_button", "Normalize") # normalize_button
                                                                            )
                                                                     )
                                                                   ),
                                                                   style = "primary"
            ),
            bsCollapsePanel("Highly Variable Genes",
                            fluidRow(
                              column(4,
                                     fluidRow(
                                       column(12,
                                              panel(heading = "Compute HVG",
                                                    selectInput(inputId = "scanpy_hvg_method", label = "Select HVG method: ", choices = c("seurat", "cell_ranger", "seurat_v3")), # hvg_method
                                                    numericInput(inputId = "scanpy_hvg_no_features", label = "Select number of features to find: ", value = "2000"), # hvg_no_features
                                                    numericInput(inputId = "scanpy_minMean", label = "Specify minMean: ", value = "0.0125"),
                                                    numericInput(inputId = "scanpy_maxMean", label = "Specify maxMean: ", value = "3"),
                                                    numericInput(inputId = "scanpy_minDisp", label = "Specify minDisp: ", value = "0.5"),
                                                    checkboxInput(inputId = "scanpy_maxDisp_Inf", label = "Set maxDisp to Infinite:", value = TRUE),
                                                    conditionalPanel(
                                                      condition = 'input.scanpy_maxDisp_Inf == false',
                                                    numericInput(inputId = "scanpy_maxDisp", label = "Specify maxDisp: ", value = 3)
                                                    ),
                                                    actionButton(inputId = "scanpy_find_hvg_button", "Find HVG") # find_hvg_button
                                              )
                                       )
                                     ),
                                     br(),
                                     fluidRow(
                                       column(12,
                                              panel(heading = "Display HVG",
                                                    numericInput(inputId = "scanpy_hvg_no_features_view", label = "Select number of features to display: ", value = 10, step = 1), # hvg_no_features_view
                                                    verbatimTextOutput(outputId = "scanpy_hvg_output", placeholder = TRUE) # hvg_output
                                              )
                                       )
                                     )
                              ),
                              column(8,
                                     fluidRow(
                                       column(12,
                                              panel(heading = "Plot",
                                                    plotOutput(outputId = "scanpy_plot_hvg") # plot_hvg
                                              )
                                       )
                                     )
                              )
                            ),
                            style = "primary"
            ),
            bsCollapsePanel("Dimensionality Reduction",
                            fluidRow(
                              column(4,
                                     fluidRow(
                                       column(12,
                                              panel(heading = "scanpy PCA",
                                                    selectInput(inputId = "scanpy_pca_method", label = "Select PCA method: ", choices = c("arpack", "randomized", "auto", "lobpcg")),
                                                    numericInput(inputId = "scanpy_pca_no_components", label = "Select number of components to compute: ", value = 50), # pca_no_components
                                                    # materialSwitch(inputId = "scanpy_pca_compute_elbow", label = "Compute ElbowPlot?", value = TRUE),
                                                    # materialSwitch(inputId = "scanpy_pca_compute_jackstraw", label = "Compute JackStrawPlot?", value = FALSE),
                                                    # materialSwitch(inputId = "scanpy_pca_compute_heatmap", label = "Compute Heatmap?", value = TRUE),
                                                    # conditionalPanel(
                                                    #   condition = 'input.pca_compute_heatmap == true',
                                                    #   numericInput(inputId = "scanpy_pca_compute_heatmap_nfeatures",
                                                    #                label = "Set number of features for heatmap:", value = 30, step = 1),
                                                    # ),
                                                    # numericInput(inputId = "scanpy_seed_PCA",
                                                    #              label = "Seed value for reproducibility of result:",
                                                    #              value = 42,
                                                    #              step = 1),
                                                    actionButton(inputId = "scanpy_run_pca_button", "Run PCA")
                                              ),
                                              panel(heading = "Select No. of Components",
                                                    htmlOutput(outputId = "scanpy_pca_significant_pc_output", inline = FALSE),
                                                    numericInput(inputId = "scanpy_pca_significant_pc_counter", label = "Select number of components for downstream analysis: ", min = 1, max = 20, value = 10)
                                              )
                                       )
                                     )
                              ),
                              column(8,
                                     fluidRow(
                                       column(12,
                                              hidden(
                                                tags$div(class = "scanpy_pca_plots", tabsetPanel(id = "scanpyPCAPlotTabset", type = "tabs"
                                                )
                                                ))
                                       )
                                       
                                     )
                              )
                            ),
                            style = "primary"
            ),
            bsCollapsePanel("tSNE/UMAP",
                            tabsetPanel(id = "tsneUmapTabsetScanpy", type = "tabs",
                                        tabPanel("scanpy tSNE",
                                                 br(),
                                                 fluidRow(
                                                   column(4,
                                                          fluidRow(
                                                            column(12,
                                                                   panel(heading = "scanpy tSNE",
                                                                         selectInput(inputId = "scanpy_reduction_tsne_method", label = "Select reduction method: ", choices = c("pca", "ica")),
                                                                         #textInput(inputId = "reduction_tsne_count", label = "Select number of reduction components: ", value = "20"),
                                                                         numericInput(inputId = "scanpy_perplexity_tsne", label = "Set perplexity:", value = 15),
                                                                         # numericInput(inputId = "scanpy_seed_TSNE",
                                                                         #              label = "Seed value for reproducibility of result:",
                                                                         #              value = 1,
                                                                         #              step = 1),
                                                                         # htmlOutput(outputId = "scanpy_display_message_tsne", inline = FALSE),
                                                                         actionButton(inputId = "scanpy_run_tsne_button", "Run tSNE")
                                                                   )
                                                            )
                                                          )
                                                   ),
                                                   column(8,
                                                          fluidRow(
                                                            panel(heading = "scanpy Plot",
                                                                  column(12,
                                                                         plotOutput(outputId = "scanpy_plot_tsne")
                                                                  )
                                                            )
                                                          )
                                                   )
                                                 )
                                        ),
                                        tabPanel("scanpy UMAP",
                                                 br(),
                                                 fluidRow(
                                                   column(4,
                                                          fluidRow(
                                                            column(12,
                                                                   panel(heading = "scanpy UMAP",
                                                                         #selectInput(inputId = "scanpy_reduction_umap_method", label = "Select reduction method: ", choices = c("pca", "ica")),
                                                                         numericInput(inputId = "scanpy_reduction_umap_count", label = "Select number of reduction components: ", value = 10),
                                                                         numericInput(inputId = "scanpy_min_dist_umap", label = "Set min.dist:", value = 0.5),
                                                                         numericInput(inputId = "scanpy_n_neighbors_umap", label = "Set n.neighbors:", value = 15, step = 1),
                                                                         numericInput(inputId = "scanpy_spread_umap", label = "Set spread:", value = 1),
                                                                         numericInput(inputId = "scanpy_spread_alpha", label = "Set alpha:", value = 1),
                                                                         numericInput(inputId = "scanpy_spread_gamma", label = "Set gamma:", value = 1),
                                                                         # numericInput(inputId = "scanpy_seed_UMAP",
                                                                         #              label = "Seed value for reproducibility of result:",
                                                                         #              value = 42,
                                                                         #              step = 1),
                                                                         # htmlOutput(outputId = "scanpy_display_message_umap", inline = FALSE),
                                                                         actionButton(inputId = "scanpy_run_umap_button", "Run UMAP")
                                                                   )
                                                            )
                                                          )
                                                   ),
                                                   column(8,
                                                          fluidRow(
                                                            panel(heading = "scanpy Plot",
                                                                  column(12,
                                                                         plotOutput(outputId = "scanpy_plot_umap")
                                                                  )
                                                            )
                                                          )
                                                   )
                                                 )
                                        )
                            ),
                            style = "primary"
            ),
            bsCollapsePanel("Clustering",
                            fluidRow(
                              column(4,
                                     fluidRow(
                                       column(12,
                                              panel(heading = "Options",
                                                    #selectInput(inputId = "scanpy_reduction_clustering_method", label = "Select reduction method: ", choices = c("pca", "ica")),
                                                    numericInput(inputId = "scanpy_reduction_clustering_count", label = "Select number of reduction components: ", value = 2),
                                                    selectInput(inputId = "scanpy_algorithm.use", label = "Select clustering algorithm: ", choices = list("louvain algorithm" = "louvain",
                                                                                                                                                   "leiden algorithm" = "leiden")),
                                                    numericInput(inputId = "scanpy_resolution_clustering", label = "Set resolution:", value = 0.8),
                                                    numericInput(inputId = "scanpy_nNeighbors", label = "Set no. of neigbours:", value = 15),
                                                    # materialSwitch(inputId = "scanpy_group.singletons", label = "Group singletons?", value = TRUE),
                                                    # htmlOutput(outputId = "scanpy_display_message_clustering", inline = FALSE),
                                                    selectInput(inputId = "scanpy_corr_method", label = "Select clustering algorithm: ", choices = list("pearson" = "pearson",
                                                                                                                                                          "kendall" = "kendall",
                                                                                                                                                        "spearman" = "spearman")),
                                                    actionButton(inputId = "scanpy_find_clusters_button", "Find Clusters")
                                              )
                                       )
                                     )
                              ),
                              column(8,
                                     fluidRow(
                                       column(12,
                                              hidden(
                                                tags$div(class = "scanpy_clustering_plots", tabsetPanel(id = "scanpyClusteringPlotTabset", type = "tabs"
                                                ))
                                              )
                                       )
                                       
                                     )
                              )
                            ),
                            style = "primary"),
            bsCollapsePanel("Find Markers",
                            fluidRow(
                              column(4,
                                     fluidRow(
                                       column(12,
                                              panel(heading = "Options",
                                                    selectInput(
                                                      inputId = "scanpyFindMarkerSelectPhenotype",
                                                      label = "Select biological phenotype:",
                                                      choices = NULL
                                                    ),
                                                    conditionalPanel(
                                                      condition = "input.scanpyFindMarkerType == 'markerDiffExp'
                                                          || input.scanpyFindMarkerType == 'markerConserved'",
                                                      selectInput(
                                                        inputId = "scanpyFindMarkerGroup1",
                                                        label = "Select first group of interest:",
                                                        choices = NULL
                                                      ),
                                                      selectInput(
                                                        inputId = "scanpyFindMarkerGroup2",
                                                        label = "Select second group of interest:",
                                                        choices = NULL
                                                      )
                                                    ),
                                                    selectInput(
                                                      inputId = "scanpyFindMarkerTest",
                                                      label = "Select test:",
                                                      choices = c("t-test", "wilcoxon", "t-test_overestim_var", "logreg")
                                                    ),
                                                    selectInput(
                                                      inputId = "scanpyFindMarkerCorrMethod",
                                                      label = "Multiple testing correction method:",
                                                      choices = c("benjamini-hochberg", "bonferroni")
                                                    ),
                                                    numericInput(
                                                      inputId = "scanpyFindMarkerNGenes",
                                                      label = "No. of genes to return per group:", 
                                                      value = 1000, 
                                                      step = 1
                                                    ),
                                                    actionButton(inputId = "scanpyFindMarkerRun", "Find Markers")
                                              )
                                       )
                                     )
                              ),
                              column(8,
                                     fluidRow(
                                       column(12,
                                              hidden(
                                                tags$div(
                                                  class = "scanpy_findmarker_table",
                                                  filterTableUI(id = "filterScanpyFindMarker")
                                                )
                                              ),
                                              br(),
                                              hidden(
                                                tags$div(class = "scanpy_findmarker_jointHeatmap",
                                                         bsCollapse(
                                                           bsCollapsePanel(
                                                             title = "scanpy Heatmap Plot",
                                                             fluidRow(
                                                               column(12, align = "center",
                                                                      panel(
                                                                        numericInput("scanpy_findMarkerHeatmapPlotFullNumeric", value = 2, max = 2000, min = 2, step = 1, label = "Select number of top genes from each cluster/group to visualize in the heatmap below based on highest average log fold change value:"),
                                                                        actionButton("scanpy_findMarkerHeatmapPlotFullNumericRun", label = "Plot"),
                                                                        hr(),
                                                                        shinyjqui::jqui_resizable(
                                                                          plotOutput(outputId = "scanpy_findMarkerHeatmapPlotFull", height = "500px")
                                                                        )
                                                                      )
                                                               )
                                                             )
                                                           )
                                                         )
                                                )
                                              ),
                                              br(),
                                              hidden(
                                                tags$div(class = "scanpy_findmarker_plots",
                                                         panel(heading = "Marker Gene Plots",
                                                               HTML("<center><h5><span style='color:red; font-weight:bold; text-align:center;'>Click on the rows of the table above to plot the selected marker genes below!</span></h5></br></center>"),
                                                               tabsetPanel(id = "scanpyFindMarkerPlotTabset", type = "tabs"))
                                                )
                                              )
                                       )
                                       
                                     )
                              )
                            ),
                            style = "primary")
       ),
    nonLinearWorkflowUI(id = "nlw-scanpy") #nlw-seurat
    )