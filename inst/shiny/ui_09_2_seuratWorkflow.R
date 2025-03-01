# User Interface for Seurat Workflow ---
shinyPanelSeurat <- fluidPage(
    tags$script("Shiny.addCustomMessageHandler('close_dropDownSeuratHM', function(x){
                  $('html').click();
                });"),
    h1("Seurat"),
    h5(tags$a(href = paste0(docs.artPath, "cnsl_seurat_curated_workflow.html"),
              "(help)", target = "_blank")),
    inlineCSS(list(".panel-danger>.panel-heading" = "background-color:#dcdcdc; color:#000000", ".panel-primary>.panel-heading" = "background-color:#f5f5f5; color:#000000; border-color:#dddddd", ".panel-primary" = "border-color:#dddddd;", ".panel-primary>.panel-heading+.panel-collapse>.panel-body" = "border-color:#dddddd;")),
    conditionalPanel(
        condition = "false",
        selectInput(
            "activePanelSelectSeurat",
            label = "Active Panel:",
            choices = c("",
                        "Normalize Data",
                        "Scale Data",
                        "Highly Variable Genes",
                        "Dimensionality Reduction",
                        "tSNE/UMAP",
                        "Clustering",
                        "Find Markers",
                        "Heatmap Plot"),
            selected = ""
        )
    ),
    bsCollapse(id = "SeuratUI", open = "Normalize Data",
            bsCollapsePanel("Normalize Data",
                fluidRow(
                    column(4,
                        panel(heading = "Options",
                              selectizeInput(
                                  inputId = "seuratSelectNormalizationAssay",
                                  label = "Select input matrix:",
                                  choices = NULL,
                                  selected = NULL,
                                  multiple = FALSE,
                                  options = NULL),
                            #uiOutput("seuratSelectNormalizationAssay"),
                            selectInput(inputId = "normalization_method", label = "Select normalization method: ", choices = c("LogNormalize", "CLR", "RC")),
                            textInput(inputId = "scale_factor", label = "Set scaling factor: ", value = "10000"),
                            actionButton(inputId = "normalize_button", "Normalize")
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
                                    selectInput(inputId = "hvg_method", label = "Select HVG method: ", choices = c("vst", "mean.var.plot", "dispersion")),
                                    textInput(inputId = "hvg_no_features", label = "Select number of features to find: ", value = "2000"),
                                    actionButton(inputId = "find_hvg_button", "Find HVG")
                                     )
                                  )
                                ),
                        br(),
                        fluidRow(
                            column(12,
                                panel(heading = "Display HVG",
                                    numericInput(inputId = "hvg_no_features_view", label = "Select number of features to display: ", value = 10, step = 1),
                                    verbatimTextOutput(outputId = "hvg_output", placeholder = TRUE)
                                     )
                                  )
                                )
                          ),
                     column(8,
                        fluidRow(
                            column(12,
                                panel(heading = "Plot",
                                    plotlyOutput(outputId = "plot_hvg")
                                     )
                                  )
                                )
                           )
                    ),
                    style = "primary"),

            # bsCollapsePanel("Scale Data",
            #                 fluidRow(
            #                   column(4,
            #                          panel(heading = "Options",
            #                                #selectInput(inputId = "model.use", label = "Select model for scaling: ", choices = c("linear", "poisson", "negbinom")),
            #                                materialSwitch(inputId = "do.scale", label = "Scale data?", value = TRUE),
            #                                materialSwitch(inputId = "do.center", label = "Center data?", value = TRUE),
            #                                textInput(inputId = "scale.max", label = "Max value for scaled data: ", value = "10"),
            #                                actionButton(inputId = "scale_button", "Scale")
            #                          )
            #                   )
            #                 ),
            #                 style = "primary"
            # ),

            bsCollapsePanel("Dimensionality Reduction",
                tabsetPanel(type = "tabs",
                    tabPanel("PCA",
                        br(),
                        fluidRow(
                            column(4,
                                fluidRow(
                                    column(12,
                                        panel(heading = "PCA",
                                            numericInput(inputId = "pca_no_components", label = "Select number of components to compute: ", value = 50),
                                            materialSwitch(inputId = "pca_compute_elbow", label = "Compute ElbowPlot?", value = TRUE),
                                            materialSwitch(inputId = "pca_compute_jackstraw", label = "Compute JackStrawPlot?", value = FALSE),
                                            materialSwitch(inputId = "pca_compute_heatmap", label = "Compute Heatmap?", value = TRUE),
                                            conditionalPanel(
                                              condition = 'input.pca_compute_heatmap == true',
                                              numericInput(inputId = "pca_compute_heatmap_nfeatures",
                                                           label = "Set number of features for heatmap:", value = 30, step = 1),
                                            ),
                                            numericInput(inputId = "seed_PCA",
                                                         label = "Seed value for reproducibility of result:",
                                                         value = 42,
                                                         step = 1),
                                            actionButton(inputId = "run_pca_button", "Run PCA")
                                             ),
                                        panel(heading = "Select No. of Components",
                                            htmlOutput(outputId = "pca_significant_pc_output", inline = FALSE),
                                            numericInput(inputId = "pca_significant_pc_counter", label = "Select number of components for downstream analysis: ", min = 1, max = 20, value = 10)
                                        )
                                          )
                                        )
                                  ),
                            column(8,
                                fluidRow(
                                    column(12,
                                           hidden(
                                        tags$div(class = "seurat_pca_plots", tabsetPanel(id = "seuratPCAPlotTabset", type = "tabs"
                                                   )
                                                 ))
                                          )

                                        )
                                  )
                               )

                            ),
                    tabPanel("ICA",
                             br(),
                             fluidRow(
                               column(4,
                                      fluidRow(
                                        column(12,
                                               panel(heading = "ICA",
                                                     textInput(inputId = "ica_no_components", label = "Select number of components to compute: ", value = "20"),
                                                     materialSwitch(inputId = "ica_compute_heatmap", label = "Compute Heatmap?", value = TRUE),
                                                     conditionalPanel(
                                                       condition = 'input.ica_compute_heatmap == true',
                                                       numericInput(inputId = "ica_compute_heatmap_nfeatures",
                                                                    label = "Set number of features for heatmap:", value = 30, step = 1),
                                                     ),
                                                     numericInput(inputId = "seed_ICA",
                                                                  label = "Seed value for reproducibility of result:",
                                                                  value = 42,
                                                                  step = 1),
                                                     actionButton(inputId = "run_ica_button", "Run ICA")
                                               ),
                                               panel(heading = "Select No. of Components",
                                                     #h5("Number of components suggested by ElbowPlot: "),
                                                     #verbatimTextOutput(outputId = "ica_significant_pc_output", placeholder = TRUE),
                                                     numericInput(inputId = "ica_significant_ic_counter", label = "Select number of components for downstream analysis: ", min = 1, max = 20, value = 10)
                                               )
                                        )
                                      )
                               ),
                               column(8,
                                      fluidRow(
                                        column(12,
                                               hidden(
                                                tags$div(class = "seurat_ica_plots", tabsetPanel(id="seuratICAPlotTabset", type = "tabs"
                                               ))
                                               )
                                               )
                                      )
                                      )
                             )

                             )
                    ),
                    style = "primary"),



            bsCollapsePanel("2D-Embedding",
                tabsetPanel(id = "tsneUmapTabsetSeurat", type = "tabs",
                            tabPanel("UMAP",
                                     br(),
                                     fluidRow(
                                       column(4,
                                              fluidRow(
                                                column(12,
                                                       panel(heading = "UMAP",
                                                             selectInput(inputId = "reduction_umap_method", label = "Select reduction method: ", choices = c("pca", "ica")),
                                                             #textInput(inputId = "reduction_umap_count", label = "Select number of reduction components: ", value = "20"),
                                                             numericInput(inputId = "min_dist_umap", label = "Set min.dist:", value = 0.3),
                                                             numericInput(inputId = "n_neighbors_umap", label = "Set n.neighbors:", value = 30, step = 1),
                                                             numericInput(inputId = "spread_umap", label = "Set spread:", value = 1),
                                                             numericInput(inputId = "seed_UMAP",
                                                                          label = "Seed value for reproducibility of result:",
                                                                          value = 42,
                                                                          step = 1),
                                                             htmlOutput(outputId = "display_message_umap", inline = FALSE),
                                                             actionButton(inputId = "run_umap_button", "Run UMAP")
                                                       )
                                                )
                                              )
                                       ),
                                       column(8,
                                              fluidRow(
                                                panel(heading = "Plot",
                                                      column(12,
                                                             plotlyOutput(outputId = "plot_umap")
                                                      )
                                                )
                                              )
                                       )
                                     )
                            ),
                    tabPanel("tSNE",
                        br(),
                        fluidRow(
                            column(4,
                                fluidRow(
                                    column(12,
                                        panel(heading = "tSNE",
                                            selectInput(inputId = "reduction_tsne_method", label = "Select reduction method: ", choices = c("pca", "ica")),
                                            #textInput(inputId = "reduction_tsne_count", label = "Select number of reduction components: ", value = "20"),
                                            numericInput(inputId = "perplexity_tsne", label = "Set perplexity:", value = 30),
                                            numericInput(inputId = "seed_TSNE",
                                                         label = "Seed value for reproducibility of result:",
                                                         value = 1,
                                                         step = 1),
                                            htmlOutput(outputId = "display_message_tsne", inline = FALSE),
                                            actionButton(inputId = "run_tsne_button", "Run tSNE")
                                             )
                                          )
                                        )
                                  ),
                            column(8,
                                fluidRow(
                                    panel(heading = "Plot",
                                        column(12,
                                            plotlyOutput(outputId = "plot_tsne")
                                              )
                                         )
                                        )
                                  )
                                )
                            )
                    ),
                    style = "primary"),

            bsCollapsePanel("Clustering",
                fluidRow(
                    column(4,
                        fluidRow(
                            column(12,
                                panel(heading = "Options",
                                    selectInput(inputId = "reduction_clustering_method", label = "Select reduction method: ", choices = c("pca", "ica")),
                                    #textInput(inputId = "reduction_clustering_count", label = "Select number of reduction components: ", value = "20"),
                                    selectInput(inputId = "algorithm.use", label = "Select clustering algorithm: ", choices = list("Original Louvain algorithm" = "louvain",
                                                                                                                               "Louvain algorithm with multilevel refinement" = "multilevel",
                                                                                                                               "SLM algorithm" = "SLM")),
                                    numericInput(inputId = "resolution_clustering", label = "Set resolution:", value = 0.8),
                                    materialSwitch(inputId = "group.singletons", label = "Group singletons?", value = TRUE),
                                    htmlOutput(outputId = "display_message_clustering", inline = FALSE),
                                    actionButton(inputId = "find_clusters_button", "Find Clusters")
                                    )
                                   )
                                )
                          ),
                    column(8,
                           fluidRow(
                             column(12,
                                    hidden(
                                    tags$div(class = "seurat_clustering_plots", tabsetPanel(id = "seuratClusteringPlotTabset", type = "tabs"
                                    ))
                                    )
                             )

                           )
                    )
                    ),
                    style = "primary"),
            bsCollapsePanel("Marker Gene Plots",
                            fluidRow(
                              column(12,
                                     fluidRow(
                                       column(12,
                                              panel(heading = "Options",
                                                    selectizeInput("selectGenesMarkerPlots", "Select a gene:", choices = NULL, multiple = TRUE, options = list('plugins' = list('remove_button'), server = TRUE)),
                                              )
                                       )
                                     )
                              ),
                              column(12,
                                     tags$div(class = "seurat_findmarker_plots",
                                              panel(heading = "",
                                                    tabsetPanel(id = "seuratFindMarkerPlotTabset", type = "tabs"))
                                     )
                              )
                            ),
                            style = "primary"
            ),
            bsCollapsePanel("Find Markers",
                            fluidRow(
                                column(4,
                                       fluidRow(
                                           column(12,
                                                  panel(heading = "Options",
                                                        h6("Compute marker genes that are either differentially expressed or conserved between selected groups and visualize them from the selected plots on right panel."),
                                                      radioButtons(
                                                          inputId = "seuratFindMarkerType",
                                                          label = "Select type of markers to identify:",
                                                          choices = c(
                                                              "markers between all groups" = "markerAll",
                                                              "markers differentially expressed between two selected groups" = "markerDiffExp",
                                                              "markers conserved between two selected groups" = "markerConserved"
                                                                      )
                                                      ),
                                                      selectInput(
                                                          inputId = "seuratFindMarkerSelectPhenotype",
                                                          label = "Select biological phenotype:",
                                                          choices = NULL
                                                      ),
                                                      conditionalPanel(
                                                          condition = "input.seuratFindMarkerType == 'markerDiffExp'
                                                          || input.seuratFindMarkerType == 'markerConserved'",
                                                          selectInput(
                                                              inputId = "seuratFindMarkerGroup1",
                                                              label = "Select first group of interest:",
                                                              choices = NULL
                                                          ),
                                                          selectInput(
                                                              inputId = "seuratFindMarkerGroup2",
                                                              label = "Select second group of interest:",
                                                              choices = NULL
                                                          )
                                                      ),
                                                          selectInput(
                                                              inputId = "seuratFindMarkerTest",
                                                              label = "Select test:",
                                                              choices = c("wilcox", "bimod",
                                                                          "t", "negbinom",
                                                                          "poisson", "LR",
                                                                          "DESeq2")
                                                          ),
                                                          selectInput(
                                                              inputId = "seuratFeatureUseReduction",
                                                              label = "Select reduction method for feature plot:",
                                                              choices = c("umap", "tsne",
                                                                          "pca", "ica")
                                                          ),
                                                          materialSwitch(
                                                              inputId = "seuratFindMarkerPosOnly",
                                                              label = "Only return positive markers?",
                                                              value = FALSE
                                                          ),
                                                          
                                                      actionButton(inputId = "seuratFindMarkerRun", "Find Markers")
                                                  )
                                           )
                                       )
                                ),
                                column(8,
                                       fluidRow(
                                           column(12,
                                                  hidden(
                                                      tags$div(
                                                          class = "seurat_findmarker_table",
                                                          filterTableUI(id = "filterSeuratFindMarker")
                                                      )
                                                  ),
                                                  br(),
                                                  hidden(
                                                      tags$div(class = "seurat_findmarker_jointHeatmap",
                                                               bsCollapse(
                                                                 bsCollapsePanel(
                                                                   title = "Heatmap Plot",
                                                                         fluidRow(
                                                                           column(12, align = "center",
                                                                                               panel(
                                                                                                   numericInput("findMarkerHeatmapPlotFullNumeric", value = 10, max = 2000, min = 2, step = 1, label = "Select number of top genes from each cluster/group to visualize in the heatmap below based on highest average log fold change value:"),
                                                                                                   actionButton("findMarkerHeatmapPlotFullNumericRun", label = "Plot"),
                                                                                                   hr(),
                                                                                                   shinyjqui::jqui_resizable(
                                                                                                       plotOutput(outputId = "findMarkerHeatmapPlotFull", height = "500px")
                                                                                                   )
                                                                                               )
                                                                           )
                                                                         )
                                                                 )
                                                               )
                                                      )
                                                  ),
                                           )

                                       )
                                )
                            ),
                            style = "primary")
       ),
    nonLinearWorkflowUI(id = "nlw-seurat")
    )

