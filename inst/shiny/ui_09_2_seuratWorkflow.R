# User Interface for Seurat Workflow ---
shinyPanelSeurat <- fluidPage(
    inlineCSS(list(".panel-danger>.panel-heading" = "background-color:#dcdcdc; color:#000000", ".panel-primary>.panel-heading" = "background-color:#f5f5f5; color:#000000; border-color:#dddddd", ".panel-primary" = "border-color:#dddddd;", ".panel-primary>.panel-heading+.panel-collapse>.panel-body" = "border-color:#dddddd;")),
        bsCollapse(id = "SeuratUI", open = "Data Input",
            bsCollapsePanel("Normalize Data",
                fluidRow(
                    column(4,
                        panel(
                            selectInput(inputId = "seuratSelectNormalizationAssay", label = "Select assay: ", choices = c()),
                            selectInput(inputId = "normalization_method", label = "Select normalization method: ", choices = c("LogNormalize", "CLR", "RC")),
                            textInput(inputId = "scale_factor", label = "Set scaling factor: ", value = "10000"),
                            actionButton(inputId = "normalize_button", "Normalize")
                             )
                          )
                        ),
                            style = "primary"
                            ),

            bsCollapsePanel("Scale Data",
                fluidRow(
                    column(4,
                        panel(
                            selectInput(inputId = "model.use", label = "Select model for scaling: ", choices = c("linear", "poisson", "negbinom")),
                            materialSwitch(inputId = "do.scale", label = "Scale data?", value = TRUE),
                            materialSwitch(inputId = "do.center", label = "Center data?", value = TRUE),
                            textInput(inputId = "scale.max", label = "Max value for scaled data: ", value = "10"),
                            actionButton(inputId = "scale_button", "Scale")
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
                                    textInput(inputId = "hvg_no_features_view", label = "Select number of features to display: ", value = "100"),
                                    verbatimTextOutput(outputId = "hvg_output", placeholder = TRUE)
                                     )
                                  )
                                )
                          ),
                     column(8,
                        fluidRow(
                            column(12,
                                panel(heading = "Plot",
                                    plotOutput(outputId = "plot_hvg")
                                     )
                                  )
                                )
                           )
                    ),
                    style = "primary"),

            bsCollapsePanel("Dimensionality Reduction",
                tabsetPanel(type = "tabs",
                    tabPanel("PCA",
                        br(),
                        fluidRow(
                            column(4,
                                fluidRow(
                                    column(12,
                                        panel(heading = "PCA",
                                            textInput(inputId = "pca_no_components", label = "Select number of components to compute: ", value = "20"),
                                            materialSwitch(inputId = "pca_compute_elbow", label = "Compute ElbowPlot?", value = TRUE),
                                            materialSwitch(inputId = "pca_compute_jackstraw", label = "Compute JackStrawPlot?", value = TRUE),
                                            materialSwitch(inputId = "pca_compute_heatmap", label = "Compute Heatmap?", value = TRUE),
                                            conditionalPanel(
                                              condition = 'input.pca_compute_heatmap == true',
                                              numericInput(inputId = "pca_compute_heatmap_nfeatures",
                                                           label = "Set number of features for heatmap:", value = 30, step = 1),
                                            ),
                                            actionButton(inputId = "run_pca_button", "Run PCA")
                                             ),
                                        panel(heading = "Select No. of Components",
                                            h5("Number of components suggested by ElbowPlot: "),
                                            verbatimTextOutput(outputId = "pca_significant_pc_output", placeholder = TRUE),
                                            sliderInput(inputId = "pca_significant_pc_slider", label = "Select number of components for downstream analysis: ", min = 1, max = 20, value = 10, round = TRUE)
                                        )
                                          )
                                        )
                                  ),
                            column(8,
                                fluidRow(
                                    column(12,
                                        tags$div(class = "seurat_pca_plots", tabsetPanel(type = "tabs",
                                            tabPanel(title = "PCA Plot",
                                                panel(heading = "PCA Plot",
                                                    plotOutput(outputId = "plot_pca")
                                                     )
                                                    ),
                                            tabPanel(title = "Elbow Plot",
                                                panel(heading = "Elbow Plot",
                                                    plotOutput(outputId = "plot_elbow_pca")
                                                     )
                                                    ),
                                            tabPanel(title = "JackStraw Plot",
                                                panel(heading = "JackStraw Plot",
                                                    plotOutput(outputId = "plot_jackstraw_pca")
                                                     )
                                                    ),
                                            tabPanel(title = "Heatmap Plot",
                                                panel(heading = "Heatmap Plot",
                                                    panel(heading = "Plot Options",
                                                        fluidRow(
                                                            column(6,
                                                                pickerInput(inputId = "picker_dimheatmap_components_pca", label = "Select principal components to plot:", choices = c(), options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), multiple = TRUE)
                                                            ),
                                                            column(6,
                                                                sliderInput(inputId = "slider_dimheatmap_pca", label = "Number of columns for the plot: ", min = 1, max = 4, value = 2)
                                                            )
                                                        ),
                                                        actionButton(inputId = "plot_heatmap_pca_button", "Plot")
                                                        ),
                                                    panel(heading = "Plot",
                                                        jqui_resizable(plotOutput(outputId = "plot_heatmap_pca"), options = list(maxWidth = 700))
                                                        )
                                                     )
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
                                                     actionButton(inputId = "run_ica_button", "Run ICA")
                                               ),
                                               panel(heading = "Select No. of Components",
                                                     #h5("Number of components suggested by ElbowPlot: "),
                                                     #verbatimTextOutput(outputId = "ica_significant_pc_output", placeholder = TRUE),
                                                     sliderInput(inputId = "ica_significant_ic_slider", label = "Select number of components for downstream analysis: ", min = 1, max = 20, value = 10, round = TRUE)
                                               )
                                        )
                                      )
                               ),
                               column(8,
                                      fluidRow(
                                        column(12,
                                                tags$div(class = "seurat_ica_plots", tabsetPanel(type = "tabs",
                                                           tabPanel(title = "ICA Plot",
                                                                    panel(heading = "ICA Plot",
                                                                          plotOutput(outputId = "plot_ica")
                                                                    )
                                                           ),
                                                           tabPanel(title = "Heatmap Plot",
                                                                    panel(heading = "Heatmap Plot",
                                                                          panel(heading = "Plot Options",
                                                                                fluidRow(
                                                                                  column(6,
                                                                                         pickerInput(inputId = "picker_dimheatmap_components_ica", label = "Select principal components to plot:", choices = c(), options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), multiple = TRUE)
                                                                                  ),
                                                                                  column(6,
                                                                                         sliderInput(inputId = "slider_dimheatmap_ica", label = "Number of columns for the plot: ", min = 1, max = 4, value = 2)
                                                                                  )
                                                                                ),
                                                                                actionButton(inputId = "plot_heatmap_ica_button", "Plot")
                                                                          ),
                                                                          panel(heading = "Plot",
                                                                                jqui_resizable(plotOutput(outputId = "plot_heatmap_ica"), options = list(maxWidth = 700))
                                                                          )
                                                                    )
                                                           )
                                               ))
                                               )
                                      )
                                      )
                             )
                             
                             )
                    ),
                    style = "primary"),



            bsCollapsePanel("tSNE/UMAP",
                tabsetPanel(type = "tabs",
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
                                            plotOutput(outputId = "plot_tsne")
                                              )
                                         )
                                        )
                                  )
                                )
                            ),
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
                                            plotOutput(outputId = "plot_umap")
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
                                panel(
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
                                    tags$div(class = "seurat_clustering_plots", tabsetPanel(type = "tabs",
                                                tabPanel(title = "PCA Plot",
                                                         panel(heading = "PCA Plot",
                                                               plotOutput(outputId = "plot_pca_clustering")
                                                         )
                                                ),
                                                tabPanel(title = "ICA Plot",
                                                         panel(heading = "ICA Plot",
                                                               plotOutput(outputId = "plot_ica_clustering")
                                                         )
                                                ),
                                                tabPanel(title = "tSNE Plot",
                                                         panel(heading = "tSNE Plot",
                                                               plotOutput(outputId = "plot_tsne_clustering")
                                                         )
                                                ),
                                                tabPanel(title = "UMAP Plot",
                                                         panel(heading = "UMAP Plot",
                                                               plotOutput(outputId = "plot_umap_clustering")
                                                         )
                                                )
                                    ))
                             )
                             
                           )
                    )
                    ),
                    style = "primary")
       )
    )
# ----
