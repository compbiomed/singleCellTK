shinyPanelQC <- fluidPage(
  #useShinyalert(),
  tags$div(
    class = "container",
    wellPanel(
      sidebarLayout(
        sidebarPanel(
          h5(tags$a(href = paste0(docs.artPath, "ui_qc.html"),
                    "(help)", target = "_blank")),
          fluidRow(
            column(12, h3("Choose which algorithms to run:"))
          ),

          h4("General"),
          # QCMetrics
          checkboxInput("QCMetrics", "QC Metrics (Number of UMIs, number of features detected, etc.)"),
          shinyjs::hidden(
            tags$style(HTML("#QCMetricsParams {margin-left:40px}")),
            tags$div(id = "QCMetricsParams",
                     actionLink("QCMhelp", "Help", icon = icon("info-circle")),
                     tags$hr(),
                     selectInput("QCMgeneSets",
                                 "collectionName - Select a Gene Set for Quality Control",
                                 c("None" = "none",
                                   "Human Mitochondrial Genes (Ensembl)" = "he",
                                   "Human Mitochondrial Genes (Symbol)" = "hs",
                                   "Mouse Mitochondrial Genes (Ensembl)" = "me",
                                   "Mouse Mitochondrial Genes (Symbol)" = "ms")),
                     actionLink("QCImportGS", "Import Gene Sets", icon = icon("upload"))
            )
          ),
          tags$hr(),
          h4("Contamination Estimation"),
          # decontX
          checkboxInput("decontX", "decontX"),
          shinyjs::hidden(
            tags$style(HTML("#decontXParams {margin-left:40px}")),
            tags$div(id = "decontXParams",
                     actionLink("DXhelp", "Help", icon = icon("info-circle")),
                     tags$hr(),
                     numericInput("DXmaxIter", "maxIter - Maximum iterations of the EM algorithm (default 500)", 500),
                     numericInput("DXnativePrior", "nativePrior - Prior for native counts (default 10)", 10),
                     numericInput("DXcontPrior", "contaminationPrior - Prior for contamination counts (default 10)", 10),
                     numericInput("DXconvergence", "convergence - Threshold difference between previous and current iterations (default 0.001)", 0.001),
                     numericInput("DXiterLogLik", "iterLogLik - Number of iterations after which to calculate the log likelihood (default 10)", 10),
                     numericInput("DXvarGenes", "varGenes - Number of variable genes to use in dimensionality reduction before clustering (default 5000)", 5000),
                     numericInput("DXdbscanEps", "dbscanEps - Clustering resolution parameter (if no cell cluster labels) (default 1)", 1),

                     checkboxInput("DXestDelta", "estimateDelta - Estimate delta?"), # T/F input
                     checkboxInput("DXverbose", "verbose - Print log messages?", value = TRUE), # T/F input
            )
          ),
          checkboxInput("soupX", "SoupX"),
          shinyjs::hidden(
            tags$style(HTML("#soupXParams {margin-left:40px}")),
            tags$div(id = "soupXParams",
                     actionLink("SoupXhelp", "Help", icon = icon("info-circle")),
                     tags$hr(),
                     selectInput("soupXCluster", "cluster - Prior knowledge of clustering labels on cells (default None)", list()),
                     numericInput("soupXTfidfMin", "tfidfMin - Minimum value of tfidf to accept for a marker gene (default 1)", 1),
                     numericInput("soupXQuantile", "soupQuantile - Only use genes that are at or above this expression quantile in the soup (default 0.9)", 0.9, min = 0, max = 1),
                     numericInput("soupXMaxMarkers", "maxMarkers - If we have heaps of good markers, keep only the best maxMarkers of them. (Default 100)", 100, min = 1),
                     p("contaminationRange - This constrains the contamination fraction to lie within this range. (default 0.01 - 0.8)"),
                     numericInput("soupXContRangeLow", "Lower range:", 0.01, min = 0, max = 1),
                     numericInput("soupXContRangeHigh", "Higher range:", 0.8, min = 0, max = 1),
                     numericInput("soupXRhoMaxFDR", "rhoMaxFDR - FDR passed to SoupX::estimateNonExpressingCells, to test if rho is less than maximumContamination (default 0.2)", 0.2, min = 0, max = 1),
                     numericInput("soupXPriorRho", "priorRho - Mode of gamma distribution prior on contamination fraction (default 0.05)", 0.05, 0, 1),
                     numericInput("soupXPriorRhoStdDev", "priorRhoStdDev - Standard deviation of gamma distribution prior on contamination fraction (default 0.1)", 0.1),
                     checkboxInput("soupXForceAccept", "forceAccept - Should we allow very high contamination fractions to be used?", FALSE),
                     selectInput("soupXAdjustMethod", "AdjustMethod - Method to use for correction (default 'subtraction')",
                                 choices = c('subtraction', 'soupOnly', 'multinomial'), selected = "subtraction"),
                     checkboxInput("soupXRoundToInt", "roundToInt - Should the resulting matrix be rounded to integers?", FALSE),
                     numericInput("soupXTol", "tol - Allowed deviation from expected number of soup counts (default 0.001)", 0.001),
                     numericInput("soupXPCut", "pCut - The p-value cut-off used when method = 'soupOnly' (default 0.01)", 0.01, 0, 1)
            )
          ),
          tags$hr(),
          h4("Doublet Detection"),
          # scDblFinder
          checkboxInput("scDblFinder", "scDblFinder"),
          shinyjs::hidden(
            tags$style(HTML("#scDblFinderParams {margin-left:40px}")),
            tags$div(id = "scDblFinderParams",
                     actionLink("DChelp", "Help", icon = icon("info-circle")),
                     tags$hr(),
                     numericInput("DCnNeighbors", "nNeighbors - Number of nearest neighbors (default 50)", 50),
                     numericInput("DCsimDoublets", "simDoublets - Number of simulated doublets (default 10000)", 10000)
            )
          ),
          # cxds
          checkboxInput("cxds", "cxds"),
          shinyjs::hidden(
            tags$style(HTML("#cxdsParams {margin-left:40px}")),
            tags$div(id = "cxdsParams",
                     actionLink("CXhelp", "Help", icon = icon("info-circle")),
                     tags$hr(),
                     numericInput("CXntop", "ntop - Number of top variance genes (default 500)", 500),
                     numericInput("CXbinThresh", "binThresh - Threshold to consider a gene 'present' (default 0)", 0),

                     checkboxInput("CXverb", "verb - Output progress messages?", value = TRUE), # T/F input
                     checkboxInput("CXretRes", "retRes - Return gene pair scores and top-scoring gene pairs?"), # T/F input
                     checkboxInput("CXestNdbl", "estNdbl - Estimate the number of doublets?"), # T/F input
            )
          ),
          # bcds
          checkboxInput("bcds", "bcds"),
          shinyjs::hidden(
            tags$style(HTML("#bcdsParams {margin-left:40px}")),
            tags$div(id = "bcdsParams",
                     actionLink("BChelp", "Help", icon = icon("info-circle")),
                     tags$hr(),
                     numericInput("BCntop", "ntop - Number of top variance genes (default 500)", 500),
                     numericInput("BCsrat", "srat - Ratio between original number of cells and simulating doublets (decimal value, default 1)", 1),

                     textInput("BCnmax", "nmax - Max number of training rounds (default 'tune')", value = "tune"),

                     checkboxInput("BCverb", "verb - Output progress messages?", value = TRUE), # T/F input
                     checkboxInput("BCretRes", "retRes - Return trained classifier?"), # T/F input
                     checkboxInput("BCvarImp", "varImp - Return variable importance?"), # T/F input
                     checkboxInput("BCestNdbl", "estNdbl - Estimate the number of doublets?"), # T/F input
            )
          ),
          # cxds_bcds_hybrid
          checkboxInput("cxds_bcds_hybrid", "cxds_bcds_hybrid"),
          shinyjs::hidden(
            tags$style(HTML("#cxds_bcds_hybridParams {margin-left:40px}")),
            tags$div(id = "cxds_bcds_hybridParams",
                     actionLink("CXBChelp", "Help", icon = icon("info-circle")),
                     tags$hr(),
                     tags$label("cxds Parameters:"),
                     numericInput("CX2ntop", "ntop - Number of top variance genes (default 500)", 500),
                     numericInput("CX2binThresh", "binThresh- Threshold to consider a gene 'present' (default 0)", 0),
                     checkboxInput("CX2retRes", "retRes - Return gene pair scores and top-scoring gene pairs?"), # T/F input

                     tags$hr(),
                     tags$label("bcds Parameters:"),
                     numericInput("BC2ntop", "ntop - Number of top variance genes  (default 500)", 500),
                     numericInput("BC2srat", "srat - Ratio between original number of cells and simulating doublets (decimal value, default 1)", 1),
                     textInput("BC2nmax", "nmax - Max number of training rounds (default 'tune')", value = "tune"),
                     checkboxInput("BC2retRes", "retRes - Return trained classifier?"), # T/F input
                     checkboxInput("BC2varImp", "varImp - Return variable importance?"), # T/F input

                     checkboxInput("CXBCverb", "verb - Output bcds progress messages?", value = TRUE), # T/F input
                     checkboxInput("CXBCestNdbl", "estNdbl - Estimate the number of doublets?", value = TRUE), # T/F input
            )
          ),
          # scrublet
          checkboxInput("scrublet", "scrublet"),
          shinyjs::hidden(
            tags$style(HTML("#scrubletParams {margin-left:40px}")),
            tags$div(id = "scrubletParams",
                     actionLink("Shelp", "Help", icon = icon("info-circle")),
                     tags$hr(),
                     numericInput("SsimDoubletRatio", "simDoubletRatio - Number of soublets to simulate (default 2.0)", 2.0),
                     numericInput("SnNeighbors", "nNeighbors - Number of nearest neighbors (default NULL)", NULL),
                     numericInput("SminDist", "minDist - Tightness of UMAP points (default 0.1)", 0.1),
                     numericInput("SexpectedDoubletRate", "expectedDoubletRate - Estimated doublet rate (default 0.1)", 0.1),
                     numericInput("SstdevDoubletRate", "stdevDoubletRate - Uncertainty in expected doublet rate (default 0.02)", 0.02),
                     numericInput("SsyntheticDoubletUmiSubsampling", "syntheticDoubletUmiSubsampling - UMI sampling rate (default 0.1)", 0.1),
                     numericInput("SminCounts", "minCounts - Prior to PCA, exclude genes with counts below (defualt 3):", 3),
                     numericInput("SminCells", "minCells - Prior to PCA, exclude genes expressed in fewer than X cells (defualt 3):", 3),
                     numericInput("SminGeneVariabilityPctl", "minGeneVariabilityPctl - Prior to PCA, keep X most highly variable genes (defualt 3):", 3),
                     numericInput("SnPrinComps", "nPrinComps - Prior to KNN graph construction, number of principle components used to embed the transcriptomes (defualt 30):", 30),
                     numericInput("StsneAngle", "tsneAngle - Angular size of distant node as measured from a point in the t-SNE plot (defualt 0.5):", 0.5),
                     numericInput("StsnePerplexity", "tsnePerplexity - Number of nearest neighbors used in other manifold learning algorithms (defualt 30):", 30),

                     textInput("SdistanceMetric", "distanceMetric - Distance metric", value = "euclidean"),

                     checkboxInput("SuseApproxNeighbors", "useApproxNeighbors - Use approximate nearest neighbor method?", value = TRUE), # T/F input
                     checkboxInput("SgetDoubletNeighborParents", "getDoubletNeighborParents - Return doublet neighbors' parent transcriptomes?"), # T/F input
                     checkboxInput("SlogTransform", "logTransform - Log transform counts matrix?"), # T/F input
                     checkboxInput("SmeanCenter", "meanCenter - Center each gene's data at zero?", value = TRUE), # T/F input
                     checkboxInput("SnormalizeVariance", "normalizeVariance - Normalize each gene's data to have a variance of 1?", value = TRUE), # T/F input
                     checkboxInput("Sverbose", "verbose - Output progress updates?", value = TRUE), # T/F input
            )
          ),
          # doubletFinder
          checkboxInput("doubletFinder", "doubletFinder"),
          shinyjs::hidden(
            tags$style(HTML("#doubletFinderParams {margin-left:40px}")),
            tags$div(id = "doubletFinderParams",
                     actionLink("DFhelp", "Help", icon = icon("info-circle")),
                     tags$hr(),
                     numericInput("DFseuratNfeatures", "seuratNfeatures - Number of highly variable genes to use (default 2000)", 2000),
                     numericInput("DFseuratRes", "seuratRes - Seurat resolution (please enter comma-separated integers, default 1.5)", 1.5),
                     numericInput("DFformationRate", "formationRate - Doublet formation rate (default 0.075)", 0.075),
                     numericInput("DFseuratPcs", "seuratPcs - PCs to determine the number of clusters (default 15)", 15),

                     checkboxInput("DFverbose", "verbose - Output log messages?", value = TRUE), # T/F input
            )
          ),
          tags$hr(),
          h4("General Parameters"),
          selectizeInput(
            inputId = "qcAssaySelect", 
            label = "Select input matrix:", 
            choices = NULL, 
            selected = NULL, 
            multiple = FALSE,
            options = NULL),
          #uiOutput("qcAssaySelect"),
          #selectInput("qcAssaySelect", "Select an Assay", list()),
          selectInput("qcSampleSelect", "Select variable containing sample labels", list()),

          tags$hr(), # UMAP params
          h4("Quick UMAP Parameters"),
          textInput("QCUMAPName", "UMAP Name (default 'QC_UMAP')", value = "QC_UMAP"),
          numericInput("UnNeighbors", "Size of local neighborhood used for manifold approximation (default 30)", 30),
          numericInput("UnIterations", "Number of iterations performed during layout optimization (default 200)", 200),
          numericInput("Ualpha", 'Initial value of "learning rate" (default 1)', 1),
          numericInput("UminDist", "Effective minimum distance between embedded points (default 0.01)", 0.01),
          numericInput("Uspread", "Effective scale of embedded points (default 1)", 1),
          numericInput("UinitialDims", "Number of dimensions from PCA to use as input (default 25)", 25),
          numericInput("Useed", "Seed value for reproducibility of UMAP result (default 12345)", 12345),


          withBusyIndicatorUI(actionButton("runQC", "Run")),
          tags$div(id = "qcPageErrors"),
        ),
        mainPanel(
          tabsetPanel(
            id = "qcResPlotTabs"
          )
        )
      )
    )
  )
)

