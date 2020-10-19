shinyPanelQC <- fluidPage(
  useShinyalert(),
  tags$div(
    class = "container",
    h1("Data QC"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v04-tab02_Data-Summary-and-Filtering.html",
              "(help)", target = "_blank")),
    wellPanel(
      sidebarLayout(
        sidebarPanel(
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
                     selectInput("QCMgeneSets", "collectionName - Select a Gene Set for Quality Control", c("None")),
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
          tags$hr(),
          h4("Doublet Detection"),
          # doubletCells
          checkboxInput("doubletCells", "doubletCells"),
          shinyjs::hidden(
            tags$style(HTML("#doubletCellsParams {margin-left:40px}")),
            tags$div(id = "doubletCellsParams",
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
          h4("General Paramters"),
          selectInput("qcAssaySelect", "Select an Assay", list()),
          selectInput("qcSampleSelect", "Select variable containing sample labels", list()),
          
          withBusyIndicatorUI(actionButton("runQC", "Run")),
          tags$div(id = "qcPageErrors"),
        ),
        mainPanel(
          shinyjs::hidden(
            tags$div(id = "qcPlotSection",
                     tags$hr(), # start plot subsection
                     h4("Plot Parameters"),
                     selectInput("qcPlotRedDim", "Select an ReducedDim obejct", list()),
                     # withBusyIndicatorUI(actionButton("plotQC", "Plot")),
            )
          ),
          tabsetPanel(
            id = "qcResPlotTabs"
          )
        )
      )
    )
  )
)

