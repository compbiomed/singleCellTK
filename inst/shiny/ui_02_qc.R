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
                     selectInput("QCMgeneSets", "Select a Gene Set for Quality Control", c("None")),
            )
          ),
          tags$hr(),
          h4("Contamination Estimation"),
          # decontX
          checkboxInput("decontX", "decontX"),
          shinyjs::hidden(
            tags$style(HTML("#decontXParams {margin-left:40px}")),
            tags$div(id = "decontXParams",
                     numericInput("DXmaxIter", "Maximum iterations of the EM algorithm (default 500)", 500),
                     numericInput("DXnativePrior", "Prior for native counts (default 10)", 10),
                     numericInput("DXcontPrior", "Prior for contamination counts (default 10)", 10),
                     numericInput("DXconvergence", "Threshold difference between previous and current iterations (default 0.001)", 0.001),
                     numericInput("DXiterLogLik", "Number of iterations after which to calculate the log likelihood (default 10)", 10),
                     numericInput("DXvarGenes", "Number of variable genes to use in dimensionality reduction before clustering (default 5000)", 5000),
                     numericInput("DXdbscanEps", "Clustering resolution parameter (if no cell cluster labels) (default 1)", 1),
                     
                     textInput("DXz", "Cell cluster labels (please enter space-separated values) (default NULL)"),
                     textInput("DXbatch", "Batch labels (please enter space-separated values) (default NULL)"),
                     
                     checkboxInput("DXestDelta", "Estimate delta?"), # T/F input
                     checkboxInput("DXverbose", "Print log messages?"), # T/F input
                     
                     fileInput('DXlogfile', 'Choose a file to write log messages to (otherwise messages will be printed to stdout)')
            )
          ),
          tags$hr(),
          h4("Doublet Detection"),
          # doubletCells
          checkboxInput("doubletCells", "doubletCells"),
          shinyjs::hidden(
            tags$style(HTML("#doubletCellsParams {margin-left:40px}")),
            tags$div(id = "doubletCellsParams",
                     numericInput("DCnNeighbors", "Number of nearest neighbors (default 50)", 50),
                     numericInput("DCsimDoublets", "Number of simulated doublets (default 10000)", 10000)
            )
          ),
          # cxds
          checkboxInput("cxds", "cxds"),
          shinyjs::hidden(
            tags$style(HTML("#cxdsParams {margin-left:40px}")),
            tags$div(id = "cxdsParams",
                     numericInput("CXntop", "Number of top variance genes (default 500)", 500),
                     numericInput("CXbinThresh", "Threshold to consider a gene 'present' (default 0)", 0),
                     
                     checkboxInput("CXverb", "Output progress messages?"), # T/F input
                     checkboxInput("CXretRes", "Return gene pair scores and top-scoring gene pairs?"), # T/F input
            )
          ),
          # bcds
          checkboxInput("bcds", "bcds"),
          shinyjs::hidden(
            tags$style(HTML("#bcdsParams {margin-left:40px}")),
            tags$div(id = "bcdsParams",
                     numericInput("BCntop", "Number of top variance genes (default 500)", 500),
                     numericInput("BCsrat", "Ratio between original number of cells and simulating doublets (decimal value, default 1)", 1),
                     
                     textInput("BCnmax", "Max number of training rounds (default 'tune')", value = "tune"),
                     
                     checkboxInput("BCverb", "Output progress messages?"), # T/F input
                     checkboxInput("BCretRes", "Return trained classifier?"), # T/F input
                     checkboxInput("BCvarImp", "Return variable importance?"), # T/F input
            )
          ),
          # cxds_bcds_hybrid
          # TODO: come back to after confirming inputs "cxdsArgs" and "bcdsArgs"
          checkboxInput("cxds_bcds_hybrid", "cxds_bcds_hybrid"),
          shinyjs::hidden(
            tags$style(HTML("#cxds_bcds_hybridParams {margin-left:40px}")),
            tags$div(id = "cxds_bcds_hybridParams",
                     tags$label("cxds Parameters:"),
                     numericInput("CX2ntop", "Number of top variance genes (default 500)", 500),
                     numericInput("CX2binThresh", "Threshold to consider a gene 'present' (default 0)", 0),
                     checkboxInput("CX2retRes", "Return gene pair scores and top-scoring gene pairs?"), # T/F input
                    
                     tags$hr(),
                     tags$label("bcds Parameters:"),
                     numericInput("BC2ntop", "Number of top variance genes  (default 500)", 500),
                     numericInput("BC2srat", "Ratio between original number of cells and simulating doublets (decimal value, default 1)", 1),
                     textInput("BC2nmax", "Max number of training rounds (default 'tune')", value = "tune"),
                     
                     checkboxInput("BC2retRes", "Return trained classifier?"), # T/F input
                     checkboxInput("BC2varImp", "Return variable importance?"), # T/F input
                     checkboxInput("CXBCverb", "Output bcds progress messages?"), # T/F input
            )
          ),
          # scrublet
          # TODO: double check these inputs too
          checkboxInput("scrublet", "scrublet"),
          shinyjs::hidden(
            tags$style(HTML("#scrubletParams {margin-left:40px}")),
            tags$div(id = "scrubletParams",
                     numericInput("SsimDoubletRatio", "Number of soublets to simulate (default 2.0)", 2.0),
                     numericInput("SnNeighbors", "Number of nearest neighbors (default NULL)", NULL),
                     numericInput("SminDist", "Tightness of UMAP points (default 0.1)", 0.1),
                     numericInput("SexpectedDoubletRate", "Estimated doublet rate (default 0.1)", 0.1),
                     numericInput("SstdevDoubletRate", "Uncertainty in expected doublet rate (default 0.02)", 0.02),
                     numericInput("SsyntheticDoubletUmiSubsampling", "UMI sampling rate (default 0.1)", 0.1),
                     numericInput("SminCounts", "Prior to PCA, exclude genes with counts below (defualt 3):", 3),
                     numericInput("SminCells", "Prior to PCA, exclude genes expressed in fewer than X cells (defualt 3):", 3),
                     numericInput("SminGeneVariabilityPctl", "Prior to PCA, keep X most highly variable genes (defualt 3):", 3),
                     numericInput("SnPrinComps", "Prior to KNN graph construction, number of principle components used to embed the transcriptomes (defualt 30):", 30),
                     numericInput("StsneAngle", "Angular size of distant node as measured from a point in the t-SNE plot (defualt 0.5):", 0.5),
                     numericInput("StsnePerplexity", "Number of nearest neighbors used in other manifold learning algorithms (defualt 30):", 30),
                     
                     textInput("SdistanceMetric", "Distance metric", value = "euclidean"),

                     checkboxInput("SuseApproxNeighbors", "Use approximate nearest neighbor method?"), # T/F input
                     checkboxInput("SgetDoubletNeighborParents", "Return doublet neighbors' parent transcriptomes?"), # T/F input
                     checkboxInput("SlogTransform", "Log transform counts matrix?"), # T/F input
                     checkboxInput("SmeanCenter", "Center each gene's data at zero?"), # T/F input
                     checkboxInput("SnormalizeVariance", "Normalize each gene's data to have a variance of 1?"), # T/F input
                     checkboxInput("Sverbose", "Output progress updates?"), # T/F input
            )
          ),
          # doubletFinder
          checkboxInput("doubletFinder", "doubletFinder"),
          shinyjs::hidden(
            tags$style(HTML("#doubletFinderParams {margin-left:40px}")),
            tags$div(id = "doubletFinderParams",
                     numericInput("DFseuratNfeatures", "Number of highly variable genes to use (default 2000)", 2000),
                     numericInput("DFseuratRes", "Seurat resolution (please enter comma-separated integers, default 1.5)", 1.5),
                     numericInput("DFformationRate", "Doublet formation rate (default 0.075)", 0.075),
                     numericInput("DFseuratPcs", "PCs to determine the number of clusters (default 15)", 15),
                     
                     checkboxInput("DFverbose", "Output log messages?"), # T/F input
            )
          ),
          tags$hr(),
          h4("General Paramters"),
          selectInput("qcAssaySelect", "Select an Assay", list()),
          selectInput("qcSampleSelect", "Select a Sample", list()),
          
          withBusyIndicatorUI(actionButton("runQC", "Run")),
          tags$div(id = "qcPageErrors"),
          
          shinyjs::hidden(
            tags$div(id = "qcPlotSection",
                     tags$hr(), # start plot subsection
                     h4("Plot Parameters"),
                     selectInput("qcPlotRedDim", "Select an ReducedDim obejct", list()),
                     withBusyIndicatorUI(actionButton("plotQC", "Plot")),
            )
          ),
          
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