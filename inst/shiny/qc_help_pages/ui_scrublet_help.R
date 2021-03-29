scrubletHelpModal <- function() {
  modalDialog(
    tags$style(HTML("
      div {
        word-wrap: break-word;
      }
      ")),
    tags$div(
      h3("Srublet Parameters"),
      fluidRow(
        column(4, tags$b("Parameter Name")),
        column(8, tags$b("Description"))
      ),
      tags$hr(),
      fluidRow(
        column(4, "simDoubletRatio"),
        column(8, "Numeric. Number of doublets to simulate relative to the number of observed transcriptomes. Default 2.0.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "nNeighbors"),
        column(8, 'Integer. Number of neighbors used to construct the KNN graph of observed transcriptomes and simulated 
               doublets. If NULL, this is set to round(0.5 * sqrt(n_cells)). Default NULL.')
      ),
      tags$hr(),
      fluidRow(
        column(4, "minDist"),
        column(8, "Float Determines how tightly UMAP packs points together. If NULL, this is set to 0.1. Default NULL.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "expectedDoubletRate"),
        column(8, "The estimated doublet rate for the experiment. Default 0.1.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "stdevDoubletRate"),
        column(8, "Uncertainty in the expected doublet rate. Default 0.02.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "syntheticDoubletUmiSubsampling"),
        column(8, "Numeric. Rate for sampling UMIs when creating synthetic doublets. If 1.0, each doublet is created by simply 
               adding the UMIs from two randomly sampled observed transcriptomes. For values less than 1, the UMI counts are 
               added and then randomly sampled at the specified rate. Defuault: 1.0.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "minCounts"),
        column(8, "Numeric. Used for gene filtering prior to PCA. Genes expressed at fewer than minCounts in fewer than minCells 
               are excluded. Default 3.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "minCells"),
        column(8, "Integer. Used for gene filtering prior to PCA. Genes expressed at fewer than minCounts in fewer than minCells 
               are excluded. Default 3.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "minGeneVariabilityPctl"),
        column(8, "Numeric. Used for gene filtering prior to PCA. Keep the most highly variable genes (in the top minGeneVariabilityPctl 
               percentile), as measured by the v-statistic (Klein et al., Cell 2015). Default 85.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "nPrinComps"),
        column(8, "Integer. Number of principal components used to embed the transcriptomes prior to k-nearest-neighbor graph
               construction. Default 30.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "tsneAngle"),
        column(8, "Float. Determines angular size of a distant node as measured from a point in the t-SNE plot. If default, it is set 
               to 0.5 Default NULL.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "tsnePerplexity"),
        column(8, "Integer. The number of nearest neighbors that is used in other manifold learning algorithms. If default, it is set
               to 30. Default NULL.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "distanceMetric"),
        column(8, 'Character. Distance metric used when finding nearest neighbors. For list of valid values, see the documentation for 
               annoy (if useApproxNeighbors is TRUE) or sklearn.neighbors.NearestNeighbors (if useApproxNeighbors is FALSE). 
               Default "euclidean".')
      ),
      tags$hr(),
      fluidRow(
        column(4, "useApproxNeighbors"),
        column(8, "Check off to use approximate nearest neighbor method (annoy) for the KNN classifier. Default TRUE.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "getDoubletNeighborParents"),
        column(8, "Check off to return the parent transcriptomes that generated the doublet neighbors of each observed transcriptome. 
               This information can be used to infer the cell states that generated a given doublet state. Default FALSE.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "logTransform"),
        column(8, "Check off to log-transform the counts matrix (log10(1+TPM)). sklearn.decomposition.TruncatedSVD will be used for 
               dimensionality reduction, unless meanCenter is TRUE. Default FALSE.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "meanCenter"),
        column(8, "Check off to center the data such that each gene has a mean of 0. sklearn.decomposition.PCA will be used for 
               dimensionality reduction. Default TRUE.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "normalizeVariance"),
        column(8, "Check off to normalize the data such that each gene has a variance of 1. sklearn.decomposition.TruncatedSVD will be 
               used for dimensionality reduction, unless meanCenter is TRUE. Default TRUE.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "verbose"),
        column(8, "Check off to print progress updates to the console. Default TRUE.")
      ),
    )
  )
}

