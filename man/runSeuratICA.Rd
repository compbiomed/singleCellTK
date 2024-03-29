% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seuratFunctions.R
\name{runSeuratICA}
\alias{runSeuratICA}
\title{runSeuratICA
Computes ICA on the input sce object and stores the calculated independent
components within the sce object}
\usage{
runSeuratICA(
  inSCE,
  useAssay = "seuratScaledData",
  useFeatureSubset = NULL,
  scale = TRUE,
  reducedDimName = "seuratICA",
  nics = 20,
  seed = 12345,
  verbose = FALSE
)
}
\arguments{
\item{inSCE}{(sce) object on which to compute ICA}

\item{useAssay}{Assay containing scaled counts to use in ICA.}

\item{useFeatureSubset}{Subset of feature to use for dimension reduction. A
character string indicating a \code{rowData} variable that stores the logical
vector of HVG selection, or a vector that can subset the rows of
\code{inSCE}. Default \code{NULL}.}

\item{scale}{Logical scalar, whether to standardize the expression values
using \code{\link[Seurat]{ScaleData}}. Default \code{TRUE}.}

\item{reducedDimName}{Name of new reducedDims object containing Seurat ICA
Default \code{seuratICA}.}

\item{nics}{Number of independent components to compute. Default \code{20}.}

\item{seed}{Random seed for reproducibility of results.
Default \code{NULL} will use global seed in use by the R environment.}

\item{verbose}{Logical value indicating if informative messages should
be displayed. Default is \code{TRUE}.}
}
\value{
Updated \code{SingleCellExperiment} object which now contains the
computed independent components
}
\description{
runSeuratICA
Computes ICA on the input sce object and stores the calculated independent
components within the sce object
}
\details{
For features used for computation, it can be controlled by \code{features} or
\code{useFeatureSubset}. When \code{features} is specified, the scaling and
dimensionality reduction will only be processed with these features. When
\code{features} is \code{NULL} but \code{useFeatureSubset} is specified, will
use the features that the HVG list points to. If both parameters are
\code{NULL}, the function will see if any Seurat's variable feature detection
has been ever performed, and use them if found. Otherwise, all features are
used.
}
\examples{
data(scExample, package = "singleCellTK")
\dontrun{
sce <- runSeuratNormalizeData(sce, useAssay = "counts")
sce <- runSeuratFindHVG(sce, useAssay = "counts")
sce <- runSeuratScaleData(sce, useAssay = "counts")
sce <- runSeuratICA(sce, useAssay = "counts")
}
}
