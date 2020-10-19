.constructSCE <- function(
  matrices,
  features,
  barcodes,
  metadata,
  reducedDims) {

  sce <- SingleCellExperiment::SingleCellExperiment(assays = matrices)
  SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(features)
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(barcodes)
  S4Vectors::metadata(sce) <- metadata
  SingleCellExperiment::reducedDims(sce) <- reducedDims
  return(sce)
}

.getDimUnion <- function(dataList){
  Row <- lapply(dataList, function(x) {rownames(x)})
  RowUnion <- base::Reduce(union, Row)
  Col <- lapply(dataList, function(x) {colnames(x)})
  ColUnion <- base::Reduce(union, Col)
  return(list(RowUnion, ColUnion))
}

.getMatUnion <- function(dimsList, x,
                         combineRow, combineCol,
                         sparse = FALSE,
                         fill = c("NA", "0")){
  row <- dimsList[[1]]
  col <- dimsList[[2]]
  matOrigin <- x
  fill <- match.arg(fill)
  if (fill == "0") {
    fill <- 0
  } else {
    fill <- NA
  }

  ### combine row
  if (isTRUE(combineRow) & (!is.null(row))) {
    missRow <- row[!row %in% rownames(x)]
    missMat <- Matrix::Matrix(fill, nrow = length(missRow), ncol = ncol(matOrigin),
                            dimnames = list(missRow, colnames(matOrigin)))
    if (!isTRUE(sparse)) {
      missMat <- as.matrix(missMat)
    }

    mat <- rbind(matOrigin, missMat)
    if (anyDuplicated(rownames(mat))) {
      mat <- mat[!duplicated(rownames(mat)), ]
    }
    matOrigin <- mat[row, ]
  }

  ### combine cols
  if (isTRUE(combineCol) & (!is.null(col))) {
    missCol <- col[!col %in% colnames(x)]
    missMat <- Matrix::Matrix(fill, nrow = nrow(matOrigin), ncol = length(missCol),
                            dimnames = list(rownames(matOrigin), missCol))
    if (!isTRUE(sparse)) {
      missMat <- as.matrix(missMat)
    }

    mat <- cbind(matOrigin, missMat)
    if (anyDuplicated(colnames(mat))) {
      mat <- mat[, !duplicated(colnames(mat))]
    }
    matOrigin <- mat[, col]
  }
  return(matOrigin)
}


.mergeRowDataSCE <- function(sceList, by.r) {
  feList <- lapply(sceList, function(x){
    rw <- SummarizedExperiment::rowData(x)
    rw[['rownames']] <- rownames(rw)
    return(rw)
  })

  ## Get merged rowData
  by.r <- unique(c('rownames', by.r))
  unionFe <- Reduce(function(r1, r2) merge(r1, r2, by=by.r, all=TRUE), feList)
  allGenes <- unique(unlist(lapply(feList, rownames)))

  ## rowData
  newFe <- unionFe
  if (nrow(newFe) != length(allGenes)) {
    warning("Conflicts were found when merging two rowData. ",
            "Resolved the conflicts by choosing the first entries.",
            "To avoid conflicts, please provide the 'by.r' arguments to ",
            "specify columns in rowData that does not have conflict between two singleCellExperiment object. ")
    newFe <- newFe[!duplicated(newFe$rownames), ]
  }
  rownames(newFe) <- newFe[['rownames']]
  newFe <- newFe[allGenes,]
  return(newFe)
}

.mergeColDataSCE <- function(sceList, by.c) {
  cbList <- lapply(sceList, function(x) {
    cD <- SummarizedExperiment::colData(x)
    cD[['rownames']] <- rownames(cD)
    return(cD)
  })

  by.c <- unique(c("rownames", by.c))
  unionCb <- Reduce(function(c1, c2) merge(c1, c2, by=by.c, all=TRUE), cbList)
  rownames(unionCb) <- unionCb[['rownames']]
  newCbList <- list()
  for (i in seq_along(sceList)) {
    newCbList[[i]] <- unionCb[colnames(sceList[[i]]),]
  }
  return(newCbList)
}

.mergeRedimSCE <- function(sceList, reduceList) {
  ## get reducedDims for each SCE SummarizedExperiment::
  reduceList <- lapply(sceList, SingleCellExperiment::reducedDims)
  ## get every reducedDim exists in at least one SCEs
  UnionReducedDims <- unique(unlist(lapply(sceList, SingleCellExperiment::reducedDimNames)))

  ## for each reducedDim, get union row/cols
  reducedDims <- list()
  for (reduceDim in UnionReducedDims) {
    x <- lapply(sceList, function(x) {if (reduceDim %in% SingleCellExperiment::reducedDimNames(x)) {SingleCellExperiment::reducedDim(x, reduceDim)}})
    reducedDims[[reduceDim]] <- .getDimUnion(x)
  }

  ## Merge reducedDim for each SCE
  redList <- list()
  for (idx in seq_along(sceList)){
    redMat <- reduceList[[idx]]

    for (DimName in UnionReducedDims) {
      if (DimName %in% names(redMat)) {
        redMat[[DimName]] <- .getMatUnion(reducedDims[[DimName]], redMat[[DimName]],
                                          combineRow = FALSE, combineCol = TRUE,
                                          sparse = FALSE, fill = "NA")
      } else {
        redMat[[DimName]] <- base::matrix(NA, nrow = ncol(sceList[[idx]]),
                                          ncol = length(reducedDims[[DimName]][[2]]),
                                          dimnames = list(colnames(sceList[[idx]]), reducedDims[[DimName]][[2]]))
      }
    }

    redList[[idx]] <- redMat
  }

  return(redList)
}

.mergeAssaySCE <- function(sceList) {
  UnionAssays <- Reduce(function(d1, d2) base::union(d1, d2),
                        lapply(sceList, SummarizedExperiment::assayNames))
  assayList <- lapply(sceList, assays)
  assayDims <- list(
    unique(unlist(lapply(sceList, rownames))),
    unique(unlist(lapply(sceList, colnames)))
  )

  asList <- list()
  for (idx in seq_along(assayList)){
    assay <- assayList[[idx]]
    for (assayName in UnionAssays) {
      if (assayName %in% names(assay)) {
        assay[[assayName]] <- .getMatUnion(assayDims, assay[[assayName]],
                                           combineRow = TRUE, combineCol = FALSE,
                                           sparse = TRUE, fill = "0")
      } else{
        assay[[assayName]] <- Matrix::Matrix(0, nrow = length(assayDims[[1]]),
                                             ncol = ncol(sceList[[idx]]),
                                             dimnames = list(assayDims[[1]], colnames(sceList[[idx]]))) #assayDims[[assayName]])
      }
    }
    asList[[idx]] <- assay
  }

  return(asList)
}

# .mergeMetaSCE <- function(sceList) {
#   metaList <- lapply(sceList, S4Vectors::metadata)
#   metaNames <- unlist(lapply(metaList, names))

#   if ("runBarcodeRanksMetaOutput" %in% metaNames) {
#     barcodeMetas <- lapply(metaList, function(x) {x[["runBarcodeRanksMetaOutput"]]})
#     barcodeMetas <- do.call(rbind, barcodeMetas)

#     for (i in seq_along(metaList)) {
#       metaList[[i]][["runBarcodeRanksMetaOutput"]] <- NULL
#     }

#     metaList[["runBarcodeRanksMetaOutput"]] <- barcodeMetas
#   }

#   return(metaList)
# }

.mergeMetaSCE <- function(SCE_list) {
  sampleMeta <- lapply(SCE_list, S4Vectors::metadata)
  metaNames <- unique(unlist(lapply(sampleMeta, names)))
  NewMeta <- list()

  for (meta in metaNames) {
    for (i in seq_along(sampleMeta)) {
      NewMeta[[meta]][[i]] <- sampleMeta[[i]][[meta]]
    }
  }

  if ("runBarcodeRanksMetaOutput" %in% metaNames) {
    NewMeta[["runBarcodeRanksMetaOutput"]] <- unlist(NewMeta[["runBarcodeRanksMetaOutput"]])
  }

  return(NewMeta)
}

#' Combine a list of SingleCellExperiment objects as one SingleCellExperiment object
#' @param sceList A list contains \link[SingleCellExperiment]{SingleCellExperiment} objects
#' @param by.r Specifications of the columns used for merging rowData. See 'Details'.
#' @param by.c Specifications of the columns used for merging colData. See 'Details'.
#' @param combined logical; if TRUE, it will combine the list of SingleCellExperiment objects. See 'Details'.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object which combines all
#' objects in sceList. The colData is merged.
#' @examples
#' combinedsce <- combineSCE(list(sce,sce), by.r = NULL, by.c = NULL, combined = TRUE)
#' @export

combineSCE <- function(sceList, by.r, by.c, combined){
  ##  rowData
  newFeList <- .mergeRowDataSCE(sceList, by.r)
  ## colData
  newCbList <- .mergeColDataSCE(sceList, by.c)
  ## reducedDim
  redMatList <- .mergeRedimSCE(sceList)
  ## assay
  assayList <- .mergeAssaySCE(sceList)

  New_SCE <- list()
  for (i in seq(length(sceList))) {
    ## create new sce
    New_SCE[[i]] <- .constructSCE(matrices = assayList[[i]], features = newFeList,
                                  barcodes = newCbList[[i]],
                                  metadata = S4Vectors::metadata(sceList[[i]]),
                                  reducedDims = redMatList[[i]])
  }

  if (isTRUE(combined)) {
    sce <- do.call(SingleCellExperiment::cbind, New_SCE)
    meta <- .mergeMetaSCE(New_SCE)
    S4Vectors::metadata(sce) <- meta
    return(sce)
  }
  return(New_SCE)
}
