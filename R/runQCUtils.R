.checkCell <- function(FilterFile, FilterDir, basepath, Reference, process) {
    if (is.null(FilterFile)) {
        if (is.null(basepath)) {
            if ((is.null(FilterDir))) {
                stop("'cellPath' need to be specified when 'basePath' is NULL.")
            }
            # message("'base_path' is NULL. Data is loaded using directories specified by '--cell_data_path' and '--raw_data_path'.")
            if (length(FilterDir) != length(sample)) {
                stop("The length of '--cellPath' should be the same as the length of '--sample'.")
            }
            if (length(FilterDir) != length(process)) {
                stop('The length of "--cellPath" should be the same as ',
                         'the length of "--preproc"!')
            }

        } else {
            if (length(basepath) != length(process)) {
                stop('The length of "--basePath" should be the same as ',
                         'the length of "--preproc"!')
            }
            if (length(basepath) != length(sample)) {
                stop('The length of "--basePath" should be the same as ',
                         'the length of "--sample"!')
            }
            if (length(Reference) != sum(process == 'CellRangerV2')) {
                stop('The length of --genome should be the same as ',
                        'the number of "CellRangerV2" in the "--preproc"!')
            }
        }

    }

    if (!is.null(FilterFile)) {
        if (length(FilterFile) != length(sample)) {
            stop("The length of '--cellData' should be the same as the length of '--sample'.")
        }
        if (length(FilterFile) != length(process)) {
            stop('The length of "--cellData" should be the same as ',
                     'the length of "--preproc"!')
        }
    }
    return(0)
}

.checkDroplet <- function(RawFile, RawDir, basepath, Reference, process) {
    if (is.null(RawFile)) {
        if (is.null(basepath)) {
            if ((is.null(RawDir))) {
                stop("'rawPath' need to be specified when 'basePath' is NULL.")
            }
            # message("'base_path' is NULL. Data is loaded using directories specified by '--cell_data_path' and '--raw_data_path'.")
            if (length(RawDir) != length(sample)) {
                stop("The length of '--rawPath' should be the same as the length of '--sample'.")
            }
            if (length(RawDir) != length(process)) {
                stop('The length of "--rawPath" should be the same as ',
                         'the length of "--preproc"!')
            }

        } else {
            if (length(basepath) != length(process)) {
                stop('The length of "--basePath" should be the same as ',
                         'the length of "--preproc"!')
            }
            if (length(basepath) != length(sample)) {
                stop('The length of "--basePath" should be the same as ',
                         'the length of "--sample"!')
            } 
            if (length(Reference) != sum(process == 'CellRangerV2')) {
                stop('The length of --genome should be the same as ',
                        'the number of "CellRangerV2" in the "--preproc"!')
            }
        }
    }

    if (!is.null(RawFile)) {
        if (length(RawFile) != length(sample)) {
            stop("The length of '--rawData' should be the same as the length of '--sample'.")
        }
        if (length(RawFile) != length(process)) {
            stop('The length of "--rawData" should be the same as ',
                     'the length of "--preproc"!')
        }
    }
    return(0)
}

.checkBoth <- function(RawFile, FilterFile, RawDir, FilterDir, basepath, Reference, process) {
    if (is.null(RawFile) && is.null(FilterFile)) {
        if (is.null(basepath)) {
            if ((is.null(FilterDir) || is.null(RawDir))) {
                stop("Both 'cellPath' and 'rawPath' need to be specified when 'basePath' is NULL.")
            } else {
                # message("'basePath' is NULL. Data is loaded using directories specified by '--cell_data_path' and '--raw_data_path'.")
                if (length(FilterDir) != length(RawDir)) {
                    stop("The length of '--cellPath' should be the same as the length of '--rawPath'.")
                }
                if (length(FilterDir) != length(sample)) {
                    stop("The length of '--cellPath' should be the same as the length of '--sample'.")
                }
                if (length(FilterDir) != length(process)) {
                    stop('The length of "--cellPath" should be the same as ',
                             'the length of "--preproc"!')
                }
            }
        } else {
            if (length(basepath) != length(process)) {
                stop('The length of "--basePath" should be the same as ',
                         'the length of "--preproc"!')
            }
            if (length(basepath) != length(sample)) {
                stop('The length of "--basePath" should be the same as ',
                         'the length of "--sample"!')

                if (length(Reference) != sum(process == 'CellRangerV2')) {
                    stop('The length of --genome should be the same as ',
                        'the number of "CellRangerV2" in the "--preproc"!')
                }
            }
        }
    }

    if (!is.null(RawFile) || !is.null(FilterFile)) {
        if (length(RawFile) != length(FilterFile)) {
             stop("The length of '--rawData' and '--cellData' should be the same when '--preproc' is SceRDS or CountMatrix.")
        }
        if (length(FilterFile) != length(sample)) {
            stop("The length of '--cellData' should be the same as the length of '--sample'.")
        }
        if (length(FilterFile) != length(process)) {
            stop('The length of "--cellData" should be the same as ',
                     'the length of "--preproc"!')
        }
    }
    return(0)
}

.cellQC <- function(cellSCE, geneSetCollection, Params, cellQCAlgos, mitoInfo, dropletSCE = NULL) {
    p <- paste0(date(), " .. Running cell QC")
    message(p)
    cellSCE <- runCellQC(inSCE = cellSCE, geneSetCollection = geneSetCollection, 
        paramsList=Params, algorithms = cellQCAlgos, background = dropletSCE,
        mitoRef = mitoInfo[['reference']], mitoIDType = mitoInfo[['id']],
        mitoGeneLocation = "rownames")
    return(cellSCE)
}

.runCell <- function(cellSCE, preproc, geneSetCollection, Params, cellQCAlgos, mitoInfo) {
    if (is.null(cellSCE) && (preproc %in% c("BUStools", "SEQC"))) {
        dropletSCE <- runDropletQC(inSCE = dropletSCE, paramsList=Params)
        ix <- !is.na(dropletSCE$dropletUtils_emptyDrops_fdr) & (dropletSCE$dropletUtils_emptyDrops_fdr < 0.01)
        cellSCE <- dropletSCE[,ix]
    }
    p <- paste0(date(), " .. Running cell QC")
    message(p)
    cellSCE <- runCellQC(inSCE = cellSCE, geneSetCollection = geneSetCollection, 
        paramsList=Params, algorithms = cellQCAlgos, background = dropletSCE,
        mitoRef = mitoInfo[['reference']], mitoIDType = mitoInfo[['id']],
        mitoGeneLocation = "rownames")
}

.runDroplet <- function(dropletSCE, geneSetCollection, Params, mitoInfo, cellQCAlgos, detectCell, cellCalling) {
    p <- paste0(date(), " .. Running droplet QC")
    message(p)
    dropletSCE <- runDropletQC(inSCE = dropletSCE, paramsList=Params)
    if (isTRUE(detectCell)) {
        if (cellCalling == "EmptyDrops") {
            ix <- !is.na(dropletSCE$dropletUtils_emptyDrops_fdr) & dropletSCE$dropletUtils_emptyDrops_fdr < 0.01
        } else if (cellCalling == "Knee") {
            ix <- dropletSCE$dropletUtils_BarcodeRank_Knee == 1
        } else {
            ix <- dropletSCE$dropletUtils_BarcodeRank_Inflection == 1
        }
        cellSCE <- dropletSCE[,ix]
        p <- paste0(date(), " .. Running cell QC")
        message(p)
        cellSCE <- runCellQC(inSCE = cellSCE, geneSetCollection = geneSetCollection, 
            paramsList=Params, algorithms = cellQCAlgos, background = dropletSCE,
            mitoRef = mitoInfo[['reference']], mitoIDType = mitoInfo[['id']],
            mitoGeneLocation = "rownames")
    }
}


.exportHTANFlat <- function(dataType, droplet = NULL, cell = NULL, samplename, directory, formats, detectCell = FALSE, level3Meta, level4Meta, i = NULL) {
    if ((dataType == "Droplet" && isTRUE(detectCell))) {
        inputType <- "Both"
    } else {
        inputType <- dataType
    }
    if ("FlatFile" %in% formats) {
        if ("HTAN" %in% formats) {
            meta <- generateHTANMeta(dropletSCE = droplet, cellSCE = cell, samplename = samplename,
                                     dir = directory, htan_biospecimen_id=samplename, inputType)
        } else {
                meta <- generateMeta(dropletSCE = droplet, cellSCE = cell, samplename = samplename,
                                     dir = directory, HTAN=FALSE, inputType)
        }
        # in place modification of level3/4Meta
        if (!is.null(i)) {
            level3Meta[[i]] <- meta[[1]]
            level4Meta[[i]] <- meta[[2]]
        }
        else {
            level3Meta <- meta[[1]]
            level4Meta <- meta[[2]]
        }
    } else {
        warning("'FlatFile' is not in output format. Skip exporting the manifest file.")
    }
    return(meta)
}