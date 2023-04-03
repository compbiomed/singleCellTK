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
}