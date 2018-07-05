#' Align Single Cell RNA-Seq Data and Create a SCtkExperiment Object
#'
#' @param inputfile1 An input file or list of files. Files can be fastq,
#' fastq.gz, or bam, but must all be of the same type. Sample names will be the
#' full file name, without _1.fastq.gz, .fastq.gz, _1.fastq, .fastq or .bam
#' endings.
#' @param inputfile2 If fastq files are provided in input list, a list of
#' corresponding paired fastq files, if applicable.
#' @param indexPath Path to the Rsubread genome index.
#' @param gtfAnnotation Path to the GTF gene annotation to use. This must
#' correspond to the genome specified in indexPath.
#' @param outputDir If saveBam or saveCountFiles is TRUE, specify a
#' directory in which to save the output files.
#' @param sampleAnnotations A data.frame of sample annotations, with samples
#' as rows and annotations in columns. The sample names must be identical to
#' and in the same order as the list of files in inputfile1. Alignment
#' statistics will be added to the annotation data frame.
#' @param featureAnnotations An optional data.frame of probe annotations, with
#' probes as rows and probe annotations in columns.
#' @param threads Number of threads to use during alignment. The default is 1.
#' @param saveBam If TRUE, bam alignment files will be saved in the outputDir.
#' The default is FALSE.
#' @param saveCountFiles If TRUE, per sample gene count files will be saved in
#' the outputDir. The default is FALSE.
#' @param isPairedEnd If input files are .bam, indicate whether the input bam
#' files are paired end.
#'
#' @return Object to import into the shiny app.
#' @export
#'
#' @examples
#' \dontrun{
#' singlecellobject <- alignSingleCellData(
#'   inputfile1 = c("/path/to/sample1_1.fastq.gz",
#'                  "/path/to/sample2_1.fastq.gz"),
#'   inputfile2 = c("/path/to/sample1_2.fastq.gz",
#'                  "/path/to/sample2_2.fastq.gz"),
#'   indexPath = "/path/to/genome/index",
#'   gtfAnnotation = "/path/to/gene/annotations.gtf",
#'   sampleAnnotations = sample.annotation.df,
#'   threads=4)}
#'
alignSingleCellData <- function(inputfile1, inputfile2=NULL, indexPath,
                                gtfAnnotation, outputDir=NULL,
                                sampleAnnotations=NULL,
                                featureAnnotations=NULL, threads=1,
                                saveBam=FALSE, saveCountFiles=FALSE,
                                isPairedEnd=FALSE) {
  if (!requireNamespace("Rsubread", quietly = TRUE)) {
    stop("The Rsubread package is required to align data with this function. ",
         "Install Rsubread to continue. NOTE: Rsubread is not supported on ",
         "Windows.", call. = FALSE)
  }

  if (any(grepl("~", c(inputfile1, inputfile2, indexPath, gtfAnnotation,
                       outputDir)), na.rm = TRUE)){
    stop("ERROR: tilde ~ filenames are not supported. Provide the full path.")
  }

  #verify that inputfile1 file(s) exist
  if (!all(file.exists(inputfile1))){
    stop("Error with inputfile1 files, verify that all file exist.")
  }

  #if inputfile2 file(s) exist, verify them
  if (!is.null(inputfile2)){
    if (!all(file.exists(inputfile2))){
      stop("Error with inputfile2 files, verify that all file exist.")
    }
    if (length(inputfile1) != length(inputfile2)){
      stop("Error. inputfile1 and inputfile2 are different lengths.")
    }
  }

  #make sure the index exists
  if (!file.exists(paste(indexPath, "00.b.array", sep = "."))){
    stop("Verify the Rsubread index is correctly specified.")
  }

  #make sure the gtf annotation exists
  if (!file.exists(gtfAnnotation)){
    stop("Error with gtf annotation file, make sure the file exists.")
  }

  #make sure the output directory exists if saveBam or saveCountFiles is set
  if (saveBam | saveCountFiles){
    if (is.null(outputDir)){
      stop("Error with output directory. Create it to continue.")
    } else if (!dir.exists(outputDir)){
      stop("Error with output directory. Create it to continue.")
    }
  } else{
    outputDir <- tempdir()
  }

  countframe <- NULL
  rsubreadStats <- NULL

  #align all of the fastq files temp or outfile, get alignment metrics
  #make feature count files temp or outfile
  for (i in seq_along(inputfile1)){
    sampleName <- gsub(
      "_1\\.fastq\\.gz$|\\.fastq\\.gz$|_1\\.fastq$|\\.fastq$|\\.bam$", "",
      basename(inputfile1[i]))
    message("Processing: ", sampleName)
    rsubLogFile <- NULL
    featureLogFile <- NULL

    if (is.null(inputfile2)){
      readfile2 <- NULL
    } else{
      readfile2 <- inputfile2[i]
    }
    bamFilePath <- paste(outputDir,
                           paste(sampleName, "bam", sep = "."), sep = "/")
    if (grepl("_1\\.fastq\\.gz$|\\.fastq\\.gz$|_1\\.fastq$|\\.fastq$",
              inputfile1[i])){
      rsubLogFile <- tempfile()
      rsubLog <- file(rsubLogFile, open = "wt")
      sink(rsubLog, type = "output")
      Rsubread::align(index = indexPath,
                      readfile1 = inputfile1[i],
                      readfile2 = readfile2,
                      output_file = bamFilePath,
                      nthreads = threads,
                      input_format = "gzFASTQ",
                      unique = TRUE,
                      indels = 5)
      sink()
      close(rsubLog)

      featureLogFile <- tempfile()
      featureLog <- file(featureLogFile, open = "wt")
      sink(featureLog, type = "output")
      fCountsList <- Rsubread::featureCounts(bamFilePath,
                                             annot.ext = gtfAnnotation,
                                             isGTFAnnotationFile = TRUE,
                                             nthreads = threads,
                                             isPairedEnd = !is.null(inputfile2))
      sink()
      close(featureLog)

    } else if (grepl("\\.bam$", inputfile1[i])){
      featureLogFile <- tempfile()
      featureLog <- file(featureLogFile, open = "wt")
      sink(featureLog, type = "output")
      fCountsList <- Rsubread::featureCounts(inputfile1[i],
                                             annot.ext = gtfAnnotation,
                                             isGTFAnnotationFile = TRUE,
                                             nthreads = threads,
                                             isPairedEnd = isPairedEnd)
      sink()
      close(featureLog)
    } else{
      stop("Input file type error. Make all files are of the supported type.")
    }

    rsubreadStatsLine <- parseRsubreadLogs(rsubLogFile, featureLogFile,
                                           sampleName)
    unlink(rsubLogFile)
    unlink(featureLogFile)

    #add the feature counts to the countframe
    if (is.null(countframe)){
      countframe <- data.frame(fCountsList$counts,
                               row.names = fCountsList$annotation[, 1])
      colnames(countframe)[i] <- sampleName
    } else {
      countframe <- cbind(countframe, data.frame(fCountsList$counts))
      colnames(countframe)[i] <- sampleName
    }

    #add readsMapped to rsubreadStats
    if (is.null(rsubreadStats)){
      rsubreadStats <- rsubreadStatsLine
    } else {
      rsubreadStats <- rbind(rsubreadStats, rsubreadStatsLine)
    }

    if (!saveBam){
      unlink(bamFilePath)
      unlink(paste(bamFilePath, "indel", sep = "."))
    }
    if (saveCountFiles){
      savecounts <- cbind(fCountsList$annotation[, 1], fCountsList$counts)
      utils::write.table(savecounts, paste(outputDir,
                                           paste(sampleName, "featureCounts",
                                                 sep = "."), sep = "/"),
                         sep = "\t", col.names = FALSE, row.names = FALSE,
                         quote = FALSE)
    }
  }

  #remove any gene names with empty gene name
  if (any(rownames(countframe) == "")){
    warning("One of the feature names is empty. This can be caused by a ",
            "problem with your GTF file. The empty feature name will be ",
            "removed.")
    countframe <- countframe[rownames(countframe) != "", , drop = FALSE]
  }

  if (!is.null(sampleAnnotations)){
    if (!(all(rownames(sampleAnnotations) == colnames(countframe)))){
      warning("Sample annotation sample names do not match the countframe ",
              "names. Sample annotations will not be added.")
      sampleAnnotations <- rsubreadStats
    } else{
      sampleAnnotations <- cbind(sampleAnnotations, rsubreadStats)
    }
  } else {
    sampleAnnotations <- rsubreadStats
  }

  if (!is.null(featureAnnotations)){
    if (!(all(rownames(featureAnnotations) == rownames(countframe)))){
      warning("Feature annotation names do not match the countframe features. ",
              "Feature annotations will not be added.")
      featureAnnotations <- NULL
    }
  }

  #createsceset from the count file, multiqcdata, and annotations if they exist
  # (validate the sample names are right)
  scobject <- createSCE(assayFile = countframe,
                        annotFile = sampleAnnotations,
                        featureFile = featureAnnotations,
                        inputDataFrames = TRUE)

  return(scobject)
}


#' Parse Rsubread Logs for Mapping and Feature Count Statistics
#'
#' @param alignLog Path to a log file created by the Rsubread align function
#' @param featurecountLog Path to a log file created by the Rsubread feature
#' count function
#' @param sampleName Sample name corresponding to the two log files
#'
#' @return A single line of a data frame with alignment and feature count
#' information
#' @export
#'
parseRsubreadLogs <- function(alignLog=NULL, featurecountLog=NULL,
                              sampleName=NULL){
  #process feature count log num reads and log num featured
  fFh <- file(featurecountLog, open = "r")
  features <- readLines(fFh)
  close(fFh)
  totalLine <- unlist(strsplit(features[grep("Total reads ", features)],
                                " +", perl = TRUE))
  featureLine <- unlist(strsplit(features[grep("Successfully assigned reads",
                                                features)], " +", perl = TRUE))
  totalReads <- as.numeric(gsub(",", "", totalLine[grep("reads",
                                                        totalLine) + 2]))
  featureReads <- as.numeric(gsub(",", "",
                                   featureLine[grep("reads",
                                                    featureLine) + 2]))
  #if not null align log
  if (!is.null(alignLog)){
    #process align log
    aFh <- file(alignLog, open = "r")
    align <- readLines(aFh)
    close(aFh)
    mapLine <- unlist(strsplit(align[grep("Mapped", align)], " +",
                                perl = TRUE))
    mappedReads <- as.numeric(gsub(",", "", mapLine[grep("reads",
                                                           mapLine) - 1]))
    return(data.frame(totalReads = totalReads,
                      readsMapped = mappedReads,
                      pctMapped = mappedReads / totalReads,
                      readsAssignedToFeatures = featureReads,
                      pctFeatures = featureReads / totalReads,
                      row.names = sampleName))
  } else{
    #return just feature stats
    return(data.frame(totalReads = totalReads,
                      readsAssignedToFeatures = featureReads,
                      pctFeatures = featureReads / totalReads,
                      row.names = sampleName))
  }
}
