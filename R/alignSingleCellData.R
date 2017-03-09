#' Align Single Cell RNA-Seq Data and Create a SCESet Object
#'
#' @param inputfile1 An input file or list of files. Files can be fastq,
#' fastq.gz, or bam, but must all be of the same type. Sample names will be the
#' full file name, without _1.fastq.gz, .fastq.gz, _1.fastq, .fastq or .bam
#' endings.
#' @param inputfile2 If fastq files are provided in input list, a list of
#' corresponding paired fastq files, if applicable. 
#' @param index_path Path to the Rsubread genome index.
#' @param gtf_annotation Path to the GTF gene annotation to use. This must
#' correspond to the genome specified in index_path.
#' @param output_dir If save_bam or save_count_files is TRUE, specify a directory
#' in which to save the output files.
#' @param sample_annotations A data.frame of sample annotations, with samples
#' as rows and annotations in columns. The sample names must be identical to
#' and in the same order as the list of files in inputfile1. Alignment
#' statistics will be added to the annotation data frame.
#' @param feature_annotations An optional data.frame of probe annotations, with
#' probes as rows and probe annoations in columns.
#' @param threads Number of threads to use during alignment. The default is 1.
#' @param save_bam If TRUE, bam alignment files will be saved in the output_dir.
#' The default is FALSE.
#' @param save_count_files If TRUE, per sample gene count files will be saved in
#' the output_dir. The default is FALSE.
#' @param isPairedEnd If input files are .bam, indicate whether the input bam
#' files are paired end.
#'
#' @return Object to import into the shiny app.
#' @export alignSingleCellData
#'
#' @examples
alignSingleCellData <- function(inputfile1, inputfile2=NULL, index_path,
                                gtf_annotation, output_dir=NULL,
                                sample_annotations=NULL,
                                feature_annotations=NULL, threads=1,
                                save_bam=FALSE, save_count_files=FALSE,
                                isPairedEnd=FALSE) {
  if (!requireNamespace("Rsubread", quietly = TRUE)) {
    stop("Rsubread package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  #verify that inputfile1 file(s) exist
  if(!all(file.exists(inputfile1))){
    stop("Error with inputfile1 files, verify that all file exist.")
  }

  #if inputfile2 file(s) exist, verify them
  if(!is.null(inputfile2)){
    if(!all(file.exists(inputfile2))){
      stop("Error with inputfile2 files, verify that all file exist.")
    }
    if(length(inputfile1) != length(inputfile2)){
      stop("Error. inputfile1 and inputfile2 are different lengths.")
    }
  }
  
  #make sure the index exists
  if(!file.exists(paste(index_path, "00.b.array", sep='.'))){
    stop("Verify the Rsubread index is correctly specified.")
  }
  
  #make sure the gtf annotation exists
  if(!file.exists(gtf_annotation)){
    stop("Error with gtf annotation file, make sure the file exists.")
  }
  
  #make sure the output directory exists if save_bam or save_count_files is set
  if(save_bam | save_count_files){
    if(is.null(output_dir)){
      stop("Error with output directory. Create it to continue.")
    } else if(!dir.exists(output_dir)){
      stop("Error with output directory. Create it to continue.")
    }
  } else{
    output_dir <- tempdir()
  }
  
  countframe <- NULL

  #align all of the fastq files temp or outfile, get alignment metrics
  #make feature count files temp or outfile
  for(i in 1:length(inputfile1)){
    if(is.null(inputfile2)){
      readfile2=NULL
    } else{
      readfile2=inputfile2[i]
    }
    bam_file_path <- paste(output_dir,
                           paste(gsub("_1\\.fastq\\.gz$|\\.fastq\\.gz$|_1\\.fastq$|\\.fastq$|\\.bam$",
                                      "",
                                      basename(inputfile1[i])),
                                 "bam", sep="."), sep="/")
    if(grepl("_1\\.fastq\\.gz$|\\.fastq\\.gz$|_1\\.fastq$|\\.fastq$", inputfile1[i])){
      Rsubread::align(index=index_path,
                      readfile1=inputfile1[i],
                      readfile2=readfile2,
                      output_file=bam_file_path,
                      nthreads=threads,
                      input_format="gzFASTQ",
                      unique=TRUE,
                      indels=5)
      fCountsList <- Rsubread::featureCounts(bam_file_path,
                                             annot.ext=gtf_annotation,
                                             isGTFAnnotationFile=TRUE,
                                             nthreads=threads,
                                             isPairedEnd=!is.null(inputfile2))
    } else if(grepl("\\.bam$", inputfile1[i])){
      fCountsList <- Rsubread::featureCounts(inputfile1[i],
                                             annot.ext=gtf_annotation,
                                             isGTFAnnotationFile=TRUE,
                                             nthreads=threads,
                                             isPairedEnd=isPairedEnd)
    } else{
      stop("Input file type error. Make all files are of the supported type.")
    }

    if(is.null(countframe)){
      countframe <- data.frame(fCountsList$counts, row.names = fCountsList$annotation[,1])
      colnames(countframe)[i] <- gsub("_1\\.fastq\\.gz|\\.fastq\\.gz|_1\\.fastq|\\.fastq|\\.bam",
                                      "",
                                      basename(inputfile1[i]))
    } else {
      countframe <- cbind(countframe, data.frame(fCountsList$counts))
      colnames(countframe)[i] <- gsub("_1\\.fastq\\.gz|\\.fastq\\.gz|_1\\.fastq|\\.fastq|\\.bam",
                                      "",
                                      basename(inputfile1[i]))
    }
    if(!save_bam){
      unlink(bam_file_path)
      unlink(paste(bam_file_path, "indel", sep="."))
    }
    if(save_count_files){
      savecounts <- cbind(fCountsList$annotation[,1], fCountsList$counts)
      write.table(savecounts,
                  paste(output_dir,
                        paste(gsub("_1\\.fastq\\.gz|\\.fastq\\.gz|_1\\.fastq|\\.fastq|\\.bam",
                                   "",
                                   basename(inputfile1[i])),
                              "featureCounts",
                              sep="."),
                        sep="/"),
                  sep="\t",
                  col.names=FALSE,
                  row.names=FALSE,
                  quote=FALSE)
    }
  }
  
  #remove any gene names with empty gene name
  if(any(rownames(countframe) == "")){
    warning("One of the feature names is empty. This can be caused by a problem with your GTF file. The empty feature name will be removed.")
    countframe <- countframe[rownames(countframe) != "",,drop=F]
  }
  
  if(!is.null(sample_annotations)){
    if(!(all(rownames(sample_annotations) == colnames(countframe)))){
      warning("Sample annotation sample names do not match the countframe names. Sample annotations will not be added.")
      sample_annotations <- NULL
    }
  }
  
  if(!is.null(feature_annotations)){
    if(!(all(rownames(feature_annotations) == rownames(countframe)))){
      warning("Feature annotation names do not match the countframe features. Feature annotations will not be added.")
      feature_annotations <- NULL
    }
  }
  
  #createsceset from the count file, multiqcdata, and annotations if they exist (validate the sample names are right)
  scobject <- createSCESet(countfile = countframe,
                           annotfile = sample_annotations,
                           featurefile = feature_annotations,
                           inputdataframes = TRUE)

  return(scobject)
}
