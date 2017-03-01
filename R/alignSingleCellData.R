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
#' @param sample_annotations A data.frame of sample annotations, with samples
#' as rows and annotations in columns. The sample names must be identical to
#' and in the same order as the list of files in inputfile1. Alignment
#' statistics will be added to the annotation data frame.
#' @param probe_annotations An optional data.frame of probe annotations, with
#' probes as rows and probe annoations in columns.
#' @param threads Number of threads to use during alignment. The default is 1.
#' @param tmp_dir Specify an optional directory to save temporary files.
#' @param save_bam If TRUE, bam alignment files will be saved in the output_dir.
#' The default is FALSE.
#' @param save_count_files If TRUE, per sample gene count files will be saved in
#' the output_dir. The default is FALSE.  
#' @param ouput_dir If save_bam or save_count_files is TRUE, specify a directory
#' in which to save the output files.
#'
#' @return Object to import into the shiny app.
#' @export
#'
#' @examples
alignSingleCellData <- function(inputfile1, inputfile2, index_path,
                                gtf_annotation, ouput_dir=NULL,
                                sample_annotations=NULL, probe_annotations=NULL,
                                threads=1, tmp_dir=NULL, save_bam=FALSE,
                                save_count_files=FALSE, multiQC_data=NULL) {
  if (!requireNamespace("Rsubread", quietly = TRUE)) {
    stop("Rsubread package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  cat("done")
}