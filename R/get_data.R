#' Get AmpSeqR extdata
#'
#' The example data that is required for AmpSeqR. The input files are the standard paired-end FASTQ files provided by the common NGS sequencing platforms, as well as sample barcodes and target amplicon details. The sample barcodes file should include sample_id, barcode_fwd, barcode_rev, sample (sample name, can be the same as sample_id), info (e.g., sample type).  The target amplicon details file should include marker_id, primer_fwd, primer_rev, seq (reference sequence), chrom (chromosome), start (reference sequence start position), end (reference sequence end position).
#'
#' @return
#'
#' List of the example data.
#'
#' @export
#'
#' @examples
#' example_data <- get_ampseqr_example_data()
get_ampseqr_example_data <- function() {
  data_dir <- system.file(package = "AmpSeqR", "extdata")
  
  reads_1 <- file.path(data_dir, "readsF.fastq.gz")
  
  reads_2 <- file.path(data_dir, "readsR.fastq.gz")
  
  sample_manifest <- readr::read_csv(file.path(data_dir, "sample_manifest.csv"), col_types = readr::cols())
  
  marker_info <- readr::read_csv(file.path(data_dir, "marker_info.csv"), col_types = readr::cols(start = readr::col_integer()))
  
  as.list(environment())
}
