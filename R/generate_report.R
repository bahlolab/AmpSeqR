#' Generate the HTML report including various summary data visualizations
#'
#' Generate the HTML report including various summary data visualizations.
#'
#' @param sample_manifest (Required). The sample barcodes file. The file should include sample_id, barcode_fwd, barcode_rev, sample (sample name, can be the same as sample_id), info (e.g., sample type).
#' @param marker_info (Required). The target amplicon file. The file should include marker_id, primer_fwd, primer_rev, seq (reference sequence), chrom (chromosome), start (reference sequence start position), end (reference sequence end position).
#' @param demultiplexed (Required). The demultiplexed read table includes: sample_id, marker_id, reads_1 (the forward read fastq file path), reads_2 (the reverse read fastq file path), n (number of demultiplexed reads), sample, info.
#' @param flt_reads (Required). The filtered read table includes: sample_id, marker_id, n (number of demultiplexed reads), sample, info, reads_1 (the filtered and trimmed forward fastq file path), reads_2 (the filtered and trimmed reverse fastq file path), n_in (number of reads before filtering and trimming), n_out (number of reads after filtering and trimming).
#' @param sub_reads (Required). The downsampled read table includes: sample_id, marker_id, sample, info, reads_1 (the downsampled forward fastq file path), reads_2 (the downsampled reverse fastq file path), n (number of reads after downsampling).
#' @param seq_ann_tbl (Required). The amplicon haplotype sequence table includes: sample_id, marker_id, sequence (haplotype sequence), count (read counts), status (haplotype status), sample, info, ident (sequence similarity), ident_z (standardized sequence similarity).
#' @param seq_flt_tbl (Required). The final amplicon haplotype sequence table includes: sample_id, marker_id, sequence (haplotype sequence), count (read counts), masked (TRUE indicates the sequence contains variants that might be errors, and replaces the nucleotide with the reference genome nucleotide), status (haplotype status, all pass), ident (sequence similarity), haplotype (a named panel of haplotypes for each marker), frequency (within-sample haplotype frequency for each marker), sample, info.
#' @param report_output_dir (Required). The path to the output report file.
#' @param browse_report (Optional). Default TRUE. Logical indicating whether the htmlFile should be opened in a web browser.
#'
#' @return
#'
#' The HTML report including various summary data visualizations.
#'
#' @export
generate_report <- function(sample_manifest,
                            marker_info,
                            demultiplexed,
                            flt_reads,
                            sub_reads,
                            seq_ann_tbl,
                            seq_flt_tbl,
                            report_output_dir,
                            browse_report = TRUE) {
  message("---- generating report ----")

  report_dir <- report_output_dir

  if (!dir.exists(report_dir)) {
    dir.create(report_dir)
  }
  run <- list(sample_manifest, marker_info, demultiplexed, flt_reads, sub_reads, seq_ann_tbl, seq_flt_tbl, report_dir)
  names(run) <- c("sample_manifest", "marker_info", "demultiplexed", "flt_reads", "sub_reads", "seq_ann_tbl", "seq_flt_tbl", "report_dir")
  rmd_file <- system.file(file.path("rmd", "report.Rmd"), package = "AmpSeqR", mustWork = TRUE)
  rmd_env <- list2env(run, envir = new.env())
  suppressMessages(suppressWarnings(rmarkdown::render(
    input = rmd_file,
    output_dir = run$report_dir,
    output_file = str_c("Report.html"),
    envir = rmd_env
  )))

  output_file <- file.path(report_dir, str_c("Report.html"))
  if (browse_report) {
    browseURL(path.expand(output_file))
  }

  message("report save to \"", output_file, "\"")

  return(invisible(output_file))
}
