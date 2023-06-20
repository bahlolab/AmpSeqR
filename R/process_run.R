#' Fast run
#'
#' High-level wrapper for all functions, can easily get all the results
#'
#' @param reads_1 (Required). The file path to the forward fastq file from paired-end sequence data. Compressed file formats such as .fastq.gz are supported.
#' @param reads_2 (Required). The file path to the reverse fastq file from paired-end sequence data corresponding to those provided to the reads_1 argument. Compressed file formats such as .fastq.gz are supported.
#' @param sample_manifest (Required). The sample barcodes file. The file should include sample_id, barcode_fwd, barcode_rev, sample (sample name, can be the same as sample_id), info (e.g., sample type).
#' @param marker_info (Required). The target amplicon file. The file should include marker_id, primer_fwd, primer_rev, seq (reference sequence), chrom (chromosome), start (reference sequence start position), end (reference sequence end position).
#' @param run_dir (Required). The path to the output file.
#' @param n_sample (Optional). Default 10000. Downsamples an exact number of reads from paired-end fastq files.
#' @param min_read_count (Optional). Default 1000. The minimum number of reads per sample per marker.
#' @param min_marker_count (Optional). Default 100. The minimum number of reads for each amplicon marker after calling haplotype.
#' @param min_asv_count (Optional). Default 5. The minimum number of reads for each haplotype.
#' @param min_asv_freq (Optional). Default 0.001. The minimum haplotype frequency.
#' @param min_ident (Optional). Default 0.75. The minimum sequence similarity.
#' @param min_ident_z (Optional). Default -3. The minimum standardized sequence similarity.
#' @param max_breakpoints (Optional). Default 3. The maximum number of breakpoints chimera reads.
#' @param min_parent_ratio (Optional). Default 1.5. The minimum ratio that the sequences greater than this value are more abundant than a sequence can be its "parents".
#' @param primer_max_mismatch (Optional). Default 3. The maximum number of mismatches of the amplicon primer sequence.
#' @param min_overlap (Optional). Default 10. The minimum length of the overlap required for merging the forward and reverse reads.
#' @param marker_trim (Optional). Default NULL. Data frame with three columns marker_id (character), trim_fwd (integer), and trim_rev (integer).
#' @param max_marker_miss (Optional). Default 0.5. The maximum fraction of missing data for amplicon markers.
#' @param max_sm_miss (Optional). Default 0.5. The maximum fraction of missing data for samples.
#' @param min_homo_rep (Optional). Default 3 (numeric vector). The minimum length of the homopolymer repeats.
#' @param terminal_region_len (Optional). Default 1 (numeric vector between 1 and 3). The terminal indel position.
#' @param sample_med_He (Optional). Default 0. The maximum median expected heterozygosity (a measure of genetic variation within populations) for the variance of all samples in the dataset.
#' @param n_alleles (Optional). Default 3. The maximum number of alleles at a locus.
#' @param var_maf (Optional). Default 0.001. The minimum minor allele frequency of the variant.
#' @param var_he (Optional). Default 0.001. The minimum expected heterozygosity of the variant.
#'
#' @return
#'
#' High-level wrapper for all functions, can easily get all the results
#'
#' @export
#'
#' @examples
#' example_data <- get_ampseqr_example_data()
#' reads_1 <- example_data$reads_1
#' reads_2 <- example_data$reads_2
#' sample_manifest <- example_data$sample_manifest
#' marker_info <- example_data$marker_info
#'
#' process_run_dir <- "runs_process"
#' dir.create("./runs_process")
#'
#' process_run(
#'   reads_1 = reads_1,
#'   reads_2 = reads_2,
#'   sample_manifest = sample_manifest,
#'   marker_info = marker_info,
#'   run_dir = process_run_dir
#' )
#'
process_run <- function(reads_1,
                        reads_2,
                        sample_manifest,
                        marker_info,
                        run_dir,
                        n_sample = 10000,
                        min_read_count = 1000,
                        min_marker_count = 100,
                        min_asv_count = 5,
                        min_asv_freq = 0.001,
                        min_ident = 0.75,
                        min_ident_z = -3,
                        max_breakpoints = 3L,
                        min_parent_ratio = 1.5,
                        primer_max_mismatch = 0L,
                        min_overlap = 10L,
                        marker_trim = NULL,
                        max_marker_miss = 0.5,
                        max_sm_miss = 0.5,
                        min_homo_rep = 3,
                        terminal_region_len = 1,
                        sample_med_He = 0,
                        n_alleles = 3,
                        var_maf = 0.001,
                        var_he = 0.001) {
  message("---- processing run ----")


  suppressWarnings(if (!dir.exists(run_dir)) {
    dir.create(run_dir)
  })

  message("---- demultiplexing----")

  demultiplexed <- demultiplex_reads(
    sample_manifest = sample_manifest,
    marker_info = marker_info,
    reads_1 = reads_1,
    reads_2 = reads_2,
    output_dir = run_dir,
    output_sub_dir = file.path(run_dir, "demultiplex")
  )

  message("---- filter and trim ----")

  flt_reads <- demultiplexed %>%
    dada_filter(
      output_dir = run_dir,
      output_sub_dir = file.path(run_dir, "filter")
    )


  message("---- downsampling ----")

  sub_reads <- flt_reads %>%
    downsample_reads(
      output_dir = run_dir,
      output_sub_dir = file.path(run_dir, "downsample"),
      n_sample = n_sample,
      min_read_count = min_read_count,
      count_col = "n_out",
      threads = 8
    )


  message("---- ASV estimation ----")

  seq_tbl <- sub_reads %>%
    dada_seq_tbl(output_dir = run_dir)

  # annotate sequence variants
  seq_ann_tbl <- seq_tbl %>%
    annotate_seq_tbl(
      marker_info = marker_info,
      output_dir = run_dir,
      min_marker_count = min_marker_count,
      min_asv_count = min_asv_count,
      min_asv_freq = min_asv_freq,
      min_ident = min_ident,
      min_ident_z = min_ident_z,
      max_breakpoints = max_breakpoints,
      min_parent_ratio = min_parent_ratio
    )

  message("---- filtering haplotype sequence ----")

  seq_flt_tbl <- sequence_filter(
    seq_ann_tbl = seq_ann_tbl,
    sample_manifest = sample_manifest,
    marker_info = marker_info,
    output_dir = run_dir,
    vcf_output_dir = file.path(run_dir, "vcf"),
    max_sm_miss = max_sm_miss,
    max_marker_miss = max_marker_miss,
    min_homo_rep = min_homo_rep,
    terminal_region_len = terminal_region_len,
    sample_med_He = sample_med_He,
    n_alleles = n_alleles,
    var_maf = var_maf,
    var_he = var_he
  )


  message("---- generating report ----")

  generate_report(
    sample_manifest = sample_manifest,
    marker_info = marker_info,
    demultiplexed = demultiplexed,
    flt_reads = flt_reads,
    sub_reads = sub_reads,
    seq_ann_tbl = seq_ann_tbl,
    seq_flt_tbl = seq_flt_tbl,
    report_output_dir = file.path(run_dir, "report")
  )
}
