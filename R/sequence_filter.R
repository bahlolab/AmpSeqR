#' Optimize amplicon haplotype calling
#'
#' This function optimize amplicon haplotype calling and generate an amplicon haplotype sequence table.
#'
#' @param seq_ann_tbl (Required). The amplicon haplotype sequence table includes: sample_id, marker_id, sequence (haplotype sequence), count (read counts), status (haplotype status), sample, info, ident (sequence similarity), ident_z (standardized sequence similarity).
#' @param sample_manifest (Required). The sample barcodes file. The file should include sample_id, barcode_fwd, barcode_rev, sample (sample name, can be the same as sample_id), info (e.g., sample type).
#' @param marker_info (Required). The target amplicon file. The file should include marker_id, primer_fwd, primer_rev, seq (reference sequence), chrom (chromosome), start (reference sequence start position), end (reference sequence end position).
#' @param output_dir (Required). The path to the output sequence table file.
#' @param vcf_output_dir (Required). The path to the output vcf files.
#' @param min_homo_rep (Optional). Default 3. The minimum size for homopolymer in reference sequence.
#' @param terminal_region_len (Optional). Default 1. The length of the terminal region. The maximum value is 3.
#' @param max_sm_miss (Optional). Default 0.5. The maximum fraction of missing data for samples.
#' @param max_marker_miss (Optional). Default 0.5. The maximum fraction of missing data for amplicon markers.
#' @param sample_med_He (Optional). Default 0. The maximum median expected heterozygosity (a measure of genetic variation within populations) for the variance of all samples in the dataset.
#' @param n_alleles (Optional). Default 3. The maximum number of alleles at a locus.
#' @param var_maf (Optional). Default 0.001. The minimum minor allele frequency (MAF) of the variant.
#' @param var_he (Optional). Default 0.001. The minimum expected heterozygosity of the variant.
#'
#' @return
#'
#' An amplicon haplotype sequence table in RDS format which includes sample_id, marker_id, sequence (haplotype sequence), count (read counts), masked (TRUE indicates the sequence contains variants that might be errors, and replaces the nucleotide with the reference genome nucleotide), status (haplotype status, all pass), ident (sequence similarity), haplotype (a named panel of haplotypes for each marker), frequency (within-sample haplotype frequency for each marker), sample, info.
#'
#' @export
#' @importFrom rlang is_scalar_integerish is_scalar_character
#' @importFrom dplyr mutate filter left_join select if_else case_when summarise group_by arrange_all anti_join semi_join row_number
#' @importFrom magrittr "%>%"
#' @importFrom Biostrings DNAStringSet
#' @importFrom DECIPHER AlignSeqs
#' @importFrom purrr map2 pmap_df
#' @importFrom stringr str_c str_remove_all
#' @importFrom tidyr pivot_longer replace_na complete
#' @importFrom withr local_options
sequence_filter <- function(seq_ann_tbl,
                            sample_manifest,
                            marker_info,
                            output_dir,
                            vcf_output_dir,
                            min_homo_rep = 3,
                            terminal_region_len = 1,
                            max_sm_miss = 0.5,
                            max_marker_miss = 0.5,
                            sample_med_He = 0,
                            n_alleles = 3,
                            var_maf = 0.001,
                            var_he = 0.001) {
  # Warining message
  withr::local_options(.new = list(warn = -1))
  # check args
  stopifnot(
    is.data.frame(seq_ann_tbl),
    is.data.frame(sample_manifest),
    is.data.frame(marker_info)
  )
  
  check_seq_table(seq_ann_tbl)
  check_sample_manifest(sample_manifest)
  check_marker_info(marker_info)
  
  # Calculate the missingness
  missigness <-
    seq_ann_tbl %>%
    select(sample_id, marker_id) %>%
    distinct() %>%
    mutate(missing = FALSE) %>%
    group_by(sample_id) %>%
    complete(marker_id = marker_info$marker_id, fill = list(missing = TRUE)) %>%
    group_by(marker_id) %>%
    complete(sample_id = sample_manifest$sample_id, fill = list(missing = TRUE)) %>%
    ungroup()
  
  percent_missing <-
    missigness$missing %>%
    (function(x) 100 * sum(x) / length(x))
  
  sample_pass <-
    missigness %>%
    group_by(sample_id) %>%
    summarise(missing_rate = sum(missing) / n()) %>%
    filter(missing_rate < max_sm_miss) %>%
    select(sample_id)
  
  marker_pass <-
    suppressMessages(missigness %>%
                       semi_join(sample_pass) %>%
                       group_by(marker_id) %>%
                       summarise(missing_rate = sum(missing) / n()) %>%
                       filter(missing_rate < max_marker_miss))
  
  
  
  # clean_homopolymers
  # minimum size for homopolymer in reference sequence
  if (is.null(min_homo_rep)) {
    seq_ann_tbl
  } else {
    seq_ann_tbl <- clean_homopolymers(
      seq_ann_tbl = seq_ann_tbl,
      marker_info = marker_info,
      min_homo_rep = min_homo_rep
    )
  }
  
  # clean terminal indels
  if (is.null(terminal_region_len)) {
    seq_ann_tbl
  } else {
    seq_ann_tbl <- clean_terminal_indels(
      seq_ann_tbl = seq_ann_tbl,
      marker_info = marker_info,
      terminal_region_len = terminal_region_len
    )
  }
  
  seq_tbl_all <- suppressMessages(seq_ann_tbl %>%
                                    filter(status == "pass") %>%
                                    semi_join(sample_pass, "sample_id") %>%
                                    semi_join(marker_pass, "marker_id") %>%
                                    mutate(seq_id = map_chr(sequence, digest::digest)))
  
  vars <- call_vars(
    seq_tbl = seq_tbl_all,
    marker_info = marker_info
  )
  
  total_var <- nrow(vars)
  var_stats <- calc_var_stats(seq_tbl_all, vars)
  
  
  # Check if there are marker all same as reference genome
  allsame_ref_marker_1 <- setdiff(marker_pass$marker_id, unique(var_stats$varwise$marker_id))
  allsame_ref_marker_2 <- setdiff(marker_pass$marker_id, unique(var_stats$samplewise$marker_id))
  allsame_ref_marker <- union(allsame_ref_marker_1, allsame_ref_marker_2)
  
  if (length(allsame_ref_marker) > 0) {
    seq_tbl_sameRef <- seq_tbl_all %>%
      filter(marker_id %in% allsame_ref_marker)
    
    haplotype_sameRef <- seq_tbl_all %>%
      filter(marker_id %in% allsame_ref_marker) %>%
      group_by(marker_id, sequence) %>%
      summarise(count = n()) %>%
      arrange(desc(count)) %>%
      ungroup()
    
    haplotype_split_sameRef <- haplotype_sameRef %>%
      split(haplotype_sameRef$marker_id)
    
    haplotype_split_sameRef <- lapply(seq_along(haplotype_split_sameRef), function(i) {
      haplotype_split_sameRef[[i]] %>%
        mutate(Number = row_number()) %>%
        mutate(haplotype = paste0(unique(marker_id), "-", Number, sep = "")) %>%
        dplyr::select(marker_id, sequence, haplotype)
    })
    
    haplotype_split_sameRef <- do.call(rbind, haplotype_split_sameRef)
    seq_tbl_masked_sameRef <- suppressMessages(seq_tbl_sameRef %>%
                                                 left_join(haplotype_split_sameRef, by = c("marker_id", "sequence")))
    
    seq_tbl_masked_sameRef <- suppressMessages(seq_tbl_masked_sameRef %>%
                                                 group_by(sample_id, marker_id) %>%
                                                 mutate(frequency = round(count / sum(count), 5)) %>%
                                                 ungroup() %>%
                                                 mutate(masked = FALSE) %>%
                                                 select(sample_id, marker_id, sequence, count, masked, status, ident, haplotype, frequency, sample, info))
  } else {
    seq_tbl_masked_sameRef <- NULL
  }
  
  
  
  if (var_maf == 0 | var_he == 0) {
    abort(str_c("The var_maf and var_he must large than 0"))
  }
  
  vars_to_mask <-
    var_stats$varwise %>%
    filter(sample_med_exp_het > sample_med_He | n_allele > n_alleles | (MAF < var_maf & MAF > 0) | (exp_het < var_he & exp_het > 0)) %>%
    select(marker_id, var_id, aln_pos, type)
  
  if (nrow(vars_to_mask) > 0 & nrow(vars_to_mask) == total_var) {
    warning("Please reset sample_med_He to a high value")
    
    seq_tbl_masked <- seq_tbl_all %>%
      mutate(masked = FALSE) %>%
      select(sample_id, marker_id, seq_id, sequence, count, masked, status, ident, ident_z)
    
    vars_masked <- call_vars(
      seq_tbl = seq_tbl_masked,
      marker_info = marker_info
    )
    var_masked_stats <- calc_var_stats(seq_tbl_masked, vars_masked)
  } else if (nrow(vars_to_mask) > 0) {
    seq_tbl_masked <-
      seq_mask_vars(
        seq_tbl = seq_tbl_all,
        marker_info = marker_info,
        vars_mask = vars_to_mask
      )
    
    vars_masked <- call_vars(
      seq_tbl = seq_tbl_masked,
      marker_info = marker_info
    )
    var_masked_stats <- calc_var_stats(seq_tbl_masked, vars_masked)
  } else {
    seq_tbl_masked <- seq_tbl_all %>%
      mutate(masked = FALSE) %>%
      select(sample_id, marker_id, seq_id, sequence, count, masked, status, ident, ident_z)
    
    vars_masked <- call_vars(
      seq_tbl = seq_tbl_masked,
      marker_info = marker_info
    )
    var_masked_stats <- calc_var_stats(seq_tbl_masked, vars_masked)
  }
  
  # export VCF, for compatibility with other tools such as SnpEff
  
  if (!dir.exists(vcf_output_dir)) {
    dir.create(vcf_output_dir)
  }
  vcf_tbl <- export_vcf(vars_masked, marker_info = marker_info, output_dir = vcf_output_dir)
  
  haplotype <- seq_tbl_masked %>%
    group_by(marker_id, sequence) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) %>%
    ungroup()
  
  haplotype_split <- haplotype %>%
    split(haplotype$marker_id)
  
  haplotype_split <- lapply(seq_along(haplotype_split), function(i) {
    haplotype_split[[i]] %>%
      mutate(Number = row_number()) %>%
      mutate(haplotype = paste0(unique(marker_id), "-", Number, sep = "")) %>%
      dplyr::select(marker_id, sequence, haplotype)
  })
  
  haplotype_split <- do.call(rbind, haplotype_split)
  seq_tbl_masked <- suppressMessages(seq_tbl_masked %>%
                                       left_join(haplotype_split, by = c("marker_id", "sequence")))
  
  seq_tbl_masked <- suppressMessages(seq_tbl_masked %>%
                                       group_by(sample_id, marker_id) %>%
                                       mutate(frequency = round(count / sum(count), 5)) %>%
                                       ungroup() %>%
                                       left_join((sample_manifest %>% select(sample_id, sample, info)), by = "sample_id") %>%
                                       select(-c(seq_id, ident_z)))
  
  
  if (length(allsame_ref_marker) > 0) {
    seq_tbl_masked <- do.call(rbind, list(seq_tbl_masked, seq_tbl_masked_sameRef))
    seq_tbl_masked <- seq_tbl_masked %>%
      arrange(sample_id, marker_id)
  } else {
    seq_tbl_masked <- seq_tbl_masked
  }
  
  write_rds(seq_tbl_masked, paste(output_dir, "/seq_flt_tbl.rds", sep = ""))
  return(seq_tbl_masked)
}
