#' dada_filter
#'
#' Translate dada2::filterAndTrim to ampseqr::dada_filter. Filter the poor-quality reads and trim low-quality bases.
#'
#' @param read_table (Required). The demultiplexed read table includes: sample_id, marker_id, reads_1 (the forward read fastq file path), reads_2 (the reverse read fastq file path), n (number of demultiplexed reads), sample, info.
#' @param output_dir (Required). The path to the output filtered table file.
#' @param output_sub_dir (Required). The path to the output filtered files.
#' @param max_exp_errors (Required). The same as in dada2::filterAndTrim().
#' @param min_len (Optional). Default 25 (numeric vector). Remove reads with length less than minLen.
#' @param marker_trim (Optional). Default NULL. If trimming is to be performed, marker_trim should be a data frame with three columns marker_id (character), trim_fwd (integer), and trim_rev (integer).
#' @param threads (Optional). Default 1. If an integer is provided, it is passed to the mc.cores argument of mclapply.
#'
#'
#' @return
#'
#' filter folder: filtered and trimmed paired-end fastq files.
#'
#' filtered_reads.rds table: the filtered table in RDS format which includes sample_id, marker_id, n (number of demultiplexed reads), sample, info, reads_1 (the filtered and trimmed forward fastq file path), reads_2 (the filtered and trimmed reverse fastq file path), n_in (number of reads before filtering and trimming), n_out (number of reads after filtering and trimming).
#'
#'
#' @export
#' @importFrom dplyr select mutate filter
#' @importFrom rlang is_scalar_integerish is_scalar_character abort
#' @importFrom magrittr "%>%"
#'
dada_filter <- function(read_table,
                        output_dir,
                        output_sub_dir,
                        max_exp_errors = 1L,
                        min_len = 25L,
                        marker_trim = NULL,
                        threads = 1L) {

  # check args
  stopifnot(
    is.data.frame(read_table),
    is_scalar_character(output_sub_dir),
    is_scalar_integerish(min_len) && min_len >= 0L,
    is_scalar_integerish(max_exp_errors) && max_exp_errors >= 0L,
    is_scalar_integerish(threads) && threads >= 1L
  )

  check_character_table(read_table, c("reads_1", "reads_2"), "read_table")
  # remove any missing reads
  read_table <- filter(
    read_table,
    !is.na(read_table$reads_1),
    !is.na(read_table$reads_2)
  ) %>%
    filter(!is.na(sample_id), !is.na(marker_id), n > 0)

  check_files_exist(read_table$reads_1, "read_table$reads_1")
  check_files_exist(read_table$reads_2, "read_table$reads_2")

  # create output directory
  if (!dir.exists(output_sub_dir)) {
    dir.create(output_sub_dir, recursive = T)
  }
  output_sub_dir <- normalizePath(output_sub_dir)

  ret_table <-
    read_table %>%
    mutate(
      reads_1_out = file.path(output_sub_dir, basename(reads_1)),
      reads_2_out = file.path(output_sub_dir, basename(reads_2))
    )

  # check distinct read names
  if (any(c(ret_table$reads_1, ret_table$reads_2) == c(ret_table$reads_1_out, ret_table$reads_2_out))) {
    abort(str_c("output_sub_dir must be different from input reads dir"))
  }

  # check read trimming to specified length by marker
  if (!is.null(marker_trim)) {
    filt_res <-
      suppressMessages(ret_table %>%
        inner_join(marker_trim) %>%
        pmap_df(function(reads_1, reads_1_out, reads_2, reads_2_out, trim_fwd, trim_rev, ...) {
          data.frame(dada2::filterAndTrim(
            fwd = reads_1, filt = reads_1_out, rev = reads_2, filt.rev = reads_2_out,
            truncQ = 0, minLen = min_len, maxEE = max_exp_errors, rm.phix = FALSE,
            truncLen = c(trim_fwd, trim_rev),
            multithread = threads
          ))
        }))
  } else {
    filt_res <-
      ret_table %>%
      with(dada2::filterAndTrim(
        fwd = reads_1, filt = reads_1_out, rev = reads_2, filt.rev = reads_2_out,
        truncQ = 0, minLen = min_len, maxEE = max_exp_errors, rm.phix = FALSE,
        multithread = threads
      ))
  }


  ret_table <-
    ret_table %>%
    select(-reads_1, -reads_2) %>%
    rename(reads_1 = reads_1_out) %>%
    rename(reads_2 = reads_2_out) %>%
    mutate(n_in = filt_res[, 1], n_out = filt_res[, 2]) %>%
    mutate(
      reads_1 = if_else(n_out > 0, reads_1, NA_character_),
      reads_2 = if_else(n_out > 0, reads_2, NA_character_)
    )

  write_rds(ret_table, paste(output_dir, "/filtered_reads.rds", sep = ""))

  return(ret_table)
}


#' Construct an amplicon sequence variant table (ASV) table
#'
#' Construct an amplicon sequence variant table (ASV) table.
#'
#' @param read_table The filtered read table includes: sample_id, marker_id, sample, info, reads_1 (the forward read fastq file path), reads_2 (the reverse read fastq file path), n_out (number of reads after filtered), n.
#' @param output_dir (Required). The path to the output file.
#' @param threads (Optional). The hostnames of workers (as a character vector) or the number of localhost workers (as a positive integer).
#' @param min_overlap (Optional). Default 10. The minimum length of the overlap required for merging the forward and reverse reads. When paired-end reads do not overlap, the min_overlap parameter should be set to -1.
#' @param max_mismatch (Optional). Default 0. The maximum mismatches allowed in the overlap region.
#' @param trim_overhang (Optional). Default FALSE. If TRUE, "overhangs" in the alignment between the forwards and reverse read are trimmed off. "Overhangs" are when the reverse read extends past the start of the forward read, and vice-versa, as can happen when reads are longer than the amplicon and read into the other-direction primer region.
#' @param homo_gap_value (Optional). Default NULL (no special homopolymer penalty). The alignment gap penalty within homopolymer regions. Should be negative.
#' @param dada_opt (Optional). setDadaOpt sets the default options used by the dada2::dada(...) function for your current session
#' @param marker_info (Optional). Default NULL. The target amplicon file. The file should include marker_id, primer_fwd, primer_rev, seq (reference sequence), chrom (chromosome), start (reference sequence start position), end (reference sequence end position).
#'
#' @return
#'
#' seq_tbl.rds table: the haplotype sequence table table in RDS format includes sample_id, marker_id, sequence (haplotype sequence), count (haplotype read counts), status (haplotype status), sample, info.
#'
#' @export
#'
#' @importFrom dplyr select mutate filter group_by ungroup summarise
#' @importFrom tidyr pivot_longer pivot_wider chop unnest_legacy unnest
#' @importFrom rlang is_scalar_integerish is_scalar_character is_scalar_double abort
#' @importFrom purrr map map2 pmap
dada_seq_tbl <- function(read_table,
                         output_dir,
                         threads = 1L,
                         min_overlap = 10L,
                         max_mismatch = 1L,
                         trim_overhang = TRUE,
                         homo_gap_value = NULL,
                         dada_opt = list(
                           HOMOPOLYMER_GAP_PENALTY = -1,
                           MIN_ABUNDANCE = 5
                         ),
                         marker_info = NULL) {
  # check args
  stopifnot(
    is.data.frame(read_table),
    is_scalar_integerish(threads) && threads >= 1L,
    is_scalar_integerish(min_overlap) && (min_overlap >= 5 || min_overlap == -1),
    is_scalar_integerish(max_mismatch) && max_mismatch >= 0L,
    is_bool(trim_overhang),
    is.null(dada_opt) || is.list(dada_opt)
  )

  check_character_table(read_table, c("sample_id", "marker_id", "sample", "info", "reads_1", "reads_2"), "read_table")
  # remove any missing reads
  read_table <- filter(
    read_table,
    !is.na(read_table$reads_1),
    !is.na(read_table$reads_2)
  ) %>%
    select(sample_id, marker_id, sample, info, reads_1, reads_2)

  # change
  sample_info <- read_table %>%
    select(sample_id, marker_id, sample, info) %>%
    distinct()

  check_files_exist(read_table$reads_1, "read_table$reads_1")
  check_files_exist(read_table$reads_2, "read_table$reads_2")


  cluster <- `if`(
    threads > 1,
    future::makeClusterPSOCK(workers = threads),
    NULL
  )

  if (!is.null(dada_opt)) {
    old_dada_opt <- dada2::getDadaOpt()
    do.call(dada2::setDadaOpt, dada_opt)
    if (threads > 1) {
      parallel::clusterExport(cluster, "dada_opt", envir = environment())
      invisible(parallel::clusterEvalQ(cluster, {
        do.call(dada2::setDadaOpt, dada_opt)
      }))
    }
  }

  # dereplicate reads with dada2::derepFastq
  derep_tbl <-
    read_table %>%
    pivot_longer(c(reads_1, reads_2),
      names_to = "read_set",
      names_prefix = "reads_",
      values_to = "read_file"
    ) %>%
    (function(x) {
      if (threads > 1) {
        mutate(x, group = seq_along(read_file) %% threads) %>%
          chop(-group) %>%
          mutate(derep = map(read_file, function(rf) {
            future::cluster(
              {
                dada2::derepFastq(rf)
              },
              workers = cluster
            )
          }) %>% map(future::value)) %>%
          unnest_legacy() %>%
          select(-group)
      } else {
        mutate(x, derep = dada2::derepFastq(read_file))
      }
    })

  # denoise amplicons with dada2::dada
  dada_tbl <-
    derep_tbl %>%
    chop(c(-marker_id, -read_set)) %>%
    (function(x) {
      if (threads > 1) {
        x %>%
          mutate(dada = map(derep, function(dr) {
            future::cluster(
              {
                err <- dada2::learnErrors(dr, verbose = 0)
                dada2::dada(dr, err, verbose = 0)
              },
              workers = cluster
            )
          }) %>% map(future::value))
      } else {
        x %>%
          mutate(dada = map(derep, function(dr) {
            err <- dada2::learnErrors(dr, verbose = 0)
            dada2::dada(dr, err, verbose = 0)
          }))
      }
    }) %>%
    mutate(dada = map(dada, ~ `if`(is(., "dada"), list(.), .)))

  # merge pairs with dada2::mergePairs
  # add status summaries
  seq_tbl <-
    dada_tbl %>%
    select(marker_id, read_set, sample_id, derep, dada) %>%
    unnest_legacy() %>%
    pivot_wider(names_from = read_set, values_from = c(derep, dada)) %>%
    mutate(merged = pmap(., function(derep_1, derep_2, dada_1, dada_2, ...) {

      # if (min_overlap != -1)
      if (min_overlap != -1) {
        suppressMessages(
          dada2::mergePairs(
            dadaF = dada_1, derepF = derep_1, dadaR = dada_2, derepR = derep_2, verbose = FALSE,
            minOverlap = min_overlap, trimOverhang = trim_overhang, maxMismatch = max_mismatch, homo_gap = homo_gap_value
          )
        ) %>%
          (function(x) {
            `if`(
              nrow(x) == 0,
              data.frame(sequence = NA_character_, status = "failed_dada"),
              mutate(x, status = "pass")
            )
          })
      } else {
        suppressMessages(
          dada2::mergePairs(
            dadaF = dada_1, derepF = derep_1, dadaR = dada_2, derepR = derep_2, verbose = FALSE,
            minOverlap = min_overlap, trimOverhang = trim_overhang, maxMismatch = max_mismatch, homo_gap = homo_gap_value, justConcatenate = TRUE
          )
        ) %>%
          (function(x) {
            `if`(
              nrow(x) == 0,
              data.frame(sequence = NA_character_, status = "failed_dada"),
              mutate(x, status = "pass")
            )
          })
      }
    })) %>%
    select(sample_id, marker_id, merged) %>%
    unnest(merged) %>%
    select(sample_id, marker_id, sequence, count = abundance, status) %>%
    group_by(sample_id, marker_id, sequence) %>%
    summarise(count = sum(count), status = status[1], .groups = "drop")

  # replace 10 Ns with number of bases missing based on marker info and missing reads based on the marker info
  # seq:            CATG------------NNNNNNNNNNCATG
  # marker:         CATGTACGATATATATATATATATATCATG
  # Replace seq as: CATGTACGATATATATATATATATATCATG
  seq_tbl <- if (min_overlap != -1) {
    seq_tbl
  } else {
    seq_tbl_marker_split <- suppressMessages(seq_tbl %>%
      left_join(select(marker_info, c(marker_id, seq)), by = "marker_id") %>%
      rename(ref = seq) %>%
      split.data.frame(.$marker_id))

    for (n in seq_along(seq_tbl_marker_split)) {
      markerseq <- DNAStringSet(unique(seq_tbl_marker_split[[n]]$ref))
      names(markerseq) <- "Marker"
      sequence <- DNAStringSet(as.character(seq_tbl_marker_split[[n]]$sequence))
      align_marker <- suppressMessages(DNAStringSet(c(sequence, markerseq)) %>%
        DECIPHER::AlignSeqs(verbose = FALSE, processors = 1) %>%
        as.matrix() %>%
        t() %>%
        as_tibble(.name_repair = "unique"))

      for (i in 1:(ncol(align_marker) - 1)) {
        select_miss <- which(align_marker[, i] == "-")
        select_N <- which(align_marker[, i] == "N")
        select_miss <- select_miss[select_miss > select_N[1] & select_miss < select_N[length(select_N)]]
        select_all <- c(select_miss, select_N)
        align_marker[select_all, i] <- align_marker[select_all, "Marker"]
      }

      align_marker <- align_marker %>%
        select(-Marker) %>%
        as.data.frame()

      align_sequence <- lapply(1:ncol(align_marker), function(i) {
        str_c(align_marker[, i][!(align_marker[, i] %in% "-")], collapse = "")
      })
      align_sequence <- as.data.frame(do.call(rbind, align_sequence))

      for (i in 1:length(seq_tbl_marker_split[[n]]$sequence)) {
        seq_tbl_marker_split[[n]]$sequence[i] <- align_sequence$V1[i]
      }
    }
    seq_tbl_marker_combine <- do.call(rbind, seq_tbl_marker_split)
    seq_tbl_marker_combine <- seq_tbl_marker_combine %>%
      select(-ref) %>%
      mutate_at("sequence", as.character)
    seq_tbl_marker_combine
  }

  if (threads > 1) parallel::stopCluster(cluster)

  if (!is.null(dada_opt)) do.call(dada2::setDadaOpt, old_dada_opt)

  # change
  seq_tbl <- seq_tbl %>%
    left_join(sample_info, by = c("sample_id", "marker_id"))

  write_rds(seq_tbl, paste(output_dir, "/seq_tbl.rds", sep = ""))

  return(seq_tbl)
}
