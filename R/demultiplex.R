#' Demultiplex paired-end FASTQ file by marker only
#'
#' Demultiplex the paired-end FASTQ reads (demultiplexed by sample) and assigns the sequence to each marker.
#'
#' @param sample_read_manifest (Required). The sample reads file. The file should include sample_id, reads_1 (Forward fastq filepaths and fastq filenames have format: SAMPLENAME_R1.fastq.gz), reads_2 (Reverse fastq filepaths and fastq filenames have format: SAMPLENAME_R2.fastq.gz), sample (sample name, can be the same as sample_id), info (e.g., sample type).
#' @param marker_info (Required). The target amplicon file. The file should include marker_id, primer_fwd, primer_rev, seq (reference sequence), chrom (chromosome), start (reference sequence start position), end (reference sequence end position).
#' @param output_dir (Required). The path to the output demultiplexed read table file.
#' @param ... Other arguments passed to AmpSeqR::demultiplex_reads.
#'
#' @return
#'
#' demultiplex folder: demultiplexed paired-end fastq files.
#'
#' demultiplex.rds: the demultiplexed table in RDS format which includes sample_id, marker_id, reads_1 (the forward read fastq file path), reads_2 (the reverse read fastq file path), n (number of demultiplexed reads), sample, info.
#'
#' @export
#'
demultiplex_marker_only <- function(sample_read_manifest,
                                    marker_info,
                                    output_dir,
                                    ...) {
  # check args
  stopifnot(
    is.data.frame(sample_read_manifest),
    is.data.frame(marker_info),
    is_string(output_dir)
  )

  check_sample_read_manifest(sample_read_manifest)
  check_marker_info(marker_info)

  ret_tbl <-
    split(sample_read_manifest, sample_read_manifest$sample_id) %>%
    map_df(function(sample_row) {
      manifest <-
        sample_row %>%
        select(-reads_1, -reads_2) %>%
        mutate(
          barcode_fwd = "",
          barcode_rev = ""
        )
      demultiplex_reads(
        sample_manifest = manifest,
        marker_info = marker_info,
        reads_1 = sample_row$reads_1,
        reads_2 = sample_row$reads_2,
        output_dir = output_dir,
        ...
      )
    }) %>%
    write_rds(file.path(output_dir, "demultiplex.rds"))

  return(ret_tbl)
}

#' Demultiplex paired-end FASTQ file
#'
#' Demultiplex the raw paired-end FASTQ reads with mixed samples and trim the sample barcodes and target amplicon primer sequences and assigns the sequence to each sample each marker.
#'
#' @param sample_manifest (Required). The sample barcodes file. The file should include sample_id, barcode_fwd, barcode_rev, sample (sample name, can be the same as sample_id), info (e.g., sample type).
#' @param marker_info (Required). The target amplicon file. The file should include marker_id, primer_fwd, primer_rev, seq (reference sequence), chrom (chromosome), start (reference sequence start position), end (reference sequence end position).
#' @param reads_1 (Required). The file path to the forward fastq file from paired-end sequence data. Compressed file formats such as .fastq.gz are supported.
#' @param reads_2 (Required). The file path to the reverse fastq file from paired-end sequence data corresponding to those provided to the reads_1 argument. Compressed file formats such as .fastq.gz are supported.
#' @param output_dir (Required). The path to the output demultiplexed read table file.
#' @param output_sub_dir (Required). The path to the output demultiplexed files.
#' @param complete_only (Optional). Default TRUE. If TRUE, only output complete fastq file.
#' @param trim_bc (Optional). Default TRUE. If TRUE, trim sample barcodes.
#' @param trim_pr (Optional). Default TRUE. If TRUE, trim amplicon primer sequence.
#' @param trim_right (Optional). Default 0. The number of nucleotides to remove from the end of each read.
#' @param chunk_n (Optional). Default 1e6. For ShortRead::FastqSampler, the size of the sample (number of records) to be drawn. For ShortRead::FastqStreamer a numeric(1) (set to 1e6 when n is missing) providing the number of successive records to be returned on each yield, or an IRanges-class delimiting the (1-based) indicies of records returned by each yield; entries in n must have non-zero width and must not overlap.
#' @param overwrite (Optional). Default TRUE. If TRUE, overwrite output files already exist.
#' @param primer_max_mismatch (Optional). Default 3. The maximum number of mismatches of the amplicon primer sequence.
#' @param max_gap_1 (Optional). Default 0. The maximum number of bases between left end of read and start of barcode sequence.
#' @param max_gap_2 (Optional). Default 0. The maximum number of bases between the barcode and the marker primer sequence.
#' @param degenerate_primers (Optional). Default TRUE. If TRUE (the default), an IUPAC ambiguity code in the pattern can only match the same code in the subject, and vice versa. If FALSE, an IUPAC ambiguity code in the pattern can match any letter in the subject that is associated with the code, and vice versa. See Biostrings::vmatchPattern for more information.
#' @param suffix_1 (Optional). Default "R1.fastq.gz". Output forward fastq file name extension.
#' @param suffix_2 (Optional). Default "R2.fastq.gz". Output reverse fastq file name extension.
#'
#' @return
#'
#' demultiplex folder: demultiplexed paired-end fastq files.
#'
#' demultiplex.rds: the demultiplexed table in RDS format which includes sample_id, marker_id, reads_1 (the forward read fastq file path), reads_2 (the reverse read fastq file path), n (number of demultiplexed reads), sample, info.
#'
#' @export
#'
#' @examples
#'
#' example_data <- get_ampseqr_example_data()
#'
#' # Set the directory
#' run_dir <- "runs"
#' dir.create(run_dir)
#'
#' demultiplexed <- demultiplex_reads(
#'   sample_manifest = example_data$sample_manifest,
#'   marker_info = example_data$marker_info,
#'   reads_1 = example_data$reads_1,
#'   reads_2 = example_data$reads_2,
#'   output_dir = run_dir,
#'   output_sub_dir = file.path(run_dir, "demultiplex")
#' )
#'
#' @importFrom rlang is_bool is_string is_scalar_integerish
#' @importFrom Biostrings width DNAStringSet
#' @importFrom dplyr select mutate left_join full_join group_by ungroup summarise filter if_else bind_rows case_when arrange pull
#' @importFrom tidyr replace_na gather chop unnest
#' @importFrom purrr map map_df map2 map_lgl walk map2_dbl
#' @importFrom parallel clusterCall clusterMap stopCluster
#' @importFrom parallelly makeClusterPSOCK
#' @importFrom stringr str_c
#' @importFrom magrittr "%>%"
#' @importFrom readr write_rds
demultiplex_reads <- function(sample_manifest,
                              marker_info,
                              reads_1,
                              reads_2,
                              output_dir,
                              output_sub_dir = file.path(output_dir, "demultiplex"),
                              complete_only = TRUE,
                              trim_bc = TRUE,
                              trim_pr = TRUE,
                              trim_right = 0L,
                              chunk_n = 1e6,
                              overwrite = TRUE,
                              primer_max_mismatch = 3L,
                              max_gap_1 = 0L,
                              max_gap_2 = 0L,
                              degenerate_primers = TRUE,
                              suffix_1 = "R1.fastq.gz",
                              suffix_2 = "R2.fastq.gz") {
  # check args
  stopifnot(
    length(reads_1) >= 1L,
    length(reads_1) == length(reads_2),
    is_scalar_character(reads_1) && file.exists(reads_1),
    is_scalar_character(reads_2) && file.exists(reads_2),
    is.data.frame(sample_manifest),
    is.data.frame(marker_info),
    is_bool(trim_bc),
    is_bool(trim_pr),
    is_scalar_integerish(trim_right) && trim_right >= 0L,
    is_bool(degenerate_primers),
    is_string(output_sub_dir),
    is_scalar_integerish(chunk_n) && chunk_n > 0L,
    is_scalar_character(suffix_1),
    is_scalar_character(suffix_2)
  )


  if (trim_pr && !trim_bc) {
    rlang::warn("trim_bc ignored as trim_pr == TRUE")
  }

  check_sample_manifest(sample_manifest)
  check_marker_info(marker_info)

  barcodes <- list(
    fwd = sample_manifest$barcode_fwd %>% unique() %>% DNAStringSet(),
    rev = sample_manifest$barcode_rev %>% unique() %>% DNAStringSet()
  )

  primers <- list(
    fwd = marker_info$primer_fwd %>% DNAStringSet(),
    rev = marker_info$primer_rev %>% DNAStringSet()
  )

  sample_manifest_2 <-
    sample_manifest %>%
    mutate(
      bc_index_1 = match(barcode_fwd, as.character(barcodes$fwd)),
      bc_index_2 = match(barcode_rev, as.character(barcodes$rev))
    ) %>%
    select(sample_id, bc_index_1, bc_index_2)

  marker_info_2 <-
    marker_info %>%
    mutate(
      pr_index_1 = match(primer_fwd, as.character(primers$fwd)),
      pr_index_2 = match(primer_rev, as.character(primers$rev))
    ) %>%
    select(marker_id, pr_index_1, pr_index_2)

  marker_trim <-
    marker_info %>%
    mutate(trim_width = nchar(seq) - trim_right) %>%
    select(marker_id, trim_width)

  # create output directory
  if (!dir.exists(output_sub_dir)) {
    dir.create(output_sub_dir, recursive = T)
  }
  output_sub_dir <- normalizePath(output_sub_dir)

  # return table
  ret_tbl <-
    tibble(
      sample_id = character(),
      marker_id = character(),
      reads_1 = character(),
      reads_2 = character(),
      n = integer()
    )

  future_rng_opt <- getOption("future.rng.onMisuse")
  options(future.rng.onMisuse = "ignore")
  on.exit(options(future.rng.onMisuse = future_rng_opt))

  workers <- list(
    future::makeClusterPSOCK(workers = 1),
    future::makeClusterPSOCK(workers = 1)
  )
  on.exit({
    walk(workers, parallel::stopCluster)
  })

  worker_args <- list(
    list(
      reads = reads_1,
      barcodes = barcodes$fwd,
      primers = primers$fwd,
      degenerate_primers = degenerate_primers,
      primer_max_mismatch = primer_max_mismatch,
      chunk_n = chunk_n,
      max_gap_1 = max_gap_1,
      max_gap_2 = max_gap_2
    ),
    list(
      reads = reads_2,
      barcodes = barcodes$rev,
      primers = primers$rev,
      degenerate_primers = degenerate_primers,
      primer_max_mismatch = primer_max_mismatch,
      chunk_n = chunk_n,
      max_gap_1 = max_gap_1,
      max_gap_2 = max_gap_2
    )
  )

  # setup workers
## ---- 1) Set up a PSOCK cluster with one node per worker_arg ----
nworkers <- length(worker_args)
stopifnot(nworkers >= 1)

cl <- parallelly::makeClusterPSOCK(nworkers)
on.exit(parallel::stopCluster(cl), add = TRUE)

## Per-worker setup (one arg list per node; keeps state on that node)
parallel::clusterMap(
  cl,
  fun = function(a) {
    suppressWarnings(do.call(AmpSeqR:::thread_setup, a))
    NULL
  },
  a = worker_args
)

## ---- 2) Processing loop with stable worker affinity ----
repeat {
  # Each node reports how many reads are queued
  nr <- unlist(parallel::clusterCall(cl, AmpSeqR:::thread_read), use.names = FALSE)
  
  # Paired-end sanity check (assumes 2 workers for R1/R2; extend as needed)
  if (length(unique(nr)) != 1L) {
    rlang::abort("number or reads in reads_1 is not equal to number of reads in reads_2")
  }
  if (nr[1] == 0L) break
  
  # Each node demultiplexes its side; returns list(dm1, dm2, ...)
  dm <- parallel::clusterCall(cl, AmpSeqR:::thread_demultiplex)
  
  # ---- Combine matched fwd/rev reads, split into sample markers ----
  dm_tbl <-
    full_join(dm[[1]], dm[[2]], by = "sr_index", suffix = c("_1", "_2")) %>%
    filter(
      replace_na(pr_end_1 < width_1, TRUE),
      replace_na(pr_end_2 < width_2, TRUE)
    ) %>%
    left_join(sample_manifest_2, by = c("bc_index_1", "bc_index_2")) %>%
    left_join(marker_info_2,    by = c("pr_index_1", "pr_index_2")) %>%
    select(sample_id, sr_index, marker_id, pr_start_1, pr_start_2, pr_end_1, pr_end_2) %>%
    chop(c(sr_index, pr_start_1, pr_start_2, pr_end_1, pr_end_2)) %>%
    mutate(
      is_complete = (!is.na(sample_id)) & (!is.na(marker_id)),
      n          = lengths(sr_index),
      prefix     = str_c("sample", replace_na(sample_id, "NA"),
                         "marker", replace_na(marker_id, "NA"), sep = "_"),
      reads_1    = file.path(output_sub_dir, str_c(prefix, "_", suffix_1)) %>%
        replace(complete_only & (!is_complete), NA_character_),
      reads_2    = file.path(output_sub_dir, str_c(prefix, "_", suffix_2)) %>%
        replace(complete_only & (!is_complete), NA_character_),
      mode       = if_else(reads_1 %in% ret_tbl$reads_1, "a", "w")
    )
  
  # ---- Check for existing files in output directory ----
  existing_files <-
    dm_tbl %>%
    filter(mode == "w") %>%
    select(reads_1, reads_2) %>%
    tidyr::pivot_longer(everything(), values_to = "value") %>%
    filter(!is.na(value), file.exists(value)) %>%
    pull(value)
  
  if (length(existing_files) > 0) {
    if (isTRUE(overwrite)) {
      invisible(file.remove(existing_files))
    } else {
      rlang::abort("some output files already exist and overwrite is set to FALSE")
    }
  }
  
  # ---- Build per-node write tables (set==1 goes to node 1, set==2 to node 2, etc.) ----
  write_table <-
    dm_tbl %>%
    filter(is_complete | (!complete_only)) %>%
    (function(x) {
      bind_rows(
        mutate(x, start = dplyr::case_when(
          trim_pr & is_complete ~ purrr::map(pr_end_1, ~ . + 1L),
          trim_bc & is_complete ~ purrr::map(pr_start_1, ~ .),
          TRUE                  ~ purrr::map(n, ~ rep(1L, .))
        )) %>%
          select(marker_id, filename = reads_1, mode, sr_index, start, pr_end = pr_end_1) %>%
          mutate(set = 1L),
        mutate(x, start = dplyr::case_when(
          trim_pr & is_complete ~ purrr::map(pr_end_2, ~ . + 1L),
          trim_bc & is_complete ~ purrr::map(pr_start_2, ~ .),
          TRUE                  ~ purrr::map(n, ~ rep(1L, .))
        )) %>%
          select(marker_id, filename = reads_2, mode, sr_index, start, pr_end = pr_end_2) %>%
          mutate(set = 2L)
      )
    }) %>%
    left_join(marker_trim, by = "marker_id") %>%
    mutate(end = map2(pr_end, trim_width, ~ .x + .y)) %>%
    select(filename, mode, sr_index, start, end, set)
  
  # Split into a list of length nworkers; pad missing sets with empty frames
  wt_split <- split(write_table, write_table$set)
  empty_df <- write_table[0, c("filename","mode","sr_index","start","end","set")]
  write_table_list <- lapply(seq_len(nworkers), function(i) wt_split[[as.character(i)]] %||% empty_df)
  
  # ---- Write on each node with its own chunk ----
  parallel::clusterMap(
    cl,
    fun = function(d) {
      AmpSeqR:::thread_write(d)
      NULL
    },
    d = write_table_list
  )
  
  # ---- Record results ----
  ret_tbl <-
    ret_tbl %>%
    bind_rows(dm_tbl %>% select(sample_id, marker_id, reads_1, reads_2, n)) %>%
    group_by(sample_id, marker_id, reads_1, reads_2) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    left_join(sample_manifest %>% select(-barcode_fwd, -barcode_rev), by = "sample_id")
}


  mutate(ret_tbl, success = !is.na(sample_id) & !is.na(marker_id)) %>%
    group_by(success) %>%
    summarise(n = sum(n, na.rm = TRUE)) %>%
    with(message(str_c("note:", n[success], "of", sum(n), "reads demultiplexed successfully.", sep = " ")))

  write_rds(ret_tbl, file.path(output_dir, "demultiplex.rds"))

  return(ret_tbl)
}


thread_setup <- function(reads, chunk_n, primers, barcodes, degenerate_primers, primer_max_mismatch, max_gap_1, max_gap_2) {
  # set arguments as options, as globals are overwritten by future::, and parallel:: doesn't support async
  options(
    ampseqr.fs = ShortRead::FastqStreamer(reads, n = chunk_n),
    ampseqr.primers = primers,
    ampseqr.barcodes = barcodes,
    ampseqr.degenerate_primers = degenerate_primers,
    ampseqr.primer_max_mismatch = primer_max_mismatch,
    ampseqr.max_gap_1 = max_gap_1,
    ampseqr.max_gap_2 = max_gap_2
  )
  return(TRUE)
}

thread_read <- function() {
  stopifnot(!is.null(getOption("ampseqr.fs")))
  options(ampseqr.sr = ShortRead::yield(getOption("ampseqr.fs")))
  return(length(getOption("ampseqr.sr")))
}

thread_demultiplex <- function() {
  stopifnot(
    !is.null(getOption("ampseqr.sr")),
    !is.null(getOption("ampseqr.primers")),
    !is.null(getOption("ampseqr.barcodes")),
    !is.null(getOption("ampseqr.degenerate_primers")),
    !is.null(getOption("ampseqr.primer_max_mismatch")),
    !is.null(getOption("ampseqr.max_gap_1")),
    !is.null(getOption("ampseqr.max_gap_2"))
  )

  dm <- match_barcode_primer(
    dss = ShortRead::sread(getOption("ampseqr.sr")),
    barcodes = getOption("ampseqr.barcodes"),
    primers = getOption("ampseqr.primers"),
    fixed = getOption("ampseqr.degenerate_primers"),
    max_mismatch = getOption("ampseqr.primer_max_mismatch"),
    max_gap_1 = getOption("ampseqr.max_gap_1"),
    max_gap_2 = getOption("ampseqr.max_gap_2")
  )
  return(dm)
}

thread_write <- function(write_tbl) {
  stopifnot(!is.null(getOption("ampseqr.sr")))
  sr <- getOption("ampseqr.sr")
  pwalk_write_reads(sr, write_tbl)
}

# return dataframe with columns sr_index, bc_index, pr_index, pr_start, pr_end, width
#' @importFrom Biostrings vcountPattern vmatchPattern width endIndex
#' @importFrom ShortRead narrow
#' @importFrom purrr map_df
#' @importFrom tibble tibble as_tibble
# Fast version: vectorized, low-allocation, Biostrings-first
match_barcode_primer <- function(dss, barcodes, primers, fixed, max_mismatch,
                                 max_gap_1, max_gap_2, offset = 0L) {
  # TODO: allow non left anchored barcodes and primers

  bc_width <- width(barcodes[1])

  res <- tibble(
    sr_index = seq_along(dss),
    bc_index = NA_integer_,
    bc_start = NA_integer_,
    bc_end = NA_integer_,
    pr_index = NA_integer_,
    pr_start = NA_integer_,
    pr_end = NA_integer_,
    width = width(dss)
  )

  if (bc_width > 0) {
    # match barcodes
    sr_index <- res$sr_index[which(width(dss) >= bc_width)]
    sub <- narrow(dss[sr_index], start = 1, end = pmin(bc_width + max_gap_1, width(dss[sr_index])))
    for (i in seq_along(barcodes)) {
      match <- which(vcountPattern(barcodes[[i]], sub, fixed = fixed) > 0)
      if (length(match) > 0) {
        match_coord <-
          Biostrings::vmatchPattern(barcodes[[i]], sub[match], fixed = fixed) %>%
          (function(x) {
            tibble(
              start = Biostrings::startIndex(x) %>% map_int(first),
              end = Biostrings::endIndex(x) %>% map_int(first),
            )
          })
        res$bc_index[sr_index[match]] <- i
        res$bc_start[sr_index[match]] <- match_coord$start
        res$bc_end[sr_index[match]] <- match_coord$end
        sub <- sub[-match]
        sr_index <- sr_index[-match]
        if (length(sr_index) == 0L) {
          break
        }
      }
    }
  } else {
    # empty barcode matches all (i.e. unbarcoded sample)
    res$bc_index <- 1L
    res$bc_start <- 0L
    res$bc_end <- 0L
  }

  # match primers
  for (i in seq_along(primers)) {
    pr_width <- width(primers[i])
    sub_range <-
      filter(res, is.na(pr_index)) %>%
      mutate(
        start = if_else(!is.na(bc_end),
                        bc_end + 1L,
                        bc_width + 1L
        ),
        end = if_else(!is.na(bc_end),
                      start + pr_width + max_gap_2,
                      start + pr_width + max_gap_1 + max_gap_2
        ),
        end = pmin(end, width)
      ) %>%
      filter(start + pr_width <= width) %>%
      select(sr_index, start, end)
    if (nrow(sub_range) == 0) {
      next
    }
    sr_index <- sub_range$sr_index
    sub <- narrow(dss[sr_index], start = sub_range$start, end = sub_range$end)
    match <- which(vcountPattern(primers[[i]], sub, fixed = fixed, max.mismatch = max_mismatch) > 0)
    if (length(match) == 0) {
      next
    }
    match_coord <-
      Biostrings::vmatchPattern(primers[[i]], sub[match], fixed = fixed, max.mismatch = max_mismatch) %>%
      (function(x) {
        tibble(
          start = Biostrings::startIndex(x) %>% map_int(first),
          end = Biostrings::endIndex(x) %>% map_int(first),
        )
      })
    res$pr_index[sr_index[match]] <- i
    res$pr_start[sr_index[match]] <- match_coord$start + sub_range$start[match] - 1L
    res$pr_end[sr_index[match]] <- match_coord$end + sub_range$start[match] - 1L
  }
  res$sr_index <- res$sr_index + offset

  return(res)
}
# match_barcode_primer <- function(dss, barcodes, primers, fixed, max_mismatch,
#                                  max_gap_1, max_gap_2, offset = 0L) {
#   # Preconditions & helpers
#   stopifnot(inherits(dss, "DNAStringSet"),
#             inherits(barcodes, "DNAStringSet"),
#             inherits(primers, "DNAStringSet"))
#   n <- length(dss)
#   w <- Biostrings::width(dss)
#   bc_width <- if (length(barcodes)) Biostrings::width(barcodes[1]) else 0L
#   
#   # Result vectors (preallocated; much faster than building tibbles and mutate)
#   sr_index <- seq_len(n)
#   bc_index <- rep(NA_integer_, n)
#   bc_start <- rep(NA_integer_, n)
#   bc_end   <- rep(NA_integer_, n)
#   pr_index <- rep(NA_integer_, n)
#   pr_start <- rep(NA_integer_, n)
#   pr_end   <- rep(NA_integer_, n)
#   
#   ## -------------------------
#   ## 1) BARCODE MATCH (exact, left-anchored window)
#   ## -------------------------
#   if (bc_width > 0L) {
#     # subset reads with enough length for barcode window
#     ok_bc <- which(w >= bc_width)
#     if (length(ok_bc)) {
#       # Views covering [1, bc_width + max_gap_1] for each ok read
#       bc_end_cap <- pmin.int(bc_width + max_gap_1, w[ok_bc])
#       bc_views <- Biostrings::Views(subject = dss[ok_bc], start = rep.int(1L, length(ok_bc)), end = bc_end_cap)
#       
#       # For exact matches, use a PDict; this vectorizes over barcodes
#       # (fast path only applies if 'fixed' implies exact matching)
#       # If 'fixed' is TRUE or c("pattern","subject") with exact semantics, this is safe.
#       # NOTE: barcodes are treated exact as in your original code.
#       pdict_bc <- Biostrings::PDict(barcodes)
#       # Count matrix: patterns x subjects (barcodes x reads)
#       # This is very fast and avoids per-pattern loops.
#       cm <- Biostrings::vcountPDict(pdict_bc, bc_views, max.mismatch = 0L, fixed = fixed)
#       
#       # First matching barcode per read (choose the first pattern with count > 0)
#       # We prefer the first barcode in the provided order, consistent with your loop.
#       has_hit <- colSums(cm > 0L) > 0L
#       if (any(has_hit)) {
#         hit_cols <- which(has_hit)
#         # vector of first barcode index per matched read
#         first_bc <- apply(cm[, hit_cols, drop = FALSE], 2L, function(z) {
#           ii <- which(z > 0L)
#           if (length(ii)) ii[1L] else NA_integer_
#         })
#         
#         # For coordinates, use vmatchPDict to get positions but only for matched barcodes
#         # This returns an MIndex list per barcode pattern over all subjects (views).
#         # We'll extract the first hit per matched subject.
#         mi <- Biostrings::vmatchPDict(pdict_bc, bc_views, max.mismatch = 0L, fixed = fixed)
#         # To avoid scanning everything, we only touch subjects that had hits.
#         for (j in seq_along(hit_cols)) {
#           col_j <- hit_cols[j]                 # subject index within ok_bc/views
#           pat_i <- first_bc[j]                 # chosen barcode index
#           if (is.na(pat_i)) next
#           hits_j <- mi[[pat_i]][[col_j]]
#           # Start/endIndex give IntegerList per subject; we take the first (your original code did the same)
#           s_idx <- Biostrings::startIndex(hits_j)
#           e_idx <- Biostrings::endIndex(hits_j)
#           if (length(s_idx)) {
#             ridx <- ok_bc[col_j]              # row in original dss
#             bc_index[ridx] <- pat_i
#             bc_start[ridx] <- s_idx[1L]
#             bc_end[ridx]   <- e_idx[1L]
#           }
#         }
#       }
#     }
#   } else {
#     # Empty barcode = match all
#     bc_index[] <- 1L
#     bc_start[] <- 0L
#     bc_end[]   <- 0L
#   }
#   
#   ## -------------------------
#   ## 2) PRIMER MATCH (may allow mismatches)
#   ## -------------------------
#   # Determine per-read search windows for the primer search.
#   # If barcode was found, start after bc_end; else start after bc_width.
#   # end = start + primer_width + max_gap_2 (+ max_gap_1 if no BC)
#   # We do per-primer passes to use vectorized vcountPattern/vmatchPattern over the candidate views.
#   # Avoid tibble/dplyr entirely.
#   # For each primer i, compute candidate subranges only for reads with no primer yet.
#   for (i in seq_along(primers)) {
#     # Skip reads already assigned a primer
#     todo <- which(is.na(pr_index))
#     if (!length(todo)) break
#     
#     pr_w <- Biostrings::width(primers[i])
#     
#     # compute starts/ends
#     # start_if_bc = bc_end + 1; if bc_end is NA, we use bc_width + 1
#     has_bc <- !is.na(bc_end[todo])
#     start_vec <- integer(length(todo))
#     start_vec[has_bc]  <- bc_end[todo][has_bc] + 1L
#     start_vec[!has_bc] <- bc_width + 1L
#     
#     # per original logic: end depends on whether bc_end is known
#     end_vec <- integer(length(todo))
#     end_vec[has_bc]  <- start_vec[has_bc] + pr_w + max_gap_2
#     end_vec[!has_bc] <- start_vec[!has_bc] + pr_w + max_gap_1 + max_gap_2
#     
#     # clamp to read width
#     # also require start + pr_w <= width to be worth searching
#     end_vec <- pmin.int(end_vec, w[todo])
#     good <- which(start_vec + pr_w <= w[todo] & start_vec <= end_vec)
#     if (!length(good)) next
#     
#     idx <- todo[good]
#     # Build views for candidate regions
#     pv <- Biostrings::Views(dss[idx], start = start_vec[good], end = end_vec[good])
#     
#     # Vectorized count across candidate views
#     vc <- Biostrings::vcountPattern(primers[[i]], pv,
#                                     fixed = fixed, max.mismatch = max_mismatch)
#     hit <- which(vc > 0L)
#     if (!length(hit)) next
#     
#     # Get coordinates for the first hit in each candidate view
#     vmi <- Biostrings::vmatchPattern(primers[[i]], pv[hit],
#                                      fixed = fixed, max.mismatch = max_mismatch)
#     # Extract first start/end per subject
#     # startIndex/endIndex return IntegerList; we take the first element
#     s_list <- Biostrings::startIndex(vmi)
#     e_list <- Biostrings::endIndex(vmi)
#     
#     # Assign back to global vectors
#     for (k in seq_along(hit)) {
#       h <- hit[k]
#       ridx <- idx[h]  # original read index
#       s1 <- s_list[[k]][1L]
#       e1 <- e_list[[k]][1L]
#       if (!is.na(s1) && !is.na(e1)) {
#         pr_index[ridx] <- i
#         # absolute coordinates: add subrange start - 1
#         pr_start[ridx] <- s1 + start_vec[good][h] - 1L
#         pr_end[ridx]   <- e1 + start_vec[good][h] - 1L
#       }
#     }
#   }
#   
#   ## -------------------------
#   ## Build result tibble
#   ## -------------------------
#   tibble::tibble(
#     sr_index = sr_index + offset,
#     bc_index = bc_index,
#     bc_start = bc_start,
#     bc_end   = bc_end,
#     pr_index = pr_index,
#     pr_start = pr_start,
#     pr_end   = pr_end,
#     width    = w
#   )
# }
# match_barcode_primer <- function(dss, barcodes, primers, fixed, max_mismatch,
#                                  max_gap_1, max_gap_2, offset = 0L) {
#   # TODO: allow non left anchored barcodes and primers
#   
#   bc_width <- width(barcodes[1])
#   
#   res <- tibble(
#     sr_index = seq_along(dss),
#     bc_index = NA_integer_,
#     bc_start = NA_integer_,
#     bc_end = NA_integer_,
#     pr_index = NA_integer_,
#     pr_start = NA_integer_,
#     pr_end = NA_integer_,
#     width = width(dss)
#   )
#   
#   if (bc_width > 0) {
#     # match barcodes
#     sr_index <- res$sr_index[which(width(dss) >= bc_width)]
#     sub <- narrow(dss[sr_index], start = 1, end = pmin(bc_width + max_gap_1, width(dss[sr_index])))
#     for (i in seq_along(barcodes)) {
#       match <- which(vcountPattern(barcodes[[i]], sub, fixed = fixed) > 0)
#       if (length(match) > 0) {
#         match_coord <-
#           Biostrings::vmatchPattern(barcodes[[i]], sub[match], fixed = fixed) %>%
#           (function(x) {
#             tibble(
#               start = Biostrings::startIndex(x) %>% map_int(first),
#               end = Biostrings::endIndex(x) %>% map_int(first),
#             )
#           })
#         res$bc_index[sr_index[match]] <- i
#         res$bc_start[sr_index[match]] <- match_coord$start
#         res$bc_end[sr_index[match]] <- match_coord$end
#         sub <- sub[-match]
#         sr_index <- sr_index[-match]
#         if (length(sr_index) == 0L) {
#           break
#         }
#       }
#     }
#   } else {
#     # empty barcode matches all (i.e. unbarcoded sample)
#     res$bc_index <- 1L
#     res$bc_start <- 0L
#     res$bc_end <- 0L
#   }
#   
#   # match primers
#   for (i in seq_along(primers)) {
#     pr_width <- width(primers[i])
#     sub_range <-
#       filter(res, is.na(pr_index)) %>%
#       mutate(
#         start = if_else(!is.na(bc_end),
#                         bc_end + 1L,
#                         bc_width + 1L
#         ),
#         end = if_else(!is.na(bc_end),
#                       start + pr_width + max_gap_2,
#                       start + pr_width + max_gap_1 + max_gap_2
#         ),
#         end = pmin(end, width)
#       ) %>%
#       filter(start + pr_width <= width) %>%
#       select(sr_index, start, end)
#     if (nrow(sub_range) == 0) {
#       next
#     }
#     sr_index <- sub_range$sr_index
#     sub <- narrow(dss[sr_index], start = sub_range$start, end = sub_range$end)
#     match <- which(vcountPattern(primers[[i]], sub, fixed = fixed, max.mismatch = max_mismatch) > 0)
#     if (length(match) == 0) {
#       next
#     }
#     match_coord <-
#       Biostrings::vmatchPattern(primers[[i]], sub[match], fixed = fixed, max.mismatch = max_mismatch) %>%
#       (function(x) {
#         tibble(
#           start = Biostrings::startIndex(x) %>% map_int(first),
#           end = Biostrings::endIndex(x) %>% map_int(first),
#         )
#       })
#     res$pr_index[sr_index[match]] <- i
#     res$pr_start[sr_index[match]] <- match_coord$start + sub_range$start[match] - 1L
#     res$pr_end[sr_index[match]] <- match_coord$end + sub_range$start[match] - 1L
#   }
#   res$sr_index <- res$sr_index + offset
#   
#   return(res)
# }
#' @importFrom purrr pwalk
pwalk_write_reads <- function(reads, data) {
  data %>%
    pwalk(function(filename, mode, sr_index, start, end, ...) {
      write_reads(
        reads = reads[sr_index],
        filename = filename,
        mode = mode,
        start = start,
        end = end
      )
    })

  return(TRUE)
}

#' @importFrom ShortRead writeFastq narrow
write_reads <- function(reads, filename, mode, start = NULL, end = NULL) {
  if (is.null(start)) {
    start <- rep(1L, length(reads))
  }
  start <- replace(start, !is.finite(start), 1L)
  if (is.null(end)) {
    end <- width(reads)
  } else {
    end <- pmin(width(reads), end)
    end <- replace(end, !is.finite(end), width(reads)[!is.finite(end)])
  }

  stopifnot(
    length(start) == length(reads),
    length(start) == length(end)
  )

  writeFastq(narrow(reads, start = start, end = end), file = filename, mode = mode)
}
