

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
        mutate(barcode_fwd = '',
               barcode_rev = '')
      demultiplex_reads(sample_manifest = manifest,
                        marker_info = marker_info,
                        reads_1 = sample_row$reads_1,
                        reads_2 = sample_row$reads_2,
                        output_dir = output_dir,
                        ...)
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
#' @param use_absolute_paths (Required). Default to TRUE. If FALSE, paths to input and output .fastq.gz files will remain relative to the current working directory (useful in case the output is shared with other users so they can possibly continue running the pipeline without running into issues related to different storage mounts).
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
#' @importFrom purrr map map_df map2 map_lgl walk
#' @importFrom stringr str_c
#' @importFrom magrittr "%>%"
#' @importFrom readr write_rds
demultiplex_reads <- function(sample_manifest,
                              marker_info,
                              reads_1,
                              reads_2,
                              output_dir,
                              output_sub_dir = file.path(output_dir, 'demultiplex'),
                              use_absolute_paths = TRUE,
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
  
  if (!use_absolute_paths) {
    output_sub_dir <- gsub(pattern = normalizePath(getwd()), replacement = ".", x = output_sub_dir, fixed = TRUE)
  }
  
  # return table
  ret_tbl <-
    tibble(
      sample_id = character(),
      marker_id = character(),
      reads_1 = character(),
      reads_2 = character(),
      n = integer()
    )
  
  future_rng_opt <- getOption('future.rng.onMisuse')
  options(future.rng.onMisuse = "ignore")
  on.exit(options(future.rng.onMisuse = future_rng_opt))
  
  workers <- list(
    future::makeClusterPSOCK(workers = floor(future::availableCores()/2)),
    future::makeClusterPSOCK(workers = 1)
  )
  on.exit({ walk(workers, parallel::stopCluster) })
  
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
  map2(workers, worker_args, function(w, a) {
    future::cluster(
      {
        suppressWarnings(do.call(AmpSeqR:::thread_setup, a))
      },
      workers = w, 
      globals = structure(TRUE, add = list(a = a, w = w))
    )
  }) %>%
    future::value() %>%
    invisible()
  
  while (TRUE) {
    nr <-
      map(workers, function(w) {
        future::cluster(
          {
            AmpSeqR:::thread_read()
          },
          workers = w, 
          globals = structure(TRUE, add = list(w = w))
        )
      }) %>%
      map_dbl(future::value)
    
    if (nr[1] != nr[2]) {
      walk(workers, parallel::stopCluster)
      rlang::abort("number or reads in reads_1 is not equal to number of reads in reads_2")
    }
    
    if (nr[1] == 0) {
      break
    }
    
    dm <-
      map(workers, function(w) {
        future::cluster(
          {
            AmpSeqR:::thread_demultiplex()
          },
          workers = w, 
          globals = structure(TRUE, add = list(w = w))
        )
      }) %>%
      future::value()
    
    # combine matched fwd and rev reads, split into sample markers
    dm_tbl <-
      full_join(dm[[1]], dm[[2]], by = "sr_index", suffix = c("_1", "_2")) %>%
      filter(
        replace_na(pr_end_1 < width_1, TRUE),
        replace_na(pr_end_2 < width_2, TRUE)
      ) %>%
      left_join(sample_manifest_2, by = c("bc_index_1", "bc_index_2")) %>%
      left_join(marker_info_2, by = c("pr_index_1", "pr_index_2")) %>%
      select(sample_id, sr_index, marker_id, pr_start_1, pr_start_2, pr_end_1, pr_end_2) %>%
      chop(c(sr_index, pr_start_1, pr_start_2, pr_end_1, pr_end_2)) %>%
      mutate(
        is_complete = (!is.na(sample_id)) & (!is.na(marker_id)),
        n = lengths(sr_index),
        prefix = str_c("sample", replace_na(sample_id, "NA"), "marker", replace_na(marker_id, "NA"), sep = "_"),
        reads_1 = file.path(output_sub_dir, str_c(prefix, "_", suffix_1)) %>% replace(complete_only & (!is_complete), NA_character_),
        reads_2 = file.path(output_sub_dir, str_c(prefix, "_", suffix_2)) %>% replace(complete_only & (!is_complete), NA_character_),
        mode = if_else(reads_1 %in% ret_tbl$reads_1, "a", "w")
      )
    
    # check for existing files in output directory
    existing_files <-
      dm_tbl %>%
      filter(mode == "w") %>%
      select(reads_1, reads_2) %>%
      na.omit() %>%
      gather() %>%
      filter(file.exists(value)) %>%
      pull(value)
    
    if (length(existing_files) > 0) {
      if (overwrite) {
        invisible(file.remove(existing_files))
      } else {
        rlang::abort("some output files already exist and overwrite is set to FALSE")
      }
    }
    
    write_table <-
      dm_tbl %>%
      filter(is_complete | (!complete_only)) %>%
      (function(x) {
        bind_rows(
          mutate(x, start = case_when(
            trim_pr & is_complete ~ map(pr_end_1, ~ . + 1),
            trim_bc & is_complete ~ map(pr_start_1, ~.),
            TRUE ~ map(n, ~ rep(1L, .))
          )) %>%
            select(marker_id, filename = reads_1, mode, sr_index, start, pr_end = pr_end_1) %>%
            mutate(set = 1L),
          mutate(x, start = case_when(
            trim_pr & is_complete ~ map(pr_end_2, ~ . + 1),
            trim_bc & is_complete ~ map(pr_start_2, ~.),
            TRUE ~ map(n, ~ rep(1L, .))
          )) %>%
            select(marker_id, filename = reads_2, mode, sr_index, start, pr_end = pr_end_2) %>%
            mutate(set = 2L),
        )
      }) %>%
      left_join(marker_trim, "marker_id") %>%
      mutate(end = map2(pr_end, trim_width, ~ .x + .y)) %>%
      select(filename, mode, sr_index, start, end, set) %>%
      split.data.frame(.$set)
    
    # write output in each thread
    map2(workers, write_table, function(w, d) {
      future::cluster(
        {
          AmpSeqR:::thread_write(d)
        },
        workers = w, 
        globals = structure(TRUE, add = list(d = d, w = w))
      )
    }) %>%
      future::value() %>%
      invisible()
    
    # record results
    ret_tbl <-
      ret_tbl %>%
      bind_rows(dm_tbl %>% select(sample_id, marker_id, reads_1, reads_2, n)) %>%
      group_by(sample_id, marker_id, reads_1, reads_2) %>%
      summarise(n = sum(n), .groups = "drop") %>%
      left_join(sample_manifest %>% select(-barcode_fwd, -barcode_rev),
                by = "sample_id")
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
          (function(x) tibble(
            start = Biostrings::startIndex(x) %>% map_int(first),
            end = Biostrings::endIndex(x) %>% map_int(first),
          ))
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
      mutate(start = if_else(!is.na(bc_end),
                             bc_end + 1L,
                             bc_width + 1L),
             end = if_else(!is.na(bc_end),
                           start + pr_width + max_gap_2,
                           start + pr_width + max_gap_1 + max_gap_2),
             end = pmin(end, width)) %>%
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
      (function(x) tibble(
        start = Biostrings::startIndex(x) %>% map_int(first),
        end = Biostrings::endIndex(x) %>% map_int(first),
      ))
    res$pr_index[sr_index[match]] <- i
    res$pr_start[sr_index[match]] <- match_coord$start + sub_range$start[match] - 1L
    res$pr_end[sr_index[match]] <- match_coord$end + sub_range$start[match] - 1L
  }
  res$sr_index <- res$sr_index + offset
  
  return(res)
}

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
