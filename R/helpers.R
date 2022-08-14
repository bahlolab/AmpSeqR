
format_int <- function(x, min_width = 1L, na_str = NA_character_) {
  stopifnot(rlang::is_integerish(x))

  width <- max(min_width, 1 + floor(log10(max(x, na.rm = T))))
  stringr::str_pad(x, pad = "0", side = "left", width = width) %>%
    tidyr::replace_na(na_str)
}


#' @importFrom rlang abort
#' @importFrom purrr map_lgl map2_lgl
#' @importFrom stringr str_c
check_table <- function(table, col_types, arg_name) {
  stopifnot(
    is.list(col_types),
    !is.null(names(col_types)),
    all(map_lgl(col_types, is_scalar_character)),
    is_scalar_character(arg_name)
  )

  if (!is.data.frame(table)) {
    abort(str_c(arg_name, " must be a data.frame or tibble"))
  }

  if (!all(names(col_types) %in% colnames(table))) {
    abort(str_c(arg_name, " is missing column(s): ", str_c(setdiff(names(col_types), colnames(table)), collapse = ", ")))
  }

  if (nrow(table) == 0) {
    abort(str_c(arg_name, " has zero rows"))
  }

  mismatch <- map2_lgl(names(col_types), col_types, function(col, type) {
    `if`(type == "integerish", !is_integerish(table[[col]]), !is(table[[col]], type))
  })

  if (any(mismatch)) {
    abort(str_c(
      "incorrect column types in ", arg_name, ": ",
      str_c(names(col_types)[mismatch], " (expected ", unlist(col_types[mismatch]), ")", collapse = ", ")
    ))
  }
}


#' @importFrom rlang abort
check_character_table <- function(table, req_cols, arg_name) {
  if (!is.data.frame(table)) {
    abort(str_c(arg_name, " must be a data.frame"))
  }

  if (!all(req_cols %in% colnames(table))) {
    abort(str_c(arg_name, " must contain columns ", str_c(req_cols, collapse = ", ")))
  }

  if (nrow(table) == 0) {
    abort(str_c(arg_name, " has zero rows"))
  }

  for (cn in req_cols) {
    if (!is.character(table[[cn]])) {
      abort(str_c(str_c(arg_name, "$", cn), " must be a character vector"))
    }
  }
}


check_seq_table <- function(seq_tbl, seq_id = FALSE) {
  col_types <- list(
    sample_id = "character",
    marker_id = "character",
    sequence = "character",
    count = "integer"
  )

  if ("status" %in% names(seq_tbl)) {
    col_types[["status"]] <- "character"
  }

  if (seq_id) {
    col_types[["seq_id"]] <- "character"
  }

  check_table(seq_tbl, col_types, "seq_tbl")
  check_is_dna_char(seq_tbl$sequence, "seq_tbl$sequence", na.rm = TRUE, IUPAC = TRUE)

  dups <- seq_tbl %>%
    group_by(sample_id, marker_id, sequence) %>%
    count() %>%
    filter(n > 1)
  if (nrow(dups)) {
    abort(str_c("duplicate entries detected for sample_id, marker_id and seq_hash in seq_tbl"))
  }
}

check_var_table <- function(var_tbl) {
  col_types <- list(
    marker_id = "character",
    var_id = "character",
    pos = "integer",
    aln_pos = "integer",
    ref = "character",
    alt = "list",
    gt_data = "list"
  )

  check_table(var_tbl, col_types, "var_tbl")
  check_is_dna_char(var_tbl$ref, "var_tbl$ref", na.rm = FALSE, IUPAC = TRUE)
  check_is_dna_char(unlist(var_tbl$alt), "var_tbl$alt", na.rm = FALSE, IUPAC = TRUE)
}

#' @importFrom rlang is_integerish
check_asv_table <- function(asv_tbl,
                            check_count = FALSE,
                            check_status = FALSE) {
  col_types <- list(
    seq_id = "character",
    marker_id = "character",
    sequence = "character"
  )

  if (check_count) {
    col_types <- c(col_types, list(count = "integer"))
  }
  if (check_status) {
    col_types <- c(col_types, list(status = "character"))
  }

  check_table(asv_tbl, col_types, "asv_tbl")
}

#' @importFrom stringr str_detect
check_is_dna_char <- function(char, arg_name, IUPAC = FALSE, na.rm = FALSE) {
  if (!is.character(char)) {
    abort(str_c(arg_name, " must be a character vector"))
  }

  if (na.rm) char <- na.omit(char)

  if (IUPAC) {
    if (any(str_detect(char, "[^ACGTMRWSYKVHDBN]"))) {
      abort(str_c(arg_name, " must only contain IUPAC nucleic acid code characters [ACGTMRWSYKVHDBN]"))
    }
  } else {
    if (any(str_detect(char, "[^ACTG]"))) {
      abort(str_c(arg_name, " must only contain nucleic acid code characters [ACTG]"))
    }
  }
}

check_files_exist <- function(filenames, arg_name, na.rm = FALSE) {
  if (!is.character(filenames)) {
    abort(str_c(arg_name, " must be a character vector"))
  }

  if (na.rm) {
    filenames <- na.omit(filenames)
  }

  if (!all(file.exists(filenames))) {
    abort(str_c("some files in ", arg_name, " do not exist"))
  }
}

#' @importFrom dplyr count pull group_by select distinct
#' @importFrom rlang abort
#' @importFrom stringr str_c
check_sample_manifest <- function(table) {
  req_cols <- c("sample_id", "barcode_fwd", "barcode_rev")
  arg_name <- "sample_manifest"
  check_character_table(table, req_cols, arg_name)

  nmax <-
    select(table, "sample_id", "barcode_fwd", "barcode_rev") %>%
    distinct() %>%
    group_by(barcode_fwd, barcode_rev) %>%
    count() %>%
    pull(n) %>%
    max(na.rm = T)

  if (nmax > 1) {
    abort(str_c(arg_name, " distinct barcode pairs must be associated with a single sample_id"))
  }

  check_is_dna_char(table$barcode_fwd, "barcode_fwd")
  check_is_dna_char(table$barcode_rev, "barcode_rev")

  if (!all(nchar(table$barcode_fwd) == nchar(table$barcode_fwd[1]))) {
    abort("all barcodes sample_manifest$barcode_fwd must be the same length")
  }
  if (!all(nchar(table$barcode_rev) == nchar(table$barcode_rev[1]))) {
    abort("all barcodes sample_manifest$barcode_rev must be the same length")
  }
}

#' @importFrom dplyr count pull group_by select distinct
#' @importFrom rlang abort
#' @importFrom stringr str_c
check_marker_info <- function(table, locus = FALSE) {
  req_cols <- c("marker_id", "primer_fwd", "primer_rev", "seq")
  arg_name <- "marker_info"
  check_character_table(table, req_cols, arg_name)

  nmax <- count(table, marker_id) %>%
    pull(n) %>%
    max(na.rm = T)

  if (nmax > 1) {
    abort(str_c(arg_name, "all entries in marker_info$marker_id must be unique"))
  }

  check_is_dna_char(table$primer_fwd, "primer_fwd", IUPAC = TRUE)
  check_is_dna_char(table$primer_rev, "primer_rev", IUPAC = TRUE)
  check_is_dna_char(table$seq, "seq")

  if (locus) {
    check_table(table, list(chrom = "character", start = "integer"), "marker_info")
  }
}

check_is_distinct <- function(x, arg_name) {
  if (length(x) != length(unique(x))) {
    abort(str_c(arg_name, " must all be unique"))
  }
}

#' @importFrom purrr map2_chr
#' @importFrom stringr str_split
add_status <- function(tbl, new_status) {
  stopifnot(
    is.data.frame(tbl),
    is.character(new_status),
    nrow(tbl) == length(new_status)
  )

  if ("status" %in% colnames(tbl)) {
    tbl %>%
      mutate(status = case_when(
        status == "pass" & new_status == "pass" ~ "pass",
        status == "pass" | is.na(status) ~ new_status,
        new_status == "pass" | is.na(new_status) ~ status,
        TRUE ~ map2_chr(status, new_status, function(s, ns) {
          c(str_split(s, "; ", simplify = TRUE), ns) %>%
            str_c(collapse = "; ")
        })
      ))
  } else {
    tbl %>% mutate(status = new_status)
  }
}


#' Downsample read counts
#'
#' This function downsamples an exact number of reads from paired end fastq files.
#'
#' @param read_table (Required). The filtered read table includes: sample_id, marker_id, n (number of demultiplexed reads), sample, info, reads_1 (the filtered and trimmed forward fastq file path), reads_2 (the filtered and trimmed reverse fastq file path), n_in (number of reads before filtering and trimming), n_out (number of reads after filtering and trimming).
#' @param output_dir (Required). The path to the output downsampled read table file.
#' @param output_sub_dir (Required). The path to the output downsampled files.
#' @param min_read_count (Optional). Default 1000. The minimum number of reads per sample per marker.
#' @param n_sample (Optional). Default 10000. Downsamples an exact number of reads from paired end fastq files.
#' @param seed (Optional). Default 1. Random seed for reproducible downsampling.
#' @param threads (Optional). Default 1. The hostnames of workers (as a character vector) or the number of localhost workers (as a positive integer).
#' @param count_col (Optional). Default NULL. The specified column to downsample reads.
#'
#' @return
#'
#' A downsampled table in RDS format includes sample_id, marker_id, sample, info, reads_1 (the downsampled forward fastq file path), reads_2 (the downsampled reverse fastq file path), n (number of reads after downsampling).
#'
#' @export
#' @importFrom purrr pmap
downsample_reads <- function(read_table,
                             output_dir,
                             output_sub_dir,
                             min_read_count = 1000,
                             n_sample = 10000,
                             seed = 1L,
                             threads = 1L,
                             count_col = NULL) {
  stopifnot(
    is_scalar_character(output_dir),
    is_scalar_character(output_sub_dir),
    is.null(count_col) | is_scalar_character(count_col),
    is_scalar_integerish(seed) && seed >= 0L,
    is_scalar_integerish(threads) && threads >= 1L,
    is_scalar_integerish(n_sample) && n_sample >= 1L
  )

  # change
  read_table <- read_table %>%
    filter(n_out >= min_read_count) %>%
    select(sample_id, marker_id, sample, info, reads_1, reads_2, n_out)

  check_character_table(read_table, c("reads_1", "reads_2"), "read_table")
  if (!is.null(count_col)) {
    check_table(read_table, setNames(list("integerish"), count_col), "read_table")
  }

  if (!is.null(count_col)) {
    ss_at <- which(pull(read_table, !!count_col) > n_sample)
  } else {
    ss <- seq_len(nrow(read_table))
  }

  if (length(ss_at)) {
    if (!dir.exists(output_sub_dir)) {
      dir.create(output_sub_dir, recursive = T)
    }
    output_sub_dir <- normalizePath(output_sub_dir)

    cluster <- `if`(
      threads > 1,
      future::makeClusterPSOCK(workers = threads),
      NULL
    )

    ss <-
      read_table %>%
      dplyr::slice(ss_at) %>%
      select(reads_1, reads_2) %>%
      mutate(
        reads_1_sub = file.path(output_sub_dir, basename(reads_1)),
        reads_2_sub = file.path(output_sub_dir, basename(reads_2))
      )

    fs <-
      ss %>%
      pmap(function(reads_1, reads_1_sub, reads_2, reads_2_sub, ...) {
        future::cluster(
          {
            ampseqr:::downsample(
              reads_1 = reads_1,
              reads_2 = reads_2,
              reads_1_sub = reads_1_sub,
              reads_2_sub = reads_2_sub,
              n_sample = n_sample,
              seed = seed
            )
          },
          workers = cluster
        )
      }) %>%
      future::value()

    read_table$reads_1[ss_at] <- ss$reads_1_sub
    read_table$reads_2[ss_at] <- ss$reads_2_sub

    if (threads > 1) {
      parallel::stopCluster(cluster)
    }
  }

  read_table <- read_table %>%
    mutate(n = pmin(n_out, n_sample))

  write_rds(read_table, paste(output_dir, "/subsampled_reads.rds", sep = ""))

  return(read_table)
}

downsample <- function(reads_1, reads_2, reads_1_sub, reads_2_sub, n_sample, seed) {
  sr1 <- ShortRead::readFastq(reads_1)
  sr2 <- ShortRead::readFastq(reads_2)
  n <- min(length(sr1), length(sr2))
  set.seed(seed)
  if (n > n_sample) {
    i <- sample(n, n_sample)
  } else {
    i <- sample(n, n)
  }
  if (file.exists(reads_1_sub)) {
    file.remove(reads_1_sub)
  }
  if (file.exists(reads_2_sub)) {
    file.remove(reads_2_sub)
  }
  ShortRead::writeFastq(sr1[i], reads_1_sub)
  ShortRead::writeFastq(sr2[i], reads_2_sub)
}
