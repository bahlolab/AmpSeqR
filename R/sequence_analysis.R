#' Amplicon haplotype inference
#'
#' Haplotype inference for amplicons.
#'
#' @param seq_tbl (Required). The haplotype sequence table includes sample_id, marker_id, sequence (haplotype sequence), count (haplotype read counts), status (haplotype status), sample, info.
#' @param marker_info (Required). The target amplicon file. The file should include marker_id, primer_fwd, primer_rev, seq (reference sequence), chrom (chromosome), start (reference sequence start position), end (reference sequence end position).
#' @param output_dir (Required). The path to the output file.
#' @param threads (Optional). Default 8. The number of processors to use, or NULL to automatically detect and use all available processors.
#' @param min_marker_count (Optional). Default 100. The minimum number of reads for each amplicon marker after calling haplotype.
#' @param min_asv_count (Optional). Default 5. The minimum number of reads for each haplotype.
#' @param min_asv_freq (Optional). Default 0.001. The minimum haplotype frequency.
#' @param min_ident (Optional). Default 0.75. The minimum sequence similarity.
#' @param min_ident_z (Optional). Default -3. The minimum standardized sequence similarity.
#' @param mark_chimeras (Optional). Default TRUE. If true, check chimeric reads.
#' @param max_breakpoints (Optional). Default 3. The maximum number of breakpoints chimera reads.
#' @param min_parent_ratio min_parent_ratio (Optional). Default 1.5. The minimum ratio that the sequences greater than this value are more abundant than a sequence can be its "parents".
#'
#' @return
#'
#' seq_ann_tbl.rds table: the annotate_seq_tbl table in RDS format includes sample_id, marker_id, sequence (haplotype sequence), count (haplotype read counts), status (haplotype status after annotation), sample, info, ident (sequence similarity), ident_z (standardized sequence similarity).
#'
#' @export
annotate_seq_tbl <- function(seq_tbl,
                             marker_info,
                             output_dir = NULL,
                             threads = 8L,
                             min_marker_count = 100L,
                             min_asv_count = 5L,
                             min_asv_freq = 1e-3,
                             min_ident = 0.75,
                             min_ident_z = -3,
                             mark_chimeras = TRUE,
                             max_breakpoints = 3L,
                             min_parent_ratio = 1.5) {
  # check args
  stopifnot(
    is.data.frame(seq_tbl),
    is.data.frame(marker_info),
    is_scalar_integerish(min_marker_count) && min_marker_count >= 0L,
    is_scalar_integerish(min_asv_count) && min_asv_count >= 0L,
    is_scalar_double(min_asv_freq) && min_asv_freq >= 0,
    is_scalar_integerish(threads) && threads >= 1L,
    is_scalar_double(min_ident) || is.null(min_ident),
    is_scalar_double(min_ident_z) || is.null(min_ident_z)
  )

  check_seq_table(seq_tbl)
  check_marker_info(marker_info)

  seq_tbl_clean <-
    seq_tbl %>%
    # annotate low_sample_count
    group_by(sample_id, marker_id) %>%
    mutate(sample_marker_count = sum(count)) %>%
    ungroup() %>%
    add_status(., if_else(.$sample_marker_count >= min_marker_count, "pass", "low_sample_count")) %>%
    select(-sample_marker_count) %>%
    # calculate seq ident, add annotations
    calc_seq_ident(
      marker_info = marker_info,
      min_ident = min_ident,
      min_ident_z = min_ident_z
    ) %>%
    # annotate low_asv_count and low_asv_freq
    group_by(sample_id, marker_id) %>%
    mutate(sample_asv_freq = count / sum(count)) %>%
    ungroup() %>%
    add_status(., if_else(.$count >= min_asv_count, "pass", "low_asv_count")) %>%
    add_status(., if_else(.$sample_asv_freq >= min_asv_freq, "pass", "low_asv_freq")) %>%
    select(-sample_asv_freq)

  if (mark_chimeras) {
    seq_tbl_clean <-
      mark_chimeras(seq_tbl_clean,
        threads = threads,
        max_breakpoints = max_breakpoints,
        min_parent_ratio = min_parent_ratio
      )
  }

  if (is.null(output_dir) == FALSE) {
    write_rds(seq_tbl_clean, paste(output_dir, "/seq_ann_tbl.rds", sep = ""))
  }

  return(seq_tbl_clean)
}

#' @export
#' @importFrom dplyr select mutate filter
#' @importFrom rlang is_scalar_integerish is_scalar_character is_scalar_double
#' @importFrom magrittr "%>%"
#' @importFrom Biostrings DNAStringSet
#' @importFrom purrr map_lgl
#' @importFrom tidyr chop unnest
calc_seq_ident <- function(seq_tbl, marker_info,
                           threads = 1L,
                           min_ident = 0.75,
                           min_ident_z = -3) {

  # check args
  stopifnot(
    is.data.frame(seq_tbl),
    is.data.frame(marker_info),
    is_scalar_integerish(threads) && threads >= 1L,
    is_scalar_double(min_ident) || is.null(min_ident),
    is_scalar_double(min_ident_z) || is.null(min_ident_z)
  )

  check_seq_table(seq_tbl)
  check_marker_info(marker_info)

  # remove any missing reads
  ident_row <-
    seq_tbl %>%
    select(marker_id, sequence) %>%
    mutate(row = seq_along(sequence)) %>%
    filter(!is.na(sequence)) %>%
    chop(row) %>%
    chop(c(sequence, row)) %>%
    left_join(select(marker_info, marker_id, marker_seq = seq), "marker_id") %>%
    mutate(ident = map2(marker_seq, sequence, function(marker_seq, sequence) {
      DNAStringSet(c(marker_seq, sequence)) %>%
        DECIPHER::AlignSeqs(verbose = FALSE, processors = threads) %>%
        DECIPHER::DistanceMatrix(verbose = FALSE, processors = threads, includeTerminalGaps = TRUE) %>%
        (function(x) 1 - x[, 1][-1])
    })) %>%
    select(row, ident) %>%
    unnest(c(row, ident)) %>%
    unnest(row)

  ident_tbl <-
    mutate(seq_tbl, ident = NA_real_) %>%
    mutate(ident = replace(ident, ident_row$row, ident_row$ident)) %>%
    group_by(marker_id) %>%
    mutate(ident_z = c(scale(replace(ident, ident < min_ident, NA_real_)))) %>%
    ungroup()

  if (!is.null(min_ident)) {
    ident_tbl <-
      ident_tbl %>%
      add_status(if_else(ident_tbl$ident < min_ident, "low_ident", "pass"))
  }

  if (!is.null(min_ident_z)) {
    ident_tbl <-
      ident_tbl %>%
      add_status(if_else(ident_tbl$ident_z < min_ident_z, "low_ident_z", "pass"))
  }


  return(ident_tbl)
}
#' @export
#' @importFrom dplyr select mutate filter slice arrange pull n
#' @importFrom rlang is_scalar_integerish is_scalar_character is_scalar_double abort
#' @importFrom magrittr "%>%"
#' @importFrom Biostrings DNAStringSet
#' @importFrom purrr map map_lgl
#' @importFrom tidyr nest unnest chop unchop
mark_chimeras <- function(seq_tbl,
                          threads = 1L,
                          max_breakpoints = 3L,
                          min_parent_ratio = 1.5,
                          pass_only = TRUE) {

  # check args
  stopifnot(
    is.data.frame(seq_tbl),
    is_scalar_integerish(threads) && threads >= 1L,
    is_scalar_integerish(max_breakpoints) && max_breakpoints >= 1L,
    is_scalar_double(min_parent_ratio) && min_parent_ratio > 1,
    is_bool(pass_only)
  )

  check_seq_table(seq_tbl)

  cluster <- `if`(
    threads > 1,
    future::makeClusterPSOCK(workers = threads),
    NULL
  )

  chimeric <-
    seq_tbl %>%
    mutate(row = seq_len(n())) %>%
    arrange(sample_id, marker_id, desc(count)) %>%
    group_by(sample_id, marker_id) %>%
    filter(!pass_only | (status == "pass")) %>%
    mutate(cand = map_lgl(seq_along(count), ~ sum(count >= min_parent_ratio * count[.]) > 1)) %>%
    filter(any(cand)) %>%
    (function(x) {
      `if`(
        nrow(x) > 0,
        nest(x) %>%
          ungroup() %>%
          (function(y) {
            `if`(
              threads > 1,
              mutate(y, group = seq_len(n()) %% threads) %>%
                chop(-group) %>%
                mutate(data = map(data, function(data) {
                  future::cluster(
                    {
                      purrr::map(data, AmpSeqR:::mark_chimeras_mapper,
                        max_breakpoints = max_breakpoints,
                        min_parent_ratio = min_parent_ratio
                      )
                    },
                    workers = cluster
                  )
                }) %>% map(future::value)) %>%
                unchop(-group) %>%
                select(-group),
              mutate(y, data = map(data, mark_chimeras_mapper,
                max_breakpoints = max_breakpoints,
                min_parent_ratio = min_parent_ratio
              ))
            )
          }) %>%
          unnest(data) %>%
          filter(is_chimeric) %>%
          pull(row),
        integer(0)
      )
    })

  if (threads > 1) parallel::stopCluster(cluster)

  chim_tbl <-
    seq_tbl %>%
    add_status(if_else(seq_len(nrow(seq_tbl)) %in% chimeric, "chimera", "pass"))

  return(chim_tbl)
}
#' @export
mark_chimeras_mapper <- function(data, max_breakpoints, min_parent_ratio) {
  data$is_chimeric <- data$cand
  for (i in which(data$cand)) {
    subject <- DNAStringSet(data$sequence[i])
    pool <-
      dplyr::slice(data, seq_len(i - 1)) %>%
      filter(count >= min_parent_ratio * data$count[i]) %>%
      with(DNAStringSet(sequence))
    data$is_chimeric[i] <- `if`(length(pool) < 2, FALSE, is_chimeric(subject, pool, max_breakpoints))
  }
  select(data, -cand)
}
#' @export
# return list of possible parents
#' @importFrom rlang is_scalar_character is_character
#' @importFrom dplyr lag lead n
#' @importFrom purrr map_chr map_lgl map_int
#' @importFrom magrittr set_names "%>%"
#' @importFrom tidyr expand_grid
#' @importFrom tibble tibble
is_chimeric <- function(subject, pool, max_breakpoints = 2L) {
  # check args
  stopifnot(
    is(subject, "DNAStringSet") && length(subject) == 1,
    is(pool, "DNAStringSet") && length(pool) > 1,
    is_scalar_integerish(max_breakpoints) && max_breakpoints >= 1L
  )

  if (any(subject == pool)) {
    rlang::warn("putative chimeric subject sequence contained in sequence pool")
    return(TRUE)
  }

  if (length(pool) > 2) {
    if (is_chimeric(subject, pool[seq_len(length(pool) - 1)])) {
      return(TRUE)
    }
  }

  n_par <- length(pool)

  aln_mat <-
    DECIPHER::AlignSeqs(c(subject, pool), verbose = FALSE) %>%
    as.matrix() %>%
    (function(x) {
      x[, !map_lgl(seq_len(ncol(x)), function(i) all(x[1, i] == x[-1, i])), drop = FALSE]
    })

  is_unique <- map_lgl(seq_len(ncol(aln_mat)), function(i) !any(aln_mat[1, i] == aln_mat[-1, i]))

  if (any(is_unique)) {
    return(FALSE)
  }

  # internal function used by is_chimera - note that aln_mat must be defined in calling env
  # depth first search for possible chimeras
  is_chimera_recursive <- function(pos, path) {
    pos <- pos + 1L
    if (pos > ncol(aln_mat)) {
      return(TRUE)
    }
    par <- which(aln_mat[1, pos] == aln_mat[-1, pos])
    breaks <- map_int(par, function(p) {
      c(path, p) %>%
        {
          sum(. != lead(.), na.rm = TRUE)
        }
    })
    for (p in par[breaks <= max_breakpoints]) {
      if (is_chimera_recursive(pos, c(path, p))) {
        return(TRUE)
      }
    }
    return(FALSE)
  }

  return(is_chimera_recursive(pos = 0L, path = integer()))
}
#' @export
#' @importFrom dplyr select mutate filter
#' @importFrom rlang is_scalar_integerish is_scalar_character abort
#' @importFrom magrittr "%>%"
#' @importFrom Biostrings DNAStringSet
#' @importFrom purrr map_df
#' @importFrom tidyr everything
# assign_asv_id <- function(seq_tbl, threads = 1L, pass_only = TRUE) {
#
#   # check args
#   stopifnot(is.data.frame(seq_tbl),
#             is_scalar_integerish(threads) && threads >= 1L,
#             is_bool(pass_only))
#
#   check_seq_table(seq_tbl)
#
#   asv_id_tbl <-
#     seq_tbl %>%
#     filter(!pass_only | replace_na(str_detect(status, 'pass'))) %>%
#     select(marker_id, sequence) %>%
#     arrange(marker_id, sequence) %>%
#     distinct() %>%
#     add_count(marker_id) %>%
#     mutate(digits = floor(log10(max(n))) + 1) %>%
#     split.data.frame(.$marker_id) %>%
#     map_df(function(x) {
#       mutate(x, id =
#                DNAStringSet(sequence) %>%
#                DECIPHER::AlignSeqs(verbose = FALSE, processors = threads) %>%
#                DECIPHER::DistanceMatrix(verbose = FALSE, processors = threads, includeTerminalGaps = TRUE) %>%
#                as.dist() %>% hclust() %>%  { .$order }) %>%
#         arrange(id) %>%
#         mutate(asv_id = str_c(marker_id, format_int(id, digits[1]), sep = '-')) %>%
#         select(marker_id, asv_id, sequence)
#     }) %>%
#     right_join(seq_tbl, by = c("marker_id", "sequence")) %>%
#     select(sample_id, marker_id, asv_id, sequence, everything()) %>%
#     arrange(sample_id, marker_id, count, sequence)
#
#   return(asv_id_tbl)
# }
