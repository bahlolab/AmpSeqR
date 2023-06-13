#' Clean homopolymers indels
#'
#' Detect indels exist in homopolymers and change the indel to be the same as in the reference genome.
#'
#'
#'
#' @param seq_ann_tbl (Required). The haplotype sequence table includes sample_id, marker_id, sequence (haplotype sequence), count (haplotype read counts).
#' @param marker_info (Required). The target amplicon table. The file should include marker_id, primer_fwd, primer_rev, seq (reference sequence).
#' @param min_homo_rep (Optional). Default 3 (numeric vector). The minimum length of the homopolymer repeats.
#'
#' @return
#'
#' The haplotype sequence with homopolymers indels are changed to be the same as in the reference genome.
#'
#' @export
#'
#' @author Jiru Han, Jacob E. Munro, Melanie Bahlo
#'
#' @examples
#' seq_ann_tbl <- data.frame(
#'   sample_id = c("S1", "S2", "S3", "S4"),
#'   marker_id = c("TRAP"),
#'   sequence = c(
#'     "TGCTAGTGTTTTTTCAAACAATGCAAGAGAAATTATTAGATTA",
#'     "TGCTAGTGTTTTCAAACAATGCAAGAGAAATTATTAGATTA",
#'     "TGCTAGTGTTTTTCAAACAATGCAAGAGAAATTATTAGATTA",
#'     "TGCTAGTGTTTTTTTTTCAAACAATGCAAGAGAAATTATTAGATTA"
#'   ),
#'   count = as.integer(10000)
#' )
#'
#' marker_info <- data.frame(
#'   marker_id = c("TRAP"),
#'   primer_fwd = c("CTTCTACGTCTTACAAAG"),
#'   primer_rev = c("CGTGATGATCCTTGCAGA"),
#'   seq = c("TGCTAGTGTTTTTTCAAACAATGCAAGAGAAATTATTAGATTA")
#' )
#'
#' clean_homopolymers(seq_ann_tbl, marker_info)
#'
clean_homopolymers <- function(seq_ann_tbl,
                               marker_info,
                               min_homo_rep = 3) {
  stopifnot(
    is.data.frame(seq_ann_tbl),
    is.data.frame(marker_info)
  )
  
  check_seq_table(seq_ann_tbl)
  check_marker_info(marker_info)
  
  homo_seq <-
    seq_ann_tbl %>%
    select(marker_id, sequence) %>%
    mutate(row = seq_along(sequence)) %>%
    filter(!is.na(sequence)) %>%
    chop(row) %>%
    chop(c(sequence, row)) %>%
    left_join(select(marker_info, marker_id, marker_seq = seq), "marker_id") %>%
    mutate(homo_seq = map2(marker_seq, sequence, function(marker_seq, sequence) {
      # Reference sequence
      ref <- DNAStringSet(marker_seq)
      ref_tbl <- ref %>%
        as.matrix() %>%
        t() %>%
        as_tibble()
      
      # Check reference homopolymer
      ref_homopolymer <- data.frame(base = rle(ref_tbl$V1)$values, width = rle(ref_tbl$V1)$lengths)
      ref_homopolymer <- ref_homopolymer %>%
        mutate(end = cumsum(width)) %>%
        mutate(start = end - width + 1) %>%
        select(start, end, width, base)
      
      # Align sequence
      check_homo_seq_tbl <- DECIPHER::AlignSeqs(c(setNames(ref, "ref"), DNAStringSet(sequence)), verbose = FALSE) %>%
        as.matrix() %>%
        t() %>%
        as_tibble() %>%
        mutate(
          pos = cumsum(ref != "-"),
          aln_pos = seq_along(pos)
        )
      
      # all homopolymer regions (based on reference sequence, >= n repeated bases)
      ref_homopolymer_select <- ref_homopolymer %>%
        filter(width > min_homo_rep)
      
      for (i in 1:nrow(ref_homopolymer_select)) {
        for (j in 2:(ncol(check_homo_seq_tbl) - 2)) {
          check_homo_seq_tbl[(check_homo_seq_tbl$pos %in% (ref_homopolymer_select[i, 1]:ref_homopolymer_select[i, 2])), j] <- check_homo_seq_tbl[(check_homo_seq_tbl$pos %in% (ref_homopolymer_select[i, 1]:ref_homopolymer_select[i, 2])), 1]
        }
      }
      
      
      # Insertion
      for (i in 1:nrow(ref_homopolymer_select)) {
        for (j in 2:(ncol(check_homo_seq_tbl) - 2)) {
          check_homo_seq_tbl[(check_homo_seq_tbl$ref == "-" & check_homo_seq_tbl$pos %in% ((ref_homopolymer_select[i, 1] - 1):(ref_homopolymer_select[i, 2] + 1))), j] <-
            check_homo_seq_tbl[(check_homo_seq_tbl$ref == "-" & check_homo_seq_tbl$pos %in% ((ref_homopolymer_select[i, 1] - 1):(ref_homopolymer_select[i, 2] + 1))), 1]
        }
      }
      
      change_homo_seq <- lapply(2:(ncol(check_homo_seq_tbl) - 2), function(i) {
        str_c(unlist(check_homo_seq_tbl[, i]), collapse = "") %>% str_remove_all("-")
      })
      
      change_homo_seq <- do.call(rbind, change_homo_seq)
    })) %>%
    select(row, homo_seq) %>%
    unnest(c(row, homo_seq)) %>%
    unnest(row)
  
  # Change the homopolymer size same as reference sequence
  seq_ann_tbl <- suppressMessages(seq_ann_tbl %>%
                                    mutate(sequence = replace(sequence, homo_seq$row, homo_seq$homo_seq)) %>%
                                    group_by(sample_id, marker_id, sequence) %>%
                                    mutate(count = sum(count)) %>%
                                    slice(1) %>%
                                    ungroup())
  
  seq_ann_tbl
}



#' Clean terminal indels
#'
#' Change terminal indels to be the same as in the reference genome.
#'
#' @param seq_ann_tbl (Required). The haplotype sequence table includes sample_id, marker_id, sequence (haplotype sequence), count (haplotype read counts).
#' @param marker_info (Required). The target amplicon table. The file should include marker_id, primer_fwd, primer_rev, seq (reference sequence).
#' @param terminal_region_len (Optional). Default 1 (numeric vector between 1 and 3). The terminal indel position.
#'
#' @return
#'
#' The haplotype sequence with terminal indels are changed to be the same as in the reference genome.
#'
#' @export
#'
#' @author Jiru Han, Jacob E. Munro, Melanie Bahlo
#'
#' @examples
#' seq_ann_tbl <- data.frame(
#'   sample_id = c("S1", "S2", "S3", "S4"),
#'   marker_id = c("TRAP"),
#'   sequence = c(
#'     "TGCTAGTGTTTTTCAAACAATGCAAGAGAAATTATTAGATT",
#'     "GCTAGTGTTTTTCAAACAATGCAAGAGAAATTATTAGATTA",
#'     "GCTAGTGTTTTTCAAACAATGCAAGAGAAATTATTAGATT",
#'     "TGCTAGTGTTTTTTCAAACAATGCAAGAGAAATTATTAGATTA"
#'   ),
#'   count = as.integer(10000)
#' )
#'
#' marker_info <- data.frame(
#'   marker_id = c("TRAP"),
#'   primer_fwd = c("CTTCTACGTCTTACAAAG"),
#'   primer_rev = c("CGTGATGATCCTTGCAGA"),
#'   seq = c("TGCTAGTGTTTTTTCAAACAATGCAAGAGAAATTATTAGATTA")
#' )
#'
#' clean_terminal_indels(seq_ann_tbl, marker_info)
#'
clean_terminal_indels <- function(seq_ann_tbl,
                                  marker_info,
                                  terminal_region_len = 1) {
  stopifnot(
    is.data.frame(seq_ann_tbl),
    is.data.frame(marker_info)
  )
  
  check_seq_table(seq_ann_tbl)
  check_marker_info(marker_info)
  
  terminal_seq <-
    seq_ann_tbl %>%
    select(marker_id, sequence) %>%
    mutate(row = seq_along(sequence)) %>%
    filter(!is.na(sequence)) %>%
    chop(row) %>%
    chop(c(sequence, row)) %>%
    left_join(select(marker_info, marker_id, marker_seq = seq), "marker_id") %>%
    mutate(terminal_seq = map2(marker_seq, sequence, function(marker_seq, sequence) {
      # Reference sequence
      ref <- DNAStringSet(marker_seq)
      
      if (terminal_region_len == 1) {
        ref_fix_pos <- c(0, 1, width(ref), (width(ref) + 1))
      } else if (terminal_region_len == 2) {
        ref_fix_pos <- c(0, 1, 2, (width(ref) - 1), width(ref), (width(ref) + 1))
      } else {
        ref_fix_pos <- c(0, 1, 2, 3, (width(ref) - 2), (width(ref) - 1), width(ref), (width(ref) + 1))
      }
      
      
      # Align sequence
      check_terminal_seq_tbl <- DECIPHER::AlignSeqs(c(setNames(ref, "ref"), DNAStringSet(sequence)), verbose = FALSE) %>%
        as.matrix() %>%
        t() %>%
        as_tibble() %>%
        mutate(
          pos = cumsum(ref != "-"),
          aln_pos = seq_along(pos)
        )
      
      # terminal regions (based on reference sequence)
      
      for (i in 2:(ncol(check_terminal_seq_tbl) - 2)) {
        check_terminal_seq_tbl[((check_terminal_seq_tbl$pos %in% ref_fix_pos) & check_terminal_seq_tbl$ref == "-"), i] <-
          check_terminal_seq_tbl[((check_terminal_seq_tbl$pos %in% ref_fix_pos) & check_terminal_seq_tbl$ref == "-"), 1]
      }
      
      # deletion
      for (i in 2:(ncol(check_terminal_seq_tbl) - 2)) {
        check_terminal_seq_tbl[((check_terminal_seq_tbl$pos %in% ref_fix_pos) & check_terminal_seq_tbl[, i] == "-"), i] <-
          check_terminal_seq_tbl[((check_terminal_seq_tbl$pos %in% ref_fix_pos) & check_terminal_seq_tbl[, i] == "-"), 1]
      }
      
      change_terminal_seq <- lapply(2:(ncol(check_terminal_seq_tbl) - 2), function(i) {
        str_c(unlist(check_terminal_seq_tbl[, i]), collapse = "") %>% str_remove_all("-")
      })
      
      change_terminal_seq <- do.call(rbind, change_terminal_seq)
    })) %>%
    select(row, terminal_seq) %>%
    unnest(c(row, terminal_seq)) %>%
    unnest(row)
  
  # Change the terminal regions indels same as reference sequence
  seq_ann_tbl <- suppressMessages(seq_ann_tbl %>%
                                    mutate(sequence = replace(sequence, terminal_seq$row, terminal_seq$terminal_seq)) %>%
                                    group_by(sample_id, marker_id, sequence) %>%
                                    mutate(count = sum(count)) %>%
                                    slice(1) %>%
                                    ungroup())
  
  seq_ann_tbl
}
