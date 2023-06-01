
#' @export
#' @importFrom rlang is_scalar_integerish is_scalar_character
#' @importFrom dplyr mutate filter left_join select add_count inner_join
#' @importFrom magrittr "%>%"
#' @importFrom Biostrings DNAStringSet
#' @importFrom purrr map2 map_dbl
#' @importFrom tidyr chop
call_vars <- function(seq_tbl, marker_info) {
  proc_vars <- process_vars(seq_tbl, marker_info)

  vars_tbl <-
    proc_vars %>%
    mutate(var_data = map(var_data, "vars")) %>%
    unnest(var_data) %>%
    mutate(var_id = str_c(marker_id, var_id, sep = "_")) %>%
    select(marker_id, var_id, everything())

  return(vars_tbl)
}

#' @export
seq_mask_vars <- function(seq_tbl, marker_info, vars_mask) {
  check_table(vars_mask,
    col_types = list(
      marker_id = "character",
      aln_pos = "integer",
      type = "character"
    ),
    arg_name = "vars_mask"
  )

  proc_vars <- process_vars(seq_tbl, marker_info)

  masked_seqs <-
    proc_vars %>%
    mutate(var_data = map(var_data, "aln_tbl")) %>%
    unnest(var_data) %>%
    filter(aln_pos > 0) %>%
    mutate(type = case_when(
      ref == asv_base ~ "ident",
      ref == "-" ~ "INS",
      asv_base == "-" ~ "DEL",
      TRUE ~ "SNP"
    )) %>%
    left_join(mutate(vars_mask, mask = TRUE),
      by = c("marker_id", "aln_pos", "type")
    ) %>%
    mutate(mask = replace_na(mask, FALSE)) %>%
    group_by(marker_id, seq_id) %>%
    summarise(
      sequence = str_c(if_else(mask, ref, asv_base), collapse = "") %>% str_remove_all("-"),
      masked = any(mask),
      .groups = "drop"
    )

  seq_tbl_masked <-
    seq_tbl %>%
    select(sample_id, marker_id, seq_id, count) %>%
    inner_join(masked_seqs, c("marker_id", "seq_id")) %>%
    group_by(sample_id, marker_id, sequence) %>%
    summarise(
      count = sum(count),
      masked = any(masked),
      seq_id = digest::digest(sequence[1]),
      .groups = "drop"
    ) %>%
    select(sample_id, marker_id, seq_id, sequence, count, masked) %>%
    annotate_seq_tbl(
      marker_info = marker_info,
      min_marker_count = 0,
      min_asv_count = 0,
      min_asv_freq = 0,
      min_ident = 0,
      min_ident_z = -Inf,
      mark_chimeras = FALSE
    )

  return(seq_tbl_masked)
}

process_vars <- function(seq_tbl, marker_info) {

  # check args
  check_seq_table(seq_tbl, seq_id = TRUE)
  check_marker_info(marker_info)
  if (!all(seq_tbl$marker_id %in% marker_info$marker_id)) {
    abort("not all marker_ids in asv_tbl present in marker_info")
  }

  seq_tbl %>%
    select(marker_id, sequence, seq_id) %>%
    na.omit() %>%
    distinct() %>%
    chop(-marker_id) %>%
    mutate(dnass = map2(seq_id, sequence, function(aid, seq) {
      DNAStringSet(seq) %>% setNames(aid)
    })) %>%
    inner_join(select(marker_info, marker_id, marker_seq = seq), "marker_id") %>%
    mutate(var_data = map2(marker_seq, dnass, function(marker_seq, seqs) {
      call_variants(DNAStringSet(marker_seq), seqs)
    })) %>%
    select(marker_id, var_data)
}



#' @importFrom rlang is_scalar_integerish is_scalar_character
#' @importFrom dplyr mutate filter left_join select if_else case_when summarise group_by arrange_all anti_join
#' @importFrom magrittr "%>%"
#' @importFrom Biostrings DNAStringSet
#' @importFrom purrr map2 pmap_df
#' @importFrom stringr str_c str_remove_all
#' @importFrom tidyr pivot_longer replace_na complete
call_variants <- function(ref_seq, var_seqs) {
  # check args
  stopifnot(
    is(ref_seq, "DNAStringSet") && length(ref_seq) == 1,
    is(var_seqs, "DNAStringSet") && length(var_seqs) > 0
  )

  n <- length(var_seqs)

  aln_tbl <-
    DECIPHER::AlignSeqs(c(setNames(ref_seq, "ref"), var_seqs), verbose = FALSE) %>%
    as.matrix() %>%
    t() %>%
    as_tibble() %>%
    mutate(
      pos = cumsum(ref != "-"),
      aln_pos = seq_along(pos)
    ) %>%
    pivot_longer(c(-ref, -pos, -aln_pos),
      names_to = "seq_id",
      values_to = "asv_base"
    ) %>%
    (function(x) {
      tibble(ref = "N", pos = 0L, aln_pos = 0L, asv_base = "N", seq_id = unique(x$seq_id)) %>%
        bind_rows(x)
    }) %>%
    arrange(seq_id, aln_pos)


  vars <- if (length(unique(aln_tbl$seq_id)) == 1 & all(aln_tbl$ref == aln_tbl$asv_base)) {
    return(NULL)
    aln_tbl <- NULL
  } else {
    aln_tbl %>%
      split.data.frame(.$seq_id) %>%
      map_df(function(data) {
        data %>%
          filter(!(ref == "-" & asv_base == "-")) %>%
          mutate(
            ins_grp_start = replace_na(lead(ref) == "-" & ref != "-", FALSE),
            ins_grp_end = replace_na(ref != "-" & lag(ref) == "-", FALSE),
            ins_id = if_else(cumsum(ins_grp_start) > cumsum(ins_grp_end),
              cumsum(ins_grp_start),
              NA_integer_
            ),
            del_grp_start = replace_na(lead(asv_base) == "-" & asv_base != "-", FALSE),
            del_grp_end = replace_na(asv_base != "-" & lag(asv_base) == "-", FALSE),
            del_id = if_else(cumsum(del_grp_start) > cumsum(del_grp_end),
              cumsum(del_grp_start),
              NA_integer_
            ),
            is_snp = is.na(ins_id) & is.na(del_id) & ref != asv_base,
            snp_id = if_else(is_snp,
              cumsum(is_snp),
              NA_integer_
            )
          ) %>%
          select(-ins_grp_start, -ins_grp_end, -del_grp_start, -del_grp_end, -is_snp) %>%
          filter(!is.na(snp_id) | !is.na(ins_id) | !is.na(del_id)) %>%
          mutate(type = case_when(
            !is.na(snp_id) ~ "SNP",
            !is.na(ins_id) ~ "INS",
            !is.na(del_id) ~ "DEL"
          )) %>%
          group_by(ins_id, del_id, snp_id) %>%
          summarise(
            aln_pos = aln_pos[1],
            pos = pos[1],
            ref = str_c(ref, collapse = "") %>% str_remove_all("-"),
            asv_base = str_c(asv_base, collapse = "") %>% str_remove_all("-"),
            type = type[1],
            .groups = "drop"
          ) %>%
          mutate(seq_id = data$seq_id[1]) %>%
          select(seq_id, aln_pos, pos, ref, alt = asv_base, type)
      }) %>%
      chop(c(seq_id)) %>%
      arrange(aln_pos, pos, ref, alt) %>%
      chop(c(-aln_pos, -pos, -ref, -type)) %>%
      mutate(genotype = map(alt, seq_along)) %>%
      unnest(c(seq_id, genotype)) %>%
      unnest(seq_id) %>%
      group_by(aln_pos, pos, ref, alt, type) %>%
      complete(seq_id = names(var_seqs), fill = list(genotype = 0L)) %>%
      ungroup() %>%
      mutate(num_allele = lengths(alt) + 1L) %>%
      nest(gt_data = c(seq_id, genotype)) %>%
      mutate(var_id = str_c(pos, ref, type, sep = "_")) %>%
      select(var_id, everything()) %>%
      arrange(pos, var_id)
  }


  # find del sites to replace allele with '*'
  if (any(vars$type == "DEL")) {
    fixed <-
      vars %>%
      filter(type == "DEL") %>%
      mutate(end = pos + nchar(ref) - 1L) %>%
      unnest(gt_data) %>%
      filter(genotype != 0) %>%
      select(pos, end, seq_id) %>%
      chop(seq_id) %>%
      pmap_df(function(pos, end, seq_id) {
        filter(vars, pos >= !!pos, pos <= !!end, type != "DEL") %>%
          (function(x) {
            if (nrow(x) > 0) {
              x %>%
                mutate(
                  alt = map(alt, ~ c(., "")) %>% vctrs::as_list_of(),
                  num_allele = num_allele + 1L,
                  gt_data = map2(gt_data, num_allele, function(data, na) {
                    mutate(data, genotype = if_else(seq_id %in% !!seq_id, na - 1L, genotype))
                  })
                )
            }
          })
      })

    if (nrow(fixed) > 0) {
      vars <-
        anti_join(vars, fixed, by = "var_id") %>%
        bind_rows(fixed) %>%
        arrange(pos, var_id)
    }
  }

  return(list(vars = vars, aln_tbl = aln_tbl))
}

#' @export
#' @importFrom dplyr mutate filter inner_join group_by summarise select n
#' @importFrom magrittr "%>%"
#' @importFrom tidyr unnest
calc_var_stats <- function(seq_tbl, var_tbl) {
  # check args
  check_var_table(var_tbl)
  check_asv_table(seq_tbl, check_count = TRUE)

  smry_1 <-
    seq_tbl %>%
    select(sample_id, marker_id, seq_id, count) %>%
    inner_join(var_tbl %>% select(marker_id, var_id, gt_data) %>% unnest(gt_data),
      by = c("marker_id", "seq_id")
    ) %>%
    group_by(sample_id, marker_id, var_id, genotype) %>%
    summarise(count = sum(count), .groups = "drop_last") %>%
    mutate(sm_freq = count / sum(count))


  samplewise <-
    smry_1 %>%
    summarise(
      n_asv = n(),
      n_read = sum(count),
      MAF = sort(sm_freq, decreasing = TRUE)[2],
      exp_het = (n_read / (n_read - 1)) * (1 - sum(sm_freq**2)),
      .groups = "drop"
    ) %>%
    mutate(MAF = replace_na(MAF, 0))

  varwise <-
    smry_1 %>%
    group_by(marker_id, var_id, genotype) %>%
    summarise(
      count = sum(sm_freq),
      .groups = "drop_last"
    ) %>%
    mutate(freq = count / sum(count)) %>%
    summarise(
      n_sample = sum(count),
      n_allele = n(),
      MAF = sort(freq, decreasing = TRUE)[2],
      exp_het = (n_sample / (n_sample - 1)) * (1 - sum(freq**2)),
      .groups = "drop"
    ) %>%
    mutate(MAF = replace_na(MAF, 0)) %>%
    left_join(samplewise %>%
      group_by(marker_id, var_id) %>%
      summarise(
        sample_med_MAF = median(MAF),
        sample_med_exp_het = median(exp_het),
        .groups = "drop"
      ),
    by = c("marker_id", "var_id")
    )

  return(list(
    varwise = inner_join(select(var_tbl, var_id, pos, aln_pos, type),
      varwise,
      by = "var_id"
    ),
    samplewise = inner_join(select(var_tbl, var_id, pos, aln_pos, type),
      samplewise,
      by = "var_id"
    )
  ))
}

remove_variants <- function(seq_tbl, var_tbl, marker_info) {
  stopifnot(
    is.data.frame(seq_tbl),
    is.data.frame(var_tbl),
    is.data.frame(marker_info),
    is_scalar_character(output_dir),
    is_bool(compress)
  )

  check_seq_table(seq_tbl)
  check_var_table(var_tbl)
  check_marker_info(marker_info)
}


#' @export
#' @importFrom dplyr mutate select if_else first slice arrange pull
#' @importFrom rlang is_scalar_character abort is_bool
#' @importFrom magrittr "%>%"
#' @importFrom Biostrings DNAStringSet
#' @importFrom purrr map2_df pmap
export_vcf <- function(var_tbl, marker_info, output_dir, use_absolute_paths = TRUE, compress = TRUE) {
  stopifnot(
    is.data.frame(var_tbl),
    is.data.frame(marker_info),
    is_scalar_character(output_dir),
    is_bool(compress)
  )

  check_var_table(var_tbl)
  check_marker_info(marker_info, locus = TRUE)

  if (!file.exists(output_dir)) {
    dir.create(output_dir)
  }
  output_dir <- normalizePath(output_dir)

  if (!use_absolute_paths) {
    output_sub_dir <- gsub(pattern = normalizePath(getwd()), replacement = ".", x = output_sub_dir, fixed = TRUE)
  }

  if (!all(var_tbl$marker_id %in% marker_info$marker_id)) {
    abort("not all marker_ids in asv_tbl present in marker_info")
  }

  res <-
    var_tbl %>%
    split.data.frame(.$marker_id) %>%
    map2_df(names(.), ., function(marker_id, data) {
      vcf_fn <- file.path(output_dir, str_c(marker_id, ".vcf", sep = ""))

      locus <- marker_info %>% dplyr::slice(first(which(marker_id == !!marker_id)))

      rowRanges <- GenomicRanges::GRanges(
        locus$chrom,
        IRanges::IRanges(start = locus$start + data$pos - 1L, width = 1L)
      )

      seq_ids <- data$gt_data %>%
        map("seq_id") %>%
        unlist() %>%
        sort() %>%
        unique()

      colData <- S4Vectors::DataFrame(Samples = seq_along(seq_ids), row.names = seq_ids)

      header <-
        VariantAnnotation::VCFHeader(
          samples = seq_ids,
          header = IRanges::DataFrameList(
            FORMAT = S4Vectors::DataFrame(
              Number = "1",
              Type = "String",
              Description = "Genotype",
              row.names = "GT"
            ),
            fileformat = S4Vectors::DataFrame(
              Value = "VCFv4.2",
              row.names = "fileformat"
            )
          )
        )

      geno <-
        data %>%
        pmap(function(var_id, gt_data, ...) {
          gt_data %>%
            arrange(match(seq_id, seq_ids)) %>%
            pull(genotype) %>%
            (function(x) {
              matrix(
                data = as.character(x), nrow = 1,
                dimnames = list(row = var_id, col = seq_ids)
              )
            })
        }) %>%
        do.call("rbind", .) %>%
        S4Vectors::SimpleList(GT = .)

      fixed <-
        with(
          data,
          S4Vectors::DataFrame(
            REF = DNAStringSet(ref),
            ALT = map(alt, DNAStringSet) %>% Biostrings::DNAStringSetList(),
            QUAL = 1,
            FILTER = "."
          )
        )

      vcf <- VariantAnnotation::VCF(
        rowRanges = rowRanges,
        colData = colData,
        exptData = list(header = header),
        fixed = fixed,
        geno = geno
      )
      VariantAnnotation::writeVcf(vcf, vcf_fn, index = compress)
      tibble(marker_id = marker_id, vcf_filename = str_c(vcf_fn, if_else(compress, ".bgz", ""), sep = ""))
    })

  return(res)
}
