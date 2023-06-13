#' @export
#' @importFrom dplyr mutate filter left_join select summarise group_by rename ungroup
#' @importFrom magrittr "%>%"
#' @importFrom plyr join_all
#' @importFrom DT datatable
read_counts_tab <- function(demultiplexed, flt_reads, sub_reads, seq_ann_tbl, seq_flt_tbl) {
  dem <- demultiplexed %>%
    filter(!(is.na(sample_id))) %>%
    filter(!(is.na(marker_id))) %>%
    select(sample_id, sample, info, marker_id, n) %>%
    dplyr::rename(demultiplexed = n)
  
  flt <- flt_reads %>%
    select(sample_id, marker_id, n_out) %>%
    dplyr::rename(flt_reads = n_out)
  
  
  downsam <- sub_reads %>%
    select(sample_id, marker_id, n) %>%
    dplyr::rename(sub_reads = n)
  
  seq1 <- suppressMessages(seq_ann_tbl %>%
                             select(sample_id, marker_id, count, status) %>%
                             mutate(newstatus = ifelse(status == "pass", "pass", "noise")) %>%
                             group_by(sample_id, marker_id, newstatus) %>%
                             summarise(count = sum(count)) %>%
                             ungroup() %>%
                             filter(newstatus == "pass") %>%
                             select(sample_id, marker_id, count) %>%
                             dplyr::rename(seq_ann_tbl = count))
  
  seq2 <- suppressMessages(seq_flt_tbl %>%
                             select(sample_id, marker_id, count) %>%
                             group_by(sample_id, marker_id) %>%
                             summarise(count = sum(count)) %>%
                             ungroup() %>%
                             dplyr::rename(seq_flt_tbl = count))
  
  all_read <- plyr::join_all(list(dem, flt, downsam, seq1, seq2), by = c("sample_id", "marker_id"), type = "left")
  
  all_read %>%
    DT::datatable(
      extensions = "Buttons",
      options = list(
        dom = "Blfrtip",
        buttons = c("csv", "excel"),
        lengthMenu = list(
          c(10, 25, 50, -1),
          c(10, 25, 50, "All")
        )
      )
    )
}
#' @importFrom dplyr mutate filter left_join select summarise group_by rename ungroup
#' @importFrom magrittr "%>%"
#' @importFrom plyr join_all
#' @importFrom plotly ggplotly
read_counts_plot <- function(demultiplexed, flt_reads, sub_reads, seq_ann_tbl, seq_flt_tbl) {
  dem <- demultiplexed %>%
    filter(!(is.na(sample_id))) %>%
    filter(!(is.na(marker_id))) %>%
    select(sample_id, sample, info, marker_id, n) %>%
    dplyr::rename(demultiplexed = n)
  
  flt <- flt_reads %>%
    select(sample_id, marker_id, n_out) %>%
    dplyr::rename(flt_reads = n_out)
  
  
  downsam <- sub_reads %>%
    select(sample_id, marker_id, n) %>%
    dplyr::rename(sub_reads = n)
  
  seq1 <- suppressMessages(seq_ann_tbl %>%
                             select(sample_id, marker_id, count, status) %>%
                             mutate(newstatus = ifelse(status == "pass", "pass", "noise")) %>%
                             group_by(sample_id, marker_id, newstatus) %>%
                             summarise(count = sum(count)) %>%
                             ungroup() %>%
                             filter(newstatus == "pass") %>%
                             select(sample_id, marker_id, count) %>%
                             dplyr::rename(seq_ann_tbl = count))
  
  
  seq2 <- suppressMessages(seq_flt_tbl %>%
                             select(sample_id, marker_id, count) %>%
                             group_by(sample_id, marker_id) %>%
                             summarise(count = sum(count)) %>%
                             ungroup() %>%
                             dplyr::rename(seq_flt_tbl = count))
  
  all_read <- plyr::join_all(list(dem, flt, downsam, seq1, seq2), by = c("sample_id", "marker_id"), type = "left")
  
  all_read_plot <- all_read %>%
    mutate(final_prop_seqhap_retained = seq_flt_tbl / demultiplexed)
  
  p <- all_read_plot %>%
    ggplot(aes(sample, final_prop_seqhap_retained, fill = marker_id)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~marker_id, nrow = 1, scales = "free_x") +
    ggtitle("Final proportion sequence retained") +
    theme_light() +
    theme(
      axis.text.x = element_text(angle = 90)
    ) +
    ylim(c(0, 1))
  
  ggplotly(p)
}
#' @importFrom dplyr select mutate distinct group_by ungroup left_join filter left_join
#' @importFrom tidyr complete
#' @importFrom viridisLite cividis
#' @importFrom stats setNames
#' @importFrom plotly ggplotly layout
#' @importFrom magrittr "%>%"
miss_visual <- function(seq_flt_tbl, sample_manifest, marker_info) {
  missigness <-
    seq_flt_tbl %>%
    select(sample_id, marker_id, sample) %>%
    distinct() %>%
    mutate(missing = FALSE) %>%
    group_by(sample_id) %>%
    complete(marker_id = marker_info$marker_id, fill = list(missing = TRUE)) %>%
    group_by(marker_id) %>%
    complete(sample_id = sample_manifest$sample_id, fill = list(missing = TRUE)) %>%
    ungroup() %>%
    select(-sample) %>%
    left_join(
      (sample_manifest %>%
         select(sample_id, sample)),
      by = "sample_id"
    )
  
  percent_missing <-
    missigness$missing %>%
    (function(x) 100 * sum(x) / length(x))
  
  p <- missigness %>%
    ggplot(aes(sample, marker_id, fill = missing)) +
    geom_tile(col = "black") +
    scale_fill_manual(values = setNames(viridisLite::cividis(2), c("TRUE", "FALSE"))) +
    theme_light() +
    theme(
      axis.text.x = element_text(angle = 90)
    )
  
  ggplotly(p) %>%
    layout(title = list(x = 0.01, title_x = 0, text = paste0(
      str_c("Missingness (", format(percent_missing, digits = 3), "%)"),
      "<br>",
      str_c("Total sample: ", length(unique(sample_manifest$sample_id))),
      str_c("Total marker: ", length(unique(marker_info$marker_id)))
    )))
}
#' @importFrom dplyr select mutate distinct group_by ungroup left_join filter left_join summarise
#' @importFrom tidyr complete
#' @importFrom flextable flextable theme_vanilla
#' @importFrom magrittr "%>%"
marker_miss <- function(seq_flt_tbl, sample_manifest, marker_info) {
  missigness <-
    seq_flt_tbl %>%
    select(sample_id, marker_id, sample) %>%
    distinct() %>%
    mutate(missing = FALSE) %>%
    group_by(sample_id) %>%
    complete(marker_id = marker_info$marker_id, fill = list(missing = TRUE)) %>%
    group_by(marker_id) %>%
    complete(sample_id = sample_manifest$sample_id, fill = list(missing = TRUE)) %>%
    ungroup() %>%
    select(-sample) %>%
    left_join(
      (sample_manifest %>%
         select(sample_id, sample)),
      by = "sample_id"
    )
  
  suppressMessages(missigness %>%
                     group_by(marker_id, missing) %>%
                     summarise(count = n()) %>%
                     ungroup() %>%
                     group_by(marker_id) %>%
                     mutate(Missingness = round((100 - count / sum(count) * 100), digits = 2)) %>%
                     ungroup() %>%
                     complete(marker_id, missing, fill = list(count = 0, Missingness = 100)) %>%
                     filter(missing == FALSE) %>%
                     select(marker_id, Missingness) %>%
                     rename("Missingness(%)" = Missingness) %>%
                     flextable() %>%
                     theme_vanilla())
}
#' @importFrom dplyr select mutate distinct group_by ungroup left_join filter left_join summarise
#' @importFrom tidyr complete pivot_wider
#' @importFrom scales hue_pal
#' @importFrom magrittr "%>%"
#' @importFrom ComplexHeatmap pheatmap
hap_read_counts_plot <- function(seq_flt_tbl) {
  suppressMessages(pheat_data <- seq_flt_tbl %>%
                     group_by(sample, marker_id) %>%
                     summarise(
                       count = log10(sum(count) + 1)
                     ) %>%
                     ungroup() %>%
                     pivot_wider(
                       names_from = marker_id,
                       values_from = count
                     ) %>% as.data.frame())
  
  rownames(pheat_data) <- pheat_data$sample
  pheat_data <- pheat_data[, -1]
  # Annotation color
  annotation <- data.frame(sample = rownames(pheat_data))
  annotation <- annotation %>%
    left_join(
      (sample_manifest %>%
         select(sample, info)),
      by = "sample"
    ) %>%
    as.data.frame()
  rownames(annotation) <- annotation$sample
  rowname_anno <- rownames(annotation)
  annotation <- data.frame(info = annotation$info)
  rownames(annotation) <- rowname_anno
  
  # Color
  suppressMessages(pal_info <- annotation %>%
                     select(info) %>%
                     distinct() %>%
                     mutate(colour = scales::hue_pal()(n())) %>%
                     with(setNames(colour, info)))
  
  pal_info <- list(
    info = pal_info
  )
  
  pheat_data[is.na(pheat_data)] <- 0
  col <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
  col <- c("#BEBEBE", col)
  
  ComplexHeatmap::pheatmap(t(pheat_data),
                           name = "Log10(reads + 1)",
                           row_title = "marker",
                           color = col,
                           column_title = "sample",
                           annotation_col = annotation,
                           annotation_colors = pal_info
  )
}
#' @importFrom dplyr select mutate distinct group_by ungroup left_join filter left_join summarise
#' @importFrom tidyr pivot_wider
#' @importFrom magrittr "%>%"
#' @importFrom heatmaply heatmaply
unique_haplotypes_heatmap <- function(data) {
  suppressMessages(uni_hap <- data %>%
                     group_by(sample, marker_id) %>%
                     summarise(uni_hap = n()) %>%
                     ungroup() %>%
                     pivot_wider(
                       names_from = marker_id,
                       values_from = uni_hap
                     ) %>%
                     as.data.frame())
  rownames(uni_hap) <- uni_hap$sample
  uni_hap <- uni_hap[, -1]
  uni_hap[is.na(uni_hap)] <- 0
  
  col <- viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis")
  col <- c("#BEBEBE", col)
  
  heatmaply(t(uni_hap),
            key.title = "Number of distinct haplotypes",
            colors = col,
            xlab = "sample_id",
            ylab = "marker_id",
            main = str_c("Sample Number: ", nrow(uni_hap)),
            column_text_angle = 90
  )
}
#' @importFrom dplyr select mutate distinct group_by ungroup left_join filter left_join summarise rename across
#' @importFrom tidyr pivot_wider
#' @importFrom magrittr "%>%"
#' @importFrom stringr str_length
#' @importFrom flextable flextable theme_vanilla
#' @importFrom purrr map_chr
#' @importFrom digest digest
hapdiv_metrics <- function(data, marker_info) {
  # Number of haplotype
  suppressMessages(haplotype <- data %>%
                     group_by(marker_id) %>%
                     summarise(haplotype = length(unique(haplotype))))
  # He
  suppressMessages(heter <- data %>%
                     group_by(marker_id, haplotype) %>%
                     summarise(n = n()) %>%
                     ungroup() %>%
                     group_by(marker_id) %>%
                     mutate(freq = n / sum(n)) %>%
                     ungroup() %>%
                     group_by(marker_id) %>%
                     mutate(Number = sum(n)) %>%
                     mutate(scale_N = Number / (Number - 1)) %>%
                     summarise(He = scale_N * (1 - sum(freq^2))) %>%
                     distinct() %>%
                     ungroup())
  
  # Marker length
  length <- marker_info
  length <- length %>%
    mutate(Size = str_length(seq)) %>%
    dplyr::select(marker_id, Size) %>%
    dplyr::rename("Size(bp)" = Size)
  
  
  # SNPs
  # Get the variant position
  all_markers <- marker_info
  all_markers <- all_markers %>%
    filter(marker_id %in% as.character(data$marker_id))
  seq_tbl_all <- data %>%
    mutate(seq_id = map_chr(sequence, digest::digest))
  vars <- call_vars(
    seq_tbl = seq_tbl_all,
    marker_info = all_markers
  )
  vars <- vars %>%
    dplyr::select(marker_id, var_id, type)
  snps <- vars
  suppressMessages(snps <- snps %>%
                     group_by(marker_id) %>%
                     summarise(SNPs = n()) %>%
                     ungroup())
  
  # Average haplotype
  suppressMessages(ave_hap <- data %>%
                     group_by(sample_id, marker_id) %>%
                     summarise(n = n()) %>%
                     ungroup() %>%
                     group_by(marker_id) %>%
                     summarise(Ave_hap = mean(n)) %>%
                     ungroup())
  
  # Full table
  Output_table <- haplotype %>%
    left_join(heter, by = "marker_id") %>%
    left_join(length, by = "marker_id") %>%
    left_join(snps, by = "marker_id") %>%
    left_join(ave_hap, by = "marker_id") %>%
    mutate(across(where(is.numeric), ~ round(., 2)))
  colnames(Output_table)[1] <- "Marker"
  # Format table
  Output_table
  flextable(Output_table) %>% theme_vanilla()
}
#' @importFrom dplyr select mutate distinct group_by ungroup left_join full_join filter summarise rename across arrange
#' @importFrom magrittr "%>%"
#' @importFrom Biostrings DNAStringSet
#' @importFrom DECIPHER BrowseSeqs AlignSeqs
#' @importFrom purrr map
#' @importFrom htmltools includeHTML tagList
haplotype_sequence_vis <- function(data, marker_info) {
  suppressMessages(plot <- data %>%
                     group_by(marker_id, sequence, haplotype) %>%
                     summarise(count = n()) %>%
                     arrange(desc(count)) %>%
                     ungroup() %>%
                     group_by(marker_id) %>%
                     mutate(frequency = round(count / sum(count) * 100, 1)) %>%
                     ungroup() %>%
                     arrange(marker_id))
  
  
  plot_split <- plot %>%
    split(plot$marker_id)
  
  plot_split <- lapply(seq_along(plot_split), function(i) {
    plot_split[[i]] %>%
      mutate(Label = paste0(haplotype, "(", frequency, "%)", sep = ""))
  })
  
  
  reference <- marker_info
  reference <- reference %>%
    dplyr::select(marker_id, seq)
  colnames(reference) <- c("marker_id", "sequence")
  reference <- reference %>%
    mutate(Label = "Reference")
  reference <- reference %>%
    filter(marker_id %in% unique(plot$marker_id))
  
  plot_split <- do.call(rbind, plot_split)
  
  plot_split <- plot_split %>%
    full_join(reference, by = c("marker_id", "sequence", "Label"))
  
  plot_split <- plot_split %>%
    split(plot_split$marker_id)
  
  
  seq_sequence <- lapply(seq_along(plot_split), function(i) {
    seq <- DNAStringSet(plot_split[[i]]$sequence)
    names(seq) <- plot_split[[i]]$Label
    seq
  })
  
  names(seq_sequence) <- names(plot_split)
  
  alignment_html <-
    map(seq_sequence, function(al) {
      fn <- tempfile(fileext = ".html")
      DECIPHER::BrowseSeqs(AlignSeqs(al, verbose = FALSE), highlight = 0, htmlFile = fn, openURL = FALSE)
      include <- htmltools::includeHTML(fn)
      file.remove(fn)
      return(include)
    })
  
  htmltools::tagList(alignment_html)
}
#' @importFrom dplyr select mutate distinct group_by ungroup left_join full_join filter summarise rename across arrange n_distinct first
#' @importFrom magrittr "%>%"
#' @importFrom Biostrings DNAStringSet
#' @importFrom DECIPHER BrowseSeqs AlignSeqs
#' @importFrom purrr map
#' @importFrom tibble as_tibble
#' @importFrom stats setNames
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stringr str_c
#' @importFrom magrittr set_rownames
#' @importFrom ComplexHeatmap pheatmap
hap_seq_sim <- function(data, marker_info) {
  suppressMessages(plot <- data %>%
                     group_by(marker_id, sequence, haplotype) %>%
                     summarise(count = n()) %>%
                     arrange(desc(count)) %>%
                     ungroup() %>%
                     group_by(marker_id) %>%
                     mutate(frequency = round(count / sum(count) * 100, 1)) %>%
                     ungroup() %>%
                     arrange(marker_id))
  
  suppressMessages(one_sequence_marker <- data %>%
                     group_by(marker_id, sequence) %>%
                     summarise(count = n()) %>%
                     ungroup() %>%
                     group_by(marker_id) %>%
                     summarise(count = n()) %>%
                     ungroup() %>%
                     filter(count > 1))
  
  plot <- plot %>%
    filter(marker_id %in% unique(one_sequence_marker$marker_id))
  
  plot_split <- plot %>%
    split(plot$marker_id)
  
  reference <- marker_info
  reference <- reference %>%
    dplyr::select(marker_id, seq)
  colnames(reference) <- c("marker_id", "sequence")
  reference <- reference %>%
    mutate(haplotype = "Reference")
  reference <- reference %>%
    filter(marker_id %in% unique(plot$marker_id))
  
  plot_split <- do.call(rbind, plot_split)
  
  plot_split <- plot_split %>%
    full_join(reference, by = c("marker_id", "sequence", "haplotype"))
  
  maj_alignments <-
    plot_split %>%
    split.data.frame(.$marker_id) %>%
    map(function(data) {
      DECIPHER::AlignSeqs(Biostrings::DNAStringSet(data$sequence), processors = 4, verbose = FALSE) %>%
        setNames(data$haplotype)
    })
  
  sample_genotypes <- lapply(seq_along(maj_alignments), function(i) {
    sample_genotypes <-
      maj_alignments[[i]] %>%
      as.matrix() %>%
      t() %>%
      as_tibble() %>%
      mutate(pos = seq_len(n())) %>%
      mutate(marker = names(maj_alignments)[[i]]) %>%
      pivot_longer(c(-pos, -marker), names_to = "sample", values_to = "base") %>%
      group_by(sample) %>%
      filter(!any(is.na(base))) %>%
      # filter for biallelic
      group_by(marker, pos) %>%
      filter(n_distinct(base, na.rm = TRUE) > 1) %>%
      mutate(genotype = if_else(base == dplyr::first(na.omit(base)), 0L, 1L)) %>%
      ungroup() %>%
      mutate(variant = str_c(marker, "-", pos)) %>%
      dplyr::select(sample, variant, genotype) %>%
      pivot_wider(names_from = variant, values_from = genotype) %>%
      (function(x) {
        select(x, -sample) %>%
          as.matrix() %>%
          set_rownames(x$sample)
      })
  })
  
  cols <- c("white", "red")
  
  for (i in seq_along(sample_genotypes)) {
    ComplexHeatmap::pheatmap(sample_genotypes[[i]], main = names(maj_alignments)[[i]], legend = FALSE, color = cols, cellheight = 10, cluster_cols = FALSE, run_draw = TRUE)
  }
}
#' @importFrom dplyr select mutate distinct group_by ungroup left_join full_join filter summarise rename across arrange n_distinct first slice
#' @importFrom magrittr "%>%"
#' @importFrom Biostrings DNAStringSet
#' @importFrom DECIPHER BrowseSeqs AlignSeqs
#' @importFrom purrr map map2_df
#' @importFrom tibble as_tibble
#' @importFrom stats setNames
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stringr str_c
#' @importFrom magrittr set_rownames
#' @importFrom ape nj dist.gene
#' @importFrom tidyr complete
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom scales hue_pal
#' @importFrom ggtree ggtree geom_tiplab facet_plot facet_widths
#' @importFrom ggstance geom_barh
hap_tree_plot <- function(data) {
  ggtree_data <- data %>%
    group_by(sample, marker_id) %>%
    arrange(desc(count)) %>%
    dplyr::slice(1) %>%
    ungroup()
  ggtree_data$count[is.na(ggtree_data$count)] <- 10000
  si_hap_table <- ggtree_data %>%
    select(sample, marker_id, haplotype, count, freq = frequency, sequence)
  
  maj_alignments <-
    si_hap_table %>%
    split.data.frame(.$marker_id) %>%
    map(function(data) {
      DECIPHER::AlignSeqs(Biostrings::DNAStringSet(data$sequence), processors = 4, verbose = FALSE) %>%
        setNames(data$sample)
    })
  
  sample_genotypes <-
    maj_alignments %>%
    map2_df(., names(.), function(aln, marker) {
      as.matrix(aln) %>%
        t() %>%
        as_tibble() %>%
        mutate(pos = seq_len(n())) %>%
        mutate(marker = marker)
    }) %>%
    pivot_longer(c(-pos, -marker), names_to = "sample", values_to = "base") %>%
    group_by(sample) %>%
    filter(!any(is.na(base))) %>%
    # filter for biallelic
    group_by(marker, pos) %>%
    filter(n_distinct(base, na.rm = TRUE) == 2) %>%
    mutate(genotype = if_else(base == dplyr::first(na.omit(base)), 0L, 1L)) %>%
    ungroup() %>%
    mutate(variant = str_c(marker, "-", pos)) %>%
    dplyr::select(sample, variant, genotype) %>%
    pivot_wider(names_from = variant, values_from = genotype) %>%
    (function(x) {
      select(x, -sample) %>%
        as.matrix() %>%
        set_rownames(x$sample)
    })
  
  if (ncol(sample_genotypes) < 2) {
    print("This analysis requires that the main haplotype sequence in the dataset contains at least 2 variants.")
  } else {
    # Nj tree
    stree <- nj(dist.gene(sample_genotypes))
    
    # ggtree haplotype plot
    ggtree_data_two <- data %>%
      dplyr::select(sample, marker_id, sequence, count, haplotype)
    ggtree_data_two <- ggtree_data_two %>%
      group_by(sample, marker_id) %>%
      arrange(desc(count)) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      dplyr::select(sample, marker_id, haplotype)
    
    # ggtree_df
    ggtree_df <- data.frame(stree$tip.label, variable = NA, value = rep(1))
    ggtree_df <- ggtree_df %>%
      group_by(stree.tip.label) %>%
      complete(variable = unique(ggtree_data$marker_id), fill = list(variable = NA))
    ggtree_df$value <- 1
    ggtree_df <- ggtree_df %>%
      na.omit()
    colnames(ggtree_data_two) <- c("stree.tip.label", "variable", "marker_name")
    ggtree_df <- ggtree_df %>%
      left_join(ggtree_data_two, by = c("stree.tip.label", "variable"))
    ggtree_df <- ggtree_df %>%
      dplyr::select(stree.tip.label, marker_name, value, variable)
    colnames(ggtree_df)[4] <- "marker"
    colnames(ggtree_df)[2] <- "variable"
    
    suppressMessages(pal <- ggtree_df %>%
                       ungroup() %>%
                       select(variable) %>%
                       distinct() %>%
                       mutate(colour = distinctColorPalette(n())) %>%
                       with(setNames(colour, variable)))
    
    suppressMessages(pal_marker <- ggtree_df %>%
                       select(marker) %>%
                       distinct() %>%
                       mutate(colour = scales::hue_pal()(n())) %>%
                       with(setNames(colour, marker)))
    p <- ggtree(stree) +
      geom_tiplab(
        size = 3, align = TRUE,
        offset = .8, hjust = .5
      ) + theme_tree()
    
    p2 <- facet_plot(p + xlim_tree(8),
                     panel = "Major Haplotype", data = ggtree_df, geom = geom_barh,
                     mapping = aes(x = value, col = marker, fill = as.factor(variable)), alpha = 0.7, show.legend = TRUE,
                     width = 1, stat = "identity"
    ) + xlim_tree(9) +
      scale_fill_manual(values = pal, guide = "none") +
      scale_color_manual(values = pal_marker)
    facet_widths(p2, widths = c(3, 2))
  }
}
#' @importFrom dplyr select mutate distinct group_by ungroup left_join full_join filter summarise rename across arrange n_distinct first slice
#' @importFrom magrittr "%>%"
#' @importFrom Biostrings DNAStringSet
#' @importFrom DECIPHER BrowseSeqs AlignSeqs
#' @importFrom purrr map map2_df map_chr
#' @importFrom tibble as_tibble
#' @importFrom stats setNames
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stringr str_c
#' @importFrom magrittr set_rownames
#' @importFrom ape nj dist.gene
#' @importFrom tidyr complete
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom scales hue_pal
#' @importFrom ggtree ggtree geom_tiplab facet_plot facet_widths theme_tree xlim_tree
#' @importFrom ggstance geom_barh
#' @importFrom digest digest
hap_tree_plot_imputed <- function(data) {
  Mode <- function(x) {
    ux <- unique(x)
    theNAs <- is.na(ux)
    ux <- ux[!theNAs]
    ux[which.max(tabulate(match(x, ux)))]
  }
  haplotype <- data %>%
    mutate(seq_id = map_chr(sequence, digest::digest)) %>%
    group_by(sample, marker_id) %>%
    arrange(desc(frequency)) %>%
    dplyr::slice(1) %>%
    group_by(sample) %>%
    complete(marker_id = unique(data$marker_id), fill = list(marker_id = NA)) %>%
    ungroup() %>%
    group_by(marker_id) %>%
    mutate(
      haplotype = replace_na(haplotype, Mode(haplotype)),
      seq_id = replace_na(seq_id, Mode(seq_id)),
      sequence = replace_na(sequence, Mode(sequence))
    ) %>%
    ungroup()
  
  haplotype$haplotype[is.na(haplotype$count) == TRUE] <- "Inputed"
  
  ggtree_data <- haplotype %>%
    group_by(sample, marker_id) %>%
    arrange(desc(count)) %>%
    dplyr::slice(1) %>%
    ungroup()
  ggtree_data$count[is.na(ggtree_data$count)] <- 10000
  
  si_hap_table <- ggtree_data %>%
    select(sample, marker_id, haplotype, count, freq = frequency, sequence)
  
  
  maj_alignments <-
    si_hap_table %>%
    split.data.frame(.$marker_id) %>%
    map(function(data) {
      DECIPHER::AlignSeqs(Biostrings::DNAStringSet(data$sequence), processors = 4, verbose = FALSE) %>%
        setNames(data$sample)
    })
  
  sample_genotypes <-
    maj_alignments %>%
    map2_df(., names(.), function(aln, marker) {
      as.matrix(aln) %>%
        t() %>%
        as_tibble() %>%
        mutate(pos = seq_len(n())) %>%
        mutate(marker = marker)
    }) %>%
    pivot_longer(c(-pos, -marker), names_to = "sample", values_to = "base") %>%
    group_by(sample) %>%
    filter(!any(is.na(base))) %>%
    # filter for biallelic
    group_by(marker, pos) %>%
    filter(n_distinct(base, na.rm = TRUE) == 2) %>%
    mutate(genotype = if_else(base == dplyr::first(na.omit(base)), 0L, 1L)) %>%
    ungroup() %>%
    mutate(variant = str_c(marker, "-", pos)) %>%
    dplyr::select(sample, variant, genotype) %>%
    pivot_wider(names_from = variant, values_from = genotype) %>%
    (function(x) {
      select(x, -sample) %>%
        as.matrix() %>%
        set_rownames(x$sample)
    })
  
  
  if (ncol(sample_genotypes) < 2) {
    print("This analysis requires that the main haplotype sequence in the dataset contains at least 2 variants.")
  } else {
    # Nj tree
    stree <- nj(dist.gene(sample_genotypes))
    
    # ggtree haplotype plot
    ggtree_data_two <- data %>%
      dplyr::select(sample, marker_id, sequence, count, haplotype)
    ggtree_data_two <- ggtree_data_two %>%
      group_by(sample, marker_id) %>%
      arrange(desc(count)) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      dplyr::select(sample, marker_id, haplotype)
    
    # try
    ggtree_df <- data.frame(stree$tip.label, variable = NA, value = rep(1))
    ggtree_df <- ggtree_df %>%
      group_by(stree.tip.label) %>%
      complete(variable = unique(ggtree_data$marker_id), fill = list(variable = NA))
    ggtree_df$value <- 1
    ggtree_df <- ggtree_df %>%
      na.omit()
    colnames(ggtree_data_two) <- c("stree.tip.label", "variable", "marker_name")
    ggtree_df <- ggtree_df %>%
      left_join(ggtree_data_two, by = c("stree.tip.label", "variable")) %>%
      dplyr::select(stree.tip.label, marker_name, value, variable)
    colnames(ggtree_df)[4] <- "marker"
    colnames(ggtree_df)[2] <- "variable"
    ggtree_df$variable[is.na(ggtree_df$variable)] <- "Inputed"
    
    
    suppressMessages(pal <- ggtree_df %>%
                       ungroup() %>%
                       select(variable) %>%
                       distinct() %>%
                       filter(variable != "Inputed") %>%
                       mutate(colour = distinctColorPalette(n())) %>%
                       add_row(variable = "Inputed", colour = "White") %>%
                       with(setNames(colour, variable)))
    
    suppressMessages(pal_marker <- ggtree_df %>%
                       select(marker) %>%
                       distinct() %>%
                       mutate(colour = scales::hue_pal()(n())) %>%
                       with(setNames(colour, marker)))
    
    
    p <- ggtree(stree) +
      geom_tiplab(
        size = 3, align = TRUE,
        offset = .8, hjust = .5
      ) + theme_tree()
    
    
    ggtree_df$variable <- factor(ggtree_df$variable,
                                 levels = names(pal)
    )
    p2 <- facet_plot(p + xlim_tree(8),
                     panel = "Major Haplotype", data = ggtree_df, geom = geom_barh,
                     mapping = aes(x = value, col = marker, fill = variable), alpha = 0.7, show.legend = TRUE,
                     width = 1, stat = "identity"
    ) + xlim_tree(9) +
      scale_fill_manual(values = pal, guide = "none") +
      scale_color_manual(values = pal_marker)
    
    facet_widths(p2, widths = c(3, 2))
  }
}
#' @importFrom dplyr select mutate distinct group_by ungroup left_join full_join filter summarise rename across arrange n_distinct first slice bind_cols
#' @importFrom magrittr "%>%"
#' @importFrom Biostrings DNAStringSet
#' @importFrom DECIPHER BrowseSeqs AlignSeqs
#' @importFrom purrr map map2_df map_chr
#' @importFrom tibble as_tibble
#' @importFrom stats setNames prcomp
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stringr str_c
#' @importFrom magrittr set_rownames
#' @importFrom ape nj
#' @importFrom tidyr complete
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom scales hue_pal
#' @importFrom plotly ggplotly
hap_pca <- function(data) {
  pca_data <- data
  sample_info <- data %>%
    select(sample_id, sample, info) %>%
    distinct()
  ggtree_data <- pca_data %>%
    group_by(sample_id, marker_id) %>%
    arrange(desc(count)) %>%
    dplyr::slice(1) %>%
    ungroup()
  ggtree_data$count[is.na(ggtree_data$count)] <- 10000
  
  si_hap_table <- ggtree_data %>%
    group_by(sample_id, marker_id) %>%
    mutate(freq = count / sum(count)) %>%
    ungroup() %>%
    select(sample_id, marker_id, haplotype, count, freq, sequence)
  
  maj_alignments <-
    si_hap_table %>%
    group_by(sample_id, marker_id) %>%
    dplyr::slice(which.max(freq)) %>%
    ungroup() %>%
    split.data.frame(.$marker_id) %>%
    map(function(data) {
      DECIPHER::AlignSeqs(Biostrings::DNAStringSet(data$sequence), processors = 4, verbose = FALSE) %>%
        setNames(data$sample_id)
    })
  
  sample_genotypes <-
    maj_alignments %>%
    map2_df(., names(.), function(aln, marker) {
      as.matrix(aln) %>%
        t() %>%
        as_tibble() %>%
        mutate(pos = seq_len(n())) %>%
        mutate(marker = marker)
    }) %>%
    pivot_longer(c(-pos, -marker), names_to = "sample", values_to = "base") %>%
    group_by(sample) %>%
    filter(!any(is.na(base))) %>%
    # filter for biallelic
    group_by(marker, pos) %>%
    filter(n_distinct(base, na.rm = TRUE) == 2) %>%
    mutate(genotype = if_else(base == dplyr::first(na.omit(base)), 0L, 1L)) %>%
    ungroup() %>%
    mutate(variant = str_c(marker, "-", pos)) %>%
    dplyr::select(sample, variant, genotype) %>%
    pivot_wider(names_from = variant, values_from = genotype) %>%
    (function(x) {
      select(x, -sample) %>%
        as.matrix() %>%
        set_rownames(x$sample)
    })
  
  if (ncol(sample_genotypes) < 2) {
    print("This analysis requires that the main haplotype sequence in the dataset contains at least 2 variants.")
  } else {
    hap_pca <-
      prcomp(sample_genotypes) %>%
      {
        as_tibble(.$x)
      } %>%
      bind_cols(sample_id = rownames(sample_genotypes), .)
    
    p <- hap_pca %>%
      left_join(sample_info, by = "sample_id") %>%
      ggplot(aes(PC1, PC2, label = sample, col = info)) +
      geom_jitter(alpha = 0.5) +
      theme_bw()
    
    ggplotly(p, tooltip = c("sample"))
  }
}
#' @importFrom dplyr select mutate distinct group_by ungroup left_join full_join filter summarise rename across arrange n_distinct first slice bind_cols
#' @importFrom magrittr "%>%"
#' @importFrom Biostrings DNAStringSet
#' @importFrom DECIPHER BrowseSeqs AlignSeqs
#' @importFrom purrr map map2_df map_chr
#' @importFrom tibble as_tibble
#' @importFrom stats setNames prcomp
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stringr str_c
#' @importFrom magrittr set_rownames
#' @importFrom ape nj
#' @importFrom tidyr complete
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom scales hue_pal
#' @importFrom plotly ggplotly
hap_pca_imputed <- function(data) {
  seq_tbl_id_input <- data
  sample_info <- data %>%
    select(sample_id, sample, info) %>%
    distinct()
  Mode <- function(x) {
    ux <- unique(x)
    theNAs <- is.na(ux)
    ux <- ux[!theNAs]
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  haplotype <- seq_tbl_id_input %>%
    mutate(seq_id = map_chr(sequence, digest::digest)) %>%
    group_by(sample_id, marker_id) %>%
    arrange(desc(frequency)) %>%
    dplyr::slice(1) %>%
    group_by(sample_id) %>%
    complete(marker_id = unique(seq_tbl_id_input$marker_id), fill = list(marker_id = NA)) %>%
    ungroup() %>%
    group_by(marker_id) %>%
    mutate(
      haplotype = replace_na(haplotype, Mode(haplotype)),
      seq_id = replace_na(seq_id, Mode(seq_id)),
      sequence = replace_na(sequence, Mode(sequence))
    ) %>%
    ungroup()
  
  inputed_sample <- haplotype %>%
    filter(is.na(count) == TRUE)
  inputed_sample <- unique(inputed_sample$sample_id)
  
  
  pca_data <- haplotype
  ggtree_data <- pca_data
  ggtree_data$count[is.na(ggtree_data$count)] <- 10000
  
  si_hap_table <- ggtree_data %>%
    select(sample_id, marker_id, haplotype, count, freq = frequency, sequence)
  
  maj_alignments <-
    si_hap_table %>%
    split.data.frame(.$marker_id) %>%
    map(function(data) {
      DECIPHER::AlignSeqs(Biostrings::DNAStringSet(data$sequence), processors = 4, verbose = FALSE) %>%
        setNames(data$sample_id)
    })
  
  sample_genotypes <-
    maj_alignments %>%
    map2_df(., names(.), function(aln, marker) {
      as.matrix(aln) %>%
        t() %>%
        as_tibble() %>%
        mutate(pos = seq_len(n())) %>%
        mutate(marker = marker)
    }) %>%
    pivot_longer(c(-pos, -marker), names_to = "sample", values_to = "base") %>%
    group_by(sample) %>%
    filter(!any(is.na(base))) %>%
    # filter for biallelic
    group_by(marker, pos) %>%
    filter(n_distinct(base, na.rm = TRUE) == 2) %>%
    mutate(genotype = if_else(base == dplyr::first(na.omit(base)), 0L, 1L)) %>%
    ungroup() %>%
    mutate(variant = str_c(marker, "-", pos)) %>%
    dplyr::select(sample, variant, genotype) %>%
    pivot_wider(names_from = variant, values_from = genotype) %>%
    (function(x) {
      select(x, -sample) %>%
        as.matrix() %>%
        set_rownames(x$sample)
    })
  
  
  if (ncol(sample_genotypes) < 2) {
    print("This analysis requires that the main haplotype sequence in the dataset contains at least 2 variants.")
  } else {
    hap_pca <-
      prcomp(sample_genotypes) %>%
      {
        as_tibble(.$x)
      } %>%
      bind_cols(sample_id = rownames(sample_genotypes), .)
    
    p <- hap_pca %>%
      mutate(Type = ifelse(sample_id %in% inputed_sample, "Missing-data-imputation", "NonMissing")) %>%
      left_join(sample_info, by = "sample_id") %>%
      mutate(info = ifelse(Type == "Missing-data-imputation", "Missing-data-imputation", info)) %>%
      ggplot(aes(PC1, PC2, label = sample_id, col = info)) +
      geom_jitter(alpha = 0.5) +
      theme_bw() +
      theme(
        legend.position = "bottom"
      )
    
    ggplotly(p, tooltip = c("sample"))
  }
}
