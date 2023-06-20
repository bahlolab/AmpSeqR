#' @export
#' @importFrom dplyr mutate_if
#' @importFrom ggplot2 aes
summarise_read_quality <- function(fwd_fqs, rev_fqs, num = 50000) {
  stopifnot(
    all(file.exists(fwd_fqs)),
    all(file.exists(rev_fqs))
  )
  smry1 <-
    bind_rows(
      ShortRead::qa(fwd_fqs, n = num)[["perCycle"]]$quality %>%
        mutate_if(is.factor, as.character) %>%
        mutate(Direction = "forward"),
      ShortRead::qa(rev_fqs, n = num)[["perCycle"]]$quality %>%
        mutate_if(is.factor, as.character) %>%
        mutate(Direction = "reverse")
    ) %>%
    group_by(Direction, Cycle, Score) %>%
    summarise(Count = sum(Count, na.rm = T)) %>%
    group_by(Direction, Cycle) %>%
    mutate(Prop = Count / sum(Count, na.rm = T)) %>%
    ungroup()

  smry2 <-
    smry1 %>%
    group_by(Direction, Cycle) %>%
    summarise(
      Total = sum(Count, na.rm = T),
      MeanScore = sum(Score * Count, na.rm = T) / Total
    ) %>%
    ungroup() %>%
    left_join(
      smry1 %>%
        filter(Cycle == 1) %>%
        group_by(Direction) %>%
        summarise(nReads = sum(Count)) %>%
        ungroup(),
      "Direction"
    ) %>%
    mutate(PropReadLen = Total / nReads)

  list(denisty_by_pos_and_qual = smry1, summary_by_pos = smry2)
}
#' @export
#' @importFrom dplyr arrange desc filter "%>%" summarise group_by mutate ungroup bind_rows left_join slice rename
#' @importFrom tidyr nest unnest spread
plot_quality_profile <- function(reads_fwd, reads_rev, title) {
  quality_smry <- summarise_read_quality(fwd_fqs = reads_fwd, rev_fqs = reads_rev)
  .plot_quality_profile(
    quality_smry$denisty_by_pos_and_qual,
    quality_smry$summary_by_pos,
    title
  )
}

#' @importFrom ggplot2 ggplot geom_tile geom_line ylab scale_fill_continuous facet_wrap theme element_blank
#' @importFrom cowplot plot_grid ggdraw draw_label
#' @importFrom dplyr "%>%"
.plot_quality_profile <- function(denisty_by_pos_and_qual, summary_by_pos, title) {
  stopifnot(
    is.data.frame(denisty_by_pos_and_qual),
    is.data.frame(summary_by_pos),
    is_string(title)
  )
  p1 <-
    denisty_by_pos_and_qual %>%
    ggplot() +
    geom_tile(aes(Cycle, Score, fill = Prop)) +
    geom_line(data = summary_by_pos, aes(Cycle, MeanScore), col = "red") +
    ylab("Quality Score") +
    scale_fill_continuous(high = "gray15", low = "gray75", limits = c(0, 1)) +
    facet_wrap(~Direction, nrow = 1) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())

  p2 <-
    summary_by_pos %>%
    ggplot() +
    geom_tile(aes(y = "PropReads", x = Cycle, fill = PropReadLen)) +
    scale_fill_continuous(high = "gray15", low = "gray75", limits = c(0, 1)) +
    facet_wrap(~Direction, nrow = 1) +
    theme(
      axis.title.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none",
      strip.background = element_blank(), strip.text = element_blank()
    )


  cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_label(paste0("Quality profile - ", title)),
    p1, p2,
    ncol = 1, rel_heights = c(0.05, 0.80, 0.15), align = "v", axis = "lr"
  )
}


#' Plot quality profile of fastq files
#'
#' Visualize the quality profiles of the forward and reverse sequence position.
#'
#' @param data (Required). The demultiplexed read table includes: sample_id, marker_id, reads_1 (the forward read fastq file path), reads_2 (the reverse read fastq file path), n (number of demultiplexed reads), sample, info.
#'
#' @return
#'
#' Generate a visual summary plot of the distribution of quality scores of forward and reverse sequence position for all fastq files.
#'
#' @export
#'
#' @importFrom dplyr arrange desc filter "%>%" summarise group_by mutate ungroup bind_rows left_join slice rename
#' @importFrom tidyr nest unnest spread
#' @importFrom furrr future_map
plot_quality <- function(data) {
  data %>%
    filter(!is.na(sample_id), !is.na(marker_id), n > 0) %>%
    split.data.frame(.$marker_id) %>%
    furrr::future_map(function(mkr_tab) {
      smry <- summarise_read_quality(
        fwd_fqs = mkr_tab$reads_1,
        rev_fqs = mkr_tab$reads_2
      )
      plot <- plot_quality_profile(mkr_tab$reads_1, mkr_tab$reads_2, unique(mkr_tab$marker_id))
      plot
    })
}
