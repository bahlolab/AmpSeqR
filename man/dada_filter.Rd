% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dada2_wrappers.R
\name{dada_filter}
\alias{dada_filter}
\title{dada_filter}
\usage{
dada_filter(
  read_table,
  output_dir,
  output_sub_dir,
  max_exp_errors = 1L,
  min_len = 25L,
  marker_trim = NULL,
  threads = 1L
)
}
\arguments{
\item{read_table}{(Required). The demultiplexed read table includes: sample_id, marker_id, reads_1 (the forward read fastq file path), reads_2 (the reverse read fastq file path), n (number of demultiplexed reads), sample, info.}

\item{output_dir}{(Required). The path to the output filtered table file.}

\item{output_sub_dir}{(Required). The path to the output filtered files.}

\item{max_exp_errors}{(Required). The same as in dada2::filterAndTrim().}

\item{min_len}{(Optional). Default 25 (numeric vector). Remove reads with length less than minLen.}

\item{marker_trim}{(Optional). Default NULL. If trimming is to be performed, marker_trim should be a data frame with three columns marker_id (character), trim_fwd (integer), and trim_rev (integer).}

\item{threads}{(Optional). Default 1. If an integer is provided, it is passed to the mc.cores argument of mclapply.}
}
\value{
filter folder: filtered and trimmed paired-end fastq files.

filtered_reads.rds table: the filtered table in RDS format which includes sample_id, marker_id, n (number of demultiplexed reads), sample, info, reads_1 (the filtered and trimmed forward fastq file path), reads_2 (the filtered and trimmed reverse fastq file path), n_in (number of reads before filtering and trimming), n_out (number of reads after filtering and trimming).
}
\description{
Translate dada2::filterAndTrim to AmpSeqR::dada_filter. Filter the poor-quality reads and trim low-quality bases.
}
