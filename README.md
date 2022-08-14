# AmpSeqR: an R package for amplicon deep sequencing data analysis
https://github.com/bahlolab/ampseqr

AmpseqR is an R package for analysis of amplicon deep sequencing (AmpSeq) data generated on the Illumina platform. The pipeline offers various useful functions including Data pre-processing, Amplicon sequence variant (ASVs) estimation, Data post-processing, and Data visualization. Additionally, AmpseqR includes several parameters to filter noise reads and improve the accuracy of the detected haplotype.

## Installation

The AmpSeqR currently available to install from Github:

```r
# install using devtools packages
# first install devtools dependencies
if (!require(devtools)) install.packages("devtools")
devtools::install_github("bahlolab/AmpSeqR")
```

## What data input does _AmpSeqR_ require?

Inputs are the standard paired-end FASTQ format provided by the common Illumina sequencing platforms (e.g., MiSeq), as well as sample barcodes and target amplicon details.

Example data:
```{r}
library(AmpSeqR)
example_data <- get_ampseqr_example_data()
```

## How do I use _AmpSeqR_?
See the [introduction vignette](http://bahlolab.github.io/AmpSeqR/vignettes/AmpSeqR.Rmd) for usage examples.
