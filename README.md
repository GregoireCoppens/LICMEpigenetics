
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LICMEpigenetics

<!-- badges: start -->
<!-- badges: end -->

Provides easy to use, objective oriented functions for preprocessing
methylation data produced by an Illumina Infinium BeadChip and detecting
differentially methylated positions and regions within the DNA.

## Installation

You can install the released version of LICMEpigenetics from
[Github](https://github.com/GregoireCoppens/LICMEpigenetics) with:

``` r
#auth_token expires 06/04/2022

library(devtools)
devtools::install_github("GregoireCoppens/LICMEpigenetics", ref="master", auth_token="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")
```

## Example

This is a basic example which shows you how to preprocess methylation data starting from idat files and find the differentially methylated positions (DMPs):

``` r
library(LICMEpigenetics)
library(dplyr)

# Setup
set.seed(99999)
base_dir <- system.file("extdata", package="minfiDataEPIC") # path to idat files or folders
cat_vars = c("cat_1") # generic categorical name

# Preprocessing methylation data
preprocessed_set <- Make_RGset(base_dir, name="test") %>% # import idats
  ProbeNormalisation() %>% # Functional normalisation using funnorm
  ProbeExclusion() %>% # Exclude CPG, SBE, sex probes and probes not exceeding background signal
  ConvertSet() # convert genomic methyl set into beta set.

# DMP calculation
Random_AdjustmentVars <- data.frame(replicate(1, sample(0:1, 3, rep=TRUE)))
colnames(Random_AdjustmentVars) <- cat_vars

LimmaOutput <- DMP_limma(preprocessed_set$mset, Random_AdjustmentVars, cat_vars=cat_vars, Group="cat_1")
```
