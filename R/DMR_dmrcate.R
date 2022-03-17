#' Detect differentially methylated positions in the DNA using the limma package.
#'
#' This function uses `limma` package to generate t-values. Then it smooths out these t-statistics using gaussian kernel smoothing.
#' The t-values that remain significant after smoothing and multiple testing, are grouped when they are within Lambda base pairs from each other.
#' Such a group is called a Differently Methylated Region (DMR).
#'
#' @param GRcset Preprocessed Green-Red methylation set
#' @param DMP Output from `DMP_limma()` function.
#' @param lambda Gaussian kernel bandwidth for smoothed-function estimation. Also informs DMR bookend definition; gaps >= lambda between significant CpG sites will be in separate DMRs. Support is truncated at 5*lambda. Default is 1000 nucleotides.
#' @param fdr false discovery rate, default = 0.05.
#'
#' @return A a S4 object with statistics for every DMR
#' \itemize{
#'    \item DMR - `dmrcate()`output, a S4 object with statistics for every DMR
#'    \item DMRRanges - `extractRanges()` output, a S4 object with Ranges from found DMRs
#'    \item annotation - `cpg.annotate()` output, a s4 object with statistics on every CpG site.
#' }
#' @export
#'
#' @examples
DMR_dmrcate <- function(GRcset, DMP, lambda=1000, fdr=0.05){
  annotation <- cpg.annotate(object = GRcset, datatype = "array", what = "M",
                             analysis.type = "differential", design = DMP$design,
                             contrasts = TRUE, cont.matrix = DMP$contMatrix,
                             coef = colnames(DMP$contMatrix)[1], arraytype = "EPIC", fdr = fdr)
  DMR <- dmrcate(annotation, lambda=lambda, C=2)
  DMRRanges <- extractRanges(DMR)

  return(list(DMR=DMR, DMRRanges=DMRRanges, annotation=annotation))
}
