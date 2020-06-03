#' Title
#'
#' @param x
#' @param badSamples
#' @param n
#' @param ylim
#' @param showPDF
#'
#' @return
#' @export
#'

FindBadSamples <- function(x, badSamples=c(), n=4, ylim=c(0,5), showPDF=TRUE){
  if(is.null(x$bset_raw)){
    betas <- preprocessRaw(x$rgset) %>% getBeta()
  }else{
    betas <- x$bset_raw
  }
  x$bset_raw <- betas
  if(!length(badSamples)){
    badSamplesNum <- 1:ncol(betas)
    pdf_filepath <- getPath(x, "FindBadSamples_All", "output/intersampleQC", ".pdf")
    x$excludedSamples$df_badSamples <- data.frame()
  }else{
    pdf_filepath <- getPath(x, "FindBadSamples_Selection", "output/intersampleQC", ".pdf")
    if(is.numeric(badSamples)) badSamplesNum <- badSamples
    if(is.character(badSamples)) badSamplesNum <- which(colnames(betas) %in% badSamples)
    badSamplesNames <- colnames(betas)[badSamplesNum]
    badSamplesSwabnames <- unlist(lapply(badSamplesNames, function(x) strsplit(x,"_")[[1]][[3]]))
    x$excludedSamples$df_badSamples <- data.frame(Number = badSamplesNum, Name = badSamplesNames, SwabName = badSamplesSwabnames)
    write.csv(x$excludedSamples$df_badSamples, getPath(x, "ExcludedSwabs_ByQC","output/intersampleQC", ".csv"))
  }

  pdf(pdf_filepath, onefile = TRUE)
  if(length(badSamplesNum) >= n){
    for(i in 1:floor(length(badSamplesNum)/n)){
      densityPlot(betas[,badSamplesNum[c((i*n-n+1):(i*n))]],
                  main=i,
                  sampGroups = badSamplesNum[c((i*n-n+1):(i*n))],
                  ylim = ylim)
      invisible(gc())
    }
  }
  if(length(badSamplesNum)%%n){
    densityPlot(betas[,badSamplesNum[c(length(badSamplesNum)-length(badSamplesNum)%%n+1):length(badSamplesNum)]],
                main="remainder",
                sampGroups = badSamplesNum[c(length(badSamplesNum)-length(badSamplesNum)%%n+1):length(badSamplesNum)],
                ylim = ylim)
  }
  dev.off()
  if(showPDF) try(openPDF(normalizePath(pdf_filepath)))
  invisible(x)
}



#' Title
#'
#' @param x
#' @param badSamplesNames
#' @param badSamplesNum
#'
#' @return
#' @export
#'
#' @examples
RemoveBadSamples <- function(x, badSamplesNames=c(), badSamplesNum=c()){
  badSamplesNum <- c(which(colnames(rgset) %in% badSamplesNames), badSamplesNum) %>% unique()
  x <- FindBadSamples(x, badSamplesNum,showPDF = FALSE)

  if(length(badSamplesNum)){
    badSamplesSwabnames <- x$excludedSamples$df_badSamples$SwabName
    pdf(getPath(x,"DensityPlots_QC", "output/intersampleQC",".pdf"), onefile = TRUE)

    sampGroups <- rep('Valid',ncol(x$bset_raw))
    sampGroups[badSamplesNum] <- "Invalid"
    densityPlot(x$bset_raw, xlim=c(-0.1,1.1), ylim = c(0,7), main="All Samples", sampGroups = sampGroups, pal=c("firebrick3", "forestgreen"))

    x$rgset <- x$rgset[, -badSamplesNum]
    x$bset_raw <- x$bset_raw[, -badSamplesNum]
    x$params <- list(interSampleQC=TRUE)
    densityPlot(x$bset_raw, main="Valid Samples", pal="darkslateblue", legend = FALSE)

    dev.off()
  }else{
    pdf(getPath(x,"DensityPlot_All", "output/intersampleQC",".pdf"), onefile = TRUE)
    densityPlot(x$bset_raw, main="Density plot of All Samples", pal="darkslateblue", legend = FALSE)
    dev.off()
  }
  invisible(x)
}
