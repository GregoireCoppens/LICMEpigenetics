#' @title FindBadSamples function
#' @description  Function to visually find bad samples with non-bimodal distributions using density plots.
#'
#' @param x Epigenetic set object.
#' @param badSamples Number or filenames of the bad idat files.
#' @param n number of density plots per plot. Max 4 is suggested.
#' @param ylim range of the y-axis
#' @param Swabname_position Position of swabname in the filename (separated by "_")
#'
#' @return Updated Epigenetic set object.
#'
#' @import minfiDataEPIC
#' @importFrom minfi densityPlot
#' @importFrom  minfi preprocessRaw
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#'
#' @export

FindBadSamples <- function(x, badSamples=c(), n=4, ylim=c(0,5), Swabname_position=3){
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

    # To do: Make more universal
    if(all(grepl("_", badSamplesNames))){
      badSamplesSwabnames <- unlist(lapply(badSamplesNames, function(x) strsplit(x,"_")[[1]][[Swabname_position]]))
    }else{
      badSamplesSwabnames <- badSamplesNames
    }
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
    if(length(badSamples)>n){
      densityPlot(betas[,badSamplesNum[c(length(badSamplesNum)-length(badSamplesNum)%%n+1):length(badSamplesNum)]],
                  main="remainder",
                  sampGroups = badSamplesNum[c(length(badSamplesNum)-length(badSamplesNum)%%n+1):length(badSamplesNum)],
                  ylim = ylim)
    } else {
      densityPlot(as.matrix(betas[,badSamplesNum]),
                  sampGroups = badSamplesNum,
                  ylim = ylim)
    }
  }
  dev.off()
  invisible(x)
}



#' @title RemoveBadSamples function
#' @description Function to remove bad samples from an epigenetic set. This function also outputs a pdf with density plots highlighting the removed samples and without the bad samples.
#' @param x Epigenetic set object.
#' @param badSamplesNames Filenames from bad samples
#' @param badSamplesNum Index numbers from bad samples.
#' @param Swabname_position Position of swabname in the filename (separated by "_")
#'
#' @return Updated Epigenetic set object.
#'
#' @importFrom minfi densityPlot
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#'
#' @export
#'
RemoveBadSamples <- function(x, badSamplesNames=c(), badSamplesNum=c(), Swabname_position=3){
  badSamplesNum <- c(which(colnames(x$rgset) %in% badSamplesNames), badSamplesNum) %>% unique()
  x <- FindBadSamples(x, badSamplesNum,Swabname_position = Swabname_position)

  if(length(badSamplesNum)){
    badSamplesSwabnames <- x$excludedSamples$df_badSamples$SwabName
    pdf(getPath(x,"DensityPlots_QC", "output/intersampleQC",".pdf"), onefile = TRUE)

    sampGroups <- rep('Valid',ncol(x$bset_raw))
    sampGroups[badSamplesNum] <- "Invalid"
    densityPlot(as.matrix(x$bset_raw), xlim=c(-0.1,1.1), ylim = c(0,7), main="All Samples", sampGroups = sampGroups, pal=c("firebrick3", "forestgreen"))

    x$rgset <- x$rgset[, -badSamplesNum]
    x$bset_raw <- x$bset_raw[, -badSamplesNum]
    x$params <- list(interSampleQC=TRUE)
    densityPlot(as.matrix(x$bset_raw), main="Valid Samples", pal="darkslateblue", legend = FALSE)

    dev.off()
  }else{
    pdf(getPath(x,"DensityPlot_All", "output/intersampleQC",".pdf"), onefile = TRUE)
    densityPlot(as.matrix(x$bset_raw), main="Density plot of All Samples", pal="darkslateblue", legend = FALSE)
    dev.off()
  }
  invisible(x)
}
