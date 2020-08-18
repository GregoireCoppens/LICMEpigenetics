#' @title BatchEffect_Combat function
#' @description Function to do a batch correction using combat on the beta set from the epigenetic set.
#'
#' @param x Epigenetic set object.
#' @param batch Non-biological factor for which the beta set needs to be corrected.
#' @param RiskFactors Risk factors for which the Combat algorithm needs to correct
#' @param M If TRUE, the combat will run using the M values instead of the beta values.
#' @param save If TRUE, the new bset will be saved into the minfi_sets directory.
#'
#' @return Updated Epigenetic set object.
#'
#' @import sva
#' @importFrom stats as.formula
#' @importFrom stats model.matrix
#'
#' @export
#'

BatchEffect_Combat <- function(x, batch, RiskFactors, M=FALSE, save=TRUE){
  if(class(RiskFactors) != "formula"){
    formula <- as.formula(paste0("~", paste(RiskFactors, collapse = "+")))
  } else {
    formula <- RiskFactors
    RiskFactors <- strsplit(as.character(RiskFactors)[[2]], " \\+ ")[[1]]
  }
  df <- x$metadata$df[,c(RiskFactors, batch)]
  df[,batch] <- as.factor(df[,batch])

  if(length(formula)>2) stop("Formula needs to be one-sided")
  # mod <-  model.matrix(formula, data=x$metadata$df)
  mod <-  model.matrix(formula, data=df)
  invisible(gc())
  if(M){
    # Beta values can not contain 0 or 1 value
    ## Zero
    min2 <- min(x$bset[x$bset>0])
    newZero <- 10**floor(log10(min2))
    x$bset[x$bset == 0] <- newZero

    ## One
    max2 <- max(x$bset[x$bset<1])
    newOne <- 1-10**floor(log10(1-max2))
    x$bset[x$bset == 1] <- newOne

    # Clean up memory
    Mset <- log2(x$bset/(1-x$bset))
    x$bset <- NULL
    invisible(gc())
    Mset <- ComBat(Mset, df[,batch], mod)
    x$bset <- 2**Mset/(2**Mset+1)
  }else{
    x$bset <- ComBat(x$bset, df[,batch], mod)
  }
  # Save
  x$params[paste0("Combat", batch)] <- TRUE
  invisible(gc())
  bset <- x$bset
  if(save) save(bset, file=getPath(x, "bset"))
  x$BatchCorrections <- c(x$BatchCorrections, batch)
  invisible(x)
}
