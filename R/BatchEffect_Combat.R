#' @title BatchEffect_Combat function
#' @description Function to do a batch correction using combat on the beta set from the epigenetic set.
#'
#' @param x Epigenetic set object.
#' @param batch Non-biological factor for which the beta set needs to be corrected.
#' @param RiskFactors Risk factors for which the Combat algorithm needs to correct
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

BatchEffect_Combat <- function(x, batch, RiskFactors, save=TRUE){
  if(class(RiskFactors) != "formula"){
    formula <- as.formula(paste0("~", paste(RiskFactors, collapse = "+")))
  } else {
    formula <- RiskFactors
  }
  x$metadata$df[,batch] <- as.factor(x$metadata$df[,batch])

  if(length(formula)>2) stop("Formula needs to be one-sided")
  mod <-  model.matrix(formula, data=x$metadata$df)
  bset <- ComBat(x$bset, x$metadata$df[,batch], mod)

  x$params[paste0("Combat_", batch)] <- TRUE
  x$bset <- bset
  if(save) save(bset, file=getPath(x, "bset"))
  invisible(x)
}
