# packages
library(sva)

# Formula: One-sided formula
#' Title
#'
#' @param x
#' @param batch
#' @param RiskFactors
#' @param save
#'
#' @return
#' @export
#'
#' @examples
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
