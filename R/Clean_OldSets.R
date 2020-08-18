#' @title Clean_OldSets function
#' @description Removes all sets except the bset to make place available in memory.
#'
#' @param x Epigenetic set object.
#' @param prc If TRUE, the prc element in the set will be removed.
#' @param SaveOld If TRUE, the whole set will be save as "AllSets" before any sets are removed.
#'
#' @return Updated Epigenetic set object.
#' @export
#'
Clean_OldSets <- function(x, prc=TRUE, SaveOld=F){
  if(SaveOld) save(x, file=getPath(x, "AllSets"))
  x$rgset <- NULL
  x$mset <- NULL
  x$Mset <- NULL
  if(prc) x$prc <- NULL
  return(x)
}
