#' @title BatchEffect_prc function
#' @description Creates principical components from the beta set to assess batch effect.
#'
#' @param x Epigenetic set object.
#' @param n_pcs numbers of principical components that need to be created. The default is 8.
#' @param filename_info A vector that explains the structure of the idat filename, such that the experiment info can be extrated using the idat filename.
#' @param filename_separator A string that identifies how the filename structure is devided. The default is "_"
#'
#' @return Updated Epigenetic set object.
#'
#' @import dplyr
#' @import factoextra
#'
#' @importFrom utils write.csv
#' @importFrom stats prcomp
#' @importFrom tidyr separate
#'
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @export
#'
BatchEffect_prc <- function(x, n_pcs=8, filename_info = c("Chip", "Well", "Sample_nr", "Sample_type", "Ptn_Cntrl", "Patnr", "Early_Late", "Timepoint"), filename_separator="_"){
  if(is.null(x$prc$n_pcs)|n_pcs != 8) x$prc$n_pcs <- n_pcs
  if(filename_separator=="" & !is.null(filename_info)) filename_separator  <-  "__/"
  if(is.null(filename_info)){
    filename_info <-  c("filenames")
    filename_separator <- "__/"
  }
  x$prc$prc_results <- prcomp(t(x$bset), scale = TRUE)
  x$prc$pca_eig <- get_eig(x$prc$prc_results)
  x$prc$pca_ind <- get_pca_ind(x$prc$prc_results)

  PC12 <- x$prc$pca_ind$coord[,1:x$prc$n_pcs] %>% as.data.frame()
  colnames(PC12) <- paste0("PC", 1:x$prc$n_pcs)
  x$prc$prc_df <- mutate(PC12, filenames=rownames(PC12)) %>%
    separate(.data$filenames,filename_info, filename_separator, convert = TRUE, remove = FALSE, fill = "right")

  write.csv(x$prc$prc_df, getPath(x, "prc_df", "output/BatchEffect", ".csv"))

  pdf(getPath(x, "prc_visualisation", "output/BatchEffect", ".pdf"), onefile = T)
  print(fviz_eig(x$prc$prc_results))
  dev.off()

  invisible(x)
}
