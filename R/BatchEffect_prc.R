#Principal Component Preprocessing

# Packages
# install.packages("factoextra")

library(factoextra)
library(dplyr)
library(tidyr)


#' Title
#'
#' @param x
#' @param n_pcas
#' @param filename_info
#' @param filename_separator
#'
#' @return
#' @export
#'
#' @examples
BatchEffect_prc <- function(x, n_pcas=8, filename_info = c("Chip", "Well", "Swabnr", "Experiment_Method", "Patient_Group", "Patient_Number", "Experiment_Group", "Time"), filename_separator="_"){
  if(is.null(x$prc$n_pcas)|n_pcas != 8) x$prc$n_pcas <- n_pcas
  if(filename_separator=="" & !is.null(filename_info)) filename_separator  <-  "__/"
  if(is.null(filename_info)){
    filename_info <-  c("filenames")
    filename_separator <- "__/"
  }
  x$prc$prc_results <- prcomp(t(x$bset), scale = TRUE)
  x$prc$pca_eig <- get_eig(x$prc$prc_results)
  x$prc$pca_ind <- get_pca_ind(x$prc$prc_results)

  PC12 <- x$prc$pca_ind$coord[,1:x$prc$n_pcas] %>% as.data.frame()
  colnames(PC12) <- paste0("PC", 1:x$prc$n_pcas)
  x$prc$prc_df <- PC12 %>%
    {mutate(., filenames=rownames(.))} %>%
    separate(filenames,filename_info, filename_separator, convert = TRUE, remove = FALSE, fill = "right")

  write.csv(x$prc$prc_df, getPath(x, "prc_df", "output/BatchEffect", ".csv"))

  pdf(getPath(x, "prc_visualisation", "output/BatchEffect", ".pdf"), onefile = T)
  print(fviz_eig(x$prc$prc_results))
  dev.off()

  invisible(x)
}
