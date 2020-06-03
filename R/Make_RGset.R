# Read idats and make RGset

# install.packages("BiocManager")
# BiocManager::install("minfi")
# BiocManager::install("IlluminaHumanMethylationEPICmanifest")
# BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library(minfi)
library(dplyr)


#' @title Make_RGset function
#' @description This function builds a red-green set (rgset) from the given idat files. The filenames from the idat files will be as unique identifiers.
#'
#' @param idat_paths Path to the idat files
#' @param name Project name (this name will be in all the following filenames as a unique identifier for the project).
#' @param subselect Select only the idat files that have a given string in their filenames.
#' @param test Make a smaller set that can be used to test all the functions.
#' @param n_testsamples The number of samples in a test set (if test = TRUE). Default: 5.
#' @param save Set to TRUE if you want to save the rgset as an .Rdata file.
#'
#' @return Epigenetic set object.
#' @export
#'

Make_RGset <- function(idat_paths="./", name="", subselect="", test=FALSE, n_testsamples=5, save=FALSE){
  if(name==""){
    name <- lapply(idat_paths, function(path){path %>% base::strsplit("/") %>% unlist() %>% tail(n=1)}) %>% paste(collapse = "_")
    if(name=="") name <- "Methylome"
  }
  name <- paste0(name, subselect)

  # Get file names
  pathNames <- lapply(idat_paths, function(idat_path){
    fileNames <- list.files(idat_path, pattern = "\\.idat") # Read all idat files
    baseNames <- fileNames[grepl("_Grn.", fileNames)] %>% # Get all filenames that contain "_Grn"
      {.[grepl(subselect,.)]} %>% # Only get those filenames that contain subselection
      {lapply(., function(i){strsplit(i, '_Grn')[[1]][[1]]})} %>% #Get the name of the file without "_Grn"
      unlist()
    paste0(idat_path,"/",baseNames)
  }) %>% unlist()

  #Import files
  ## Force is set to TRUE because we are working with EPIC data.
  rgset <- read.metharray(pathNames, force=TRUE, verbose = TRUE)

  if(test){
    # randomsamples <- sample(ncol(rgset), floor(ncol(rgset)/20))
    randomsamples <- sample(ncol(rgset), n_testsamples)
    rgset_test <- rgset[,randomsamples]
    if(save){
      save_path <- getPath(name, subName = "rgset_test", "minfi_sets")
      saveRDS(rgset_test, save_path)
    }
    return(list(rgset = rgset_test, name = name))
  }
  # Save rgset
  if(save){
    save_path <- getPath(name, subName = "rgset", "minfi_sets")
    saveRDS(rgset, save_path)
  }
  return(list(rgset = rgset, name = name))
}
