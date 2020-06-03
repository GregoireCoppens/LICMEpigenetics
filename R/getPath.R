#' @title GetPath function
#' @description Function to automatically make a unique and descriptive path name based on the properties of the epigenetic set.
#' @param x Epigenetic set object.
#' @param subName Name of the file.
#' @param subDir Directory of the file. The default is the "minfi_sets" directory.
#' @param extension Extension of the file. Default is ".Rdata".
#' @param mainDir Main directory to save the file to. The default is the directory upstream from the working directory.
#'
#' @return Returns a filepath.
#' @export
#'
#' @examples
#' getPath(list(name="test", params=list(combatChip = TRUE, combatRow=TRUE, recalibrate=FALSE)))
#' getPath("rgset", "plot3", "output")
#'
getPath <- function(x, subName="", subDir = "minfi_sets", extension=".Rdata",  mainDir=".."){
  if(is.list(x)){
    main_name <- x$name
    params <- x$params
  }
  if(is.character(x)){
    main_name <- x
    params <- list()
  }

  # Create subdir if it doesn't exist
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

  # Create path
  ## Create substring using all the params that are TRUE to make filename unique
  params_str <- paste0(names(unlist(lapply(params, function(i){if(i)i}))),collapse="_")

  # Remove any possible extensions in the name
  main_name <- base::strsplit(main_name, "\\.")[[1]][1]

  # Remove any spaces, forward slashes or points from the subname
  subName <- gsub(" |/|\\.", "_", subName)

  # Remove any strings without contents
  filename_vector <- c(main_name, params_str, subName)[base::which(c(main_name, params_str, subName) != "")]

  return(file.path(mainDir, subDir, paste0(paste(filename_vector, collapse = "_"), extension)))
}
# getPath(list(name="test", params=list(combatChip = TRUE, combatRow=TRUE, recalibrate=FALSE)))
# getPath("rgset", "plot3", "output")
