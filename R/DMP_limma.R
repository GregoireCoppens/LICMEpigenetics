#' Detect differentially methylated positions in the DNA using the limma package.
#'
#' This function first uses the `lmfit` function from the `limma` package with methylset and the generated design matrix as its inputs.
#' Next it uses the `contrast.fit` function from the `limma` package with the `lmfit` output and the generated contrast matrix as its inputs.
#' The contrast matrix specifies the comparison of interest (e.g.: Patients vs Controls).
#' To compute a moderated t-statistic of differential expression, it uses empirical Bayes moderation of the standard error.
#'
#' @param set Preprocessed methylation set, can be genomic methyl set or genomic ratio set.
#' @param Patient_info Data frame with adjustment variables.
#' @param cat_vars The variable names in Patient_info that are categorical variables
#' @param cont_vars The variable names in Patient_info that are continuous variables
#' @param Group The variable in Patient_info that represents which observations belong to which groups
#' @param confint To enable confidence interval for the logFC, set to TRUE or change to CI range.
#' 
#' @return A list with limma output variables and temporary results
#' \itemize{
#'    \item DMPresult - Data table with statistics per CpG sites
#'    \item model - eBayes fitted linear model
#'    \item design - Design matrix
#'    \item contMatrix - Contrast matrix
#'    \item Sign_CpG - CpG names of DMPs
#' }
#'
#' @export
#'
#'
#' @importFrom limma makeContrasts
#' @importFrom limma contrasts.fit
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @importFrom minfi getM
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter
#'
DMP_limma <- function(set, Patient_info, cat_vars = NULL, cont_vars=NULL, Group=NULL, confint=TRUE){
  # Define inputs
  if(class(set) =="GenomicRatioSet") GRcset <- set
  if(class(set) =="GenomicMethylSet") GRcset <- minfi::ratioConvert(set)

  if(is.null(cat_vars) & !is.null(Group)) cat_vars <- c(Group)
  if(is.null(Group)) Group <- cat_vars[1]
  if(Group %in% cat_vars) cat_vars <- cat_vars[!cat_vars %in% Group]

  # Prepare Design
  if(!is.null(cont_vars)){
    if(length(cat_vars)>0){
      design <- model.matrix(
        as.formula(
          paste("~0",
                paste("factor(",Group,")", collapse = "+"),
                paste(
                  "factor(", cat_vars, ", levels=sort(unique(",cat_vars ,"), decreasing = T))",
                  sep = "",
                  collapse = " + "),
                paste(cont_vars, collapse = "+"), sep = "+" )),
        Patient_info)
    } else{
      design <- model.matrix(
        as.formula(
          paste("~0",
                paste("factor(",Group,")", collapse = "+"),
                paste(cont_vars, collapse = "+"), sep = "+" )),
        Patient_info)
    }
  } else {
    if(length(cat_vars)>0){
      design <- model.matrix(
        as.formula(
          paste("~0",
                paste("factor(",Group,")", collapse = "+"),
                paste(
                  "factor(", cat_vars, ", levels=sort(unique(",cat_vars ,"), decreasing = T))",
                  sep = "",
                  collapse = " + "), sep = "+" )),
        Patient_info)
    } else {
      design <- model.matrix(
        as.formula(
          paste("~0",
                paste("factor(",Group,")", collapse = "+"),
                sep = "+" )),
        Patient_info)
    }
  }

  colnames(design) <-   gsub("\\)", "_", gsub("\\,.+\\)", "_", gsub(pattern = "factor\\(", replacement = "", x = colnames(design))))

  ContrastFormula <- paste(colnames(design)[c(2,1)], collapse = "-") %>% as.vector()
  contMatrix <- limma::makeContrasts(contrasts = ContrastFormula, levels=design)

  # Calculate
  fit <- limma::lmFit(minfi::getM(GRcset), design)
  fit2 <- limma::contrasts.fit(fit, contMatrix)
  fit_EB <- limma::eBayes(fit2)

  DMPresult <- limma::topTable(fit_EB,
                      n=Inf,
                      coef=1,
                      confint=confint)

  Sign_CpG_df <- DMPresult %>%
    tibble::rownames_to_column("CpG_name") %>%
    dplyr::filter(.data$adj.P.Val<=0.05)
  Sign_CpG <- Sign_CpG_df$CpG_name

  DMPresult_sig <-DMPresult %>% dplyr::filter(.data$adj.P.Val<=0.05)

  return(list(
    DMPresult=DMPresult,
    model = fit_EB,
    design = design, contMatrix = contMatrix,
    Sign_CpG = Sign_CpG))
}
