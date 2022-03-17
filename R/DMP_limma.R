#' Detect differentially methylated positions in the DNA using the limma package.
#'
#' This function first uses the `lmfit` function from the `limma` package with methylset and the generated design matrix as its inputs.
#' Next it uses the `contrast.fit` function from the `limma` package with the `lmfit` output and the generated contrast matrix as its inputs.
#' The contrast matrix specifies the comparison of interest (e.g.: Patients vs Controls).
#' To compute a moderated t-statistic of differential expression, it uses empirical Bayes moderation of the standard error.
#'
#' @param GRcset Preprocessed Green-Red methylation set
#' @param Patient_info Data frame with adjustment variables.
#' @param cat_vars The variable names in Patient_info that are categorical variables
#' @param cont_vars The variable names in Patient_info that are continuous variables
#' @param Group The variable in Patient_info that represents which observations belong to which groups
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
#' @examples
#'
DMP_limma <- function(GRcset, Patient_info, cat_vars, cont_vars=c(), Group=NULL){
  # Define inputs
  if(is.null(Group)) Group <- cat_vars[1]
  if(Group %in% cat_vars) cat_vars <- cat_vars[!cat_vars %in% Group]

  # Prepare Design
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
  colnames(design) <- colnames(design) %>% gsub(pattern = "factor\\(", replacement = "", x = .) %>%  gsub("\\,.+\\)", "_", .) %>% gsub("\\)", "_", .)

  ContrastFormula <- paste(colnames(design)[c(2,1)], collapse = "-") %>% as.vector()
  contMatrix <- makeContrasts(contrasts = ContrastFormula, levels=design)

  # Calculate
  fit <- lmFit(getM(GRcset), design)
  fit2 <- contrasts.fit(fit, contMatrix)
  fit_EB <- eBayes(fit2)

  DMPresult <- topTable(fit_EB,
                      n=Inf,
                      coef=1)

  Sign_CpG <- DMPresult %>%
    rownames_to_column("CpG_name") %>%
    filter(adj.P.Val<=0.05) %>%
    .$CpG_name

  DMPresult_sig <-DMPresult %>% filter(adj.P.Val<=0.05)

  return(list(
    DMPresult=DMPresult,
    model = fit_EB,
    design = design, contMatrix = contMatrix,
    Sign_CpG = Sign_CpG))
}
