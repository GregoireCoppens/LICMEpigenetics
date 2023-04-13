#' @title ProbeBackgroundCheck function
#' @description A function that detects which probes from the rgset do not differ from the background probes.
#'
#' @param x Epigenetic set object.
#'
#' @return an updated Epigenetic set object with the background check included.
#'
#' @import dplyr
#'
#' @importFrom minfi detectionP
#' @importFrom rlang .data
#'
#' @export
#'
ProbeBackgroundCheck <- function(x){
  detP <- detectionP(x$rgset, type = "m+u")
  failed_01<-detP > 0.01

  df_rs <- tibble(names = rownames(failed_01), n_nonsig = rowSums(1*failed_01), n_cols = ncol(detP), perc = .data$n_nonsig/.data$n_cols*100) %>% arrange(.data$perc)
  DetectionP_failed_s50 <- filter(df_rs, .data$n_nonsig>ncol(x$rgset)/2)

  x$excludedProbes$FailedBackgroundProbes <- DetectionP_failed_s50
  invisible(x)
}

# ERROR with openblast, solution: https://support.bioconductor.org/p/122925/
#' @title ProbeNormalisation function
#' @description Functional or Quantile Normalisation of the CpG probes.
#'
#' @param x Epigenetic set object.
#' @param method Type of normalisation algorithm. Options: funnorm, quantile, raw.
#' @param nPCs Number of principal components that funnorm uses.
#' @param save If TRUE, the methylset will be saved to the "minfi_sets" sub-directory.
#' @param verbose When set to 1 or higher more messages will be shown.
#'
#' @return Updated Epigenetic set object.
#'
#' @importFrom minfi preprocessFunnorm
#' @importFrom minfi preprocessQuantile
#'
#' @export
#'
ProbeNormalisation <- function(x, method="funnorm", nPCs = 3, save=FALSE, verbose=0){
  if(ncol(x$rgset)>=500) message("Warning: This set might be too large to fit in your RAM memory")
  if(method=="funnorm"){
    x$mset <- preprocessFunnorm(x$rgset, ratioConvert = FALSE, nPCs = nPCs, verbose = verbose)
    x$params$funnorm <- TRUE
  }
  if(method=="quantile"){
    x$mset <- preprocessQuantile(x$rgset, verbose = verbose)
    x$params$quant <- TRUE
  }
  if(method=="raw"){
    x$mset <- preprocessRaw(x$rgset)
    x$params$preprocessRaw <- TRUE
  }
  x$normalisationMethod <- method
  if(save){
    mset <- x$mset
    save(mset, file=getPath(x, "mset"))
  }
  invisible(gc())
  invisible(x)
}


#' @title ProbeExclusion function
#' @description Removes probes on SNPS, Sex chromosomes or that do not differ from the background noise.
#'
#' @param x Epigenetic set object.
#' @param snps Which SNPs position should considered when excluding probes on SNPs.
#' @param sex If TRUE, probes on the sex chromosomes are excluded.
#' @param background If TRUE, a background check will be completed and those probes which signals do not differ from the background noise, will be removed.
#' @param ExportExcludedProbes if TRUE, the locations of all the removed probes will be saved in the Epigenetic set and exported to a '.csv' file
#' @param delay If TRUE, the excluded probes will be saved, but not yet excluded.
#'
#' @return Updated Epigenetic set object.
#'
#' @importFrom minfi getAnnotation
#' @importFrom minfi addSnpInfo
#' @importFrom minfi dropLociWithSnps
#' @importFrom utils write.csv
#' @importFrom stringr str_detect
#' @export
#'
ProbeExclusion <- function(x, snps = c("CpG", "SBE"), sex = TRUE, background = TRUE, ExportExcludedProbes = FALSE, delay = FALSE){
  # Removing CpGs
  # Source: http://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0194938.s001
  if(delay) x$mset_old <- x$mset # Temp save original mset

  message("CpGs before removal: ",nrow(x$mset))

  # Remove probes that do not exceed background
  if(background){
    x <- ProbeBackgroundCheck(x)
    badCpGs <- x$excludedProbes$FailedBackgroundProbes$names

    if(ExportExcludedProbes){
      x$excludedProbes$backgroundCheck <- data.frame("CpG"=badCpGs, "Type"="FailedBackgroundCheck")
      write.csv(x$excludedProbes$backgroundCheck, getPath(x,"ExcludedProbes_BackgroundCheck", "output/ProbeSelection", ".csv"))
    }

    # goodCpGs <- rownames(x$mset)[-grep(paste0(badCpGs, collapse = "|"), rownames(x$mset))]
    goodCpGs <- rownames(x$mset)[!stringr::str_detect(rownames(x$mset), paste0(badCpGs, collapse = "|"))]
    x$mset <- x$mset[goodCpGs,]

    message("Number of CpGs that didn't pass the background check: ",length(badCpGs))
    message("CpGs after background check: ",nrow(x$mset))
  }

  # Remove Sex
  annotation <- getAnnotation(x$rgset)
  sex_cpgs <- c()
  SNP_cpgs <- c()
  if(sex){
    autosomes <-  annotation[!annotation$chr %in% c("chrX","chrY"), ]
    sex_cpgs <-  rownames(annotation[annotation$chr %in% c("chrX","chrY"), ])
    x$mset <- x$mset[rownames(x$mset) %in% row.names(autosomes),]
    x$params$sexExcl = TRUE
    message("CpGs after sex removal: ",nrow(x$mset))
  } else {
    x$params$sexExcl = FALSE
    message("No CpG sites on Sex Chromosomes removed.")
  }
  if(!is.null(snps)){
    # More info on SNPs: https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html#snps
    # Options snps: Probe, CpG, SBE
    # Remove SNP
    x$mset <- addSnpInfo(x$mset)
    mset_temp <- dropLociWithSnps(x$mset, snps=snps, maf=0)
    SNP_cpgs <- rownames(x$mset[!rownames(x$mset) %in% rownames(mset_temp)])

    x$mset <- mset_temp
    x$params$SnpExcl = TRUE
    message("CpGs after SNP removal: ",nrow(x$mset))
  } else {
    x$params$SnpExcl = FALSE
    message("No CpG sites on SNPs removed.")
  }

  if(x$params$SnpExcl & x$params$sexExcl){
    x$params$sexExcl <- NULL
    x$params$SnpExcl <- NULL
    x$params$SnpSexExcl <- TRUE
  }

  if(delay){ # Reset all params
    x$params$sexExcl <- NULL
    x$params$SnpExcl <- NULL
    x$params$SnpSexExcl <- NULL
    x$mset <- x$mset_old
    x$mset_old <- NULL
  }

  if(ExportExcludedProbes){
    if(sex) df_sex <- data.frame("CpG" = sex_cpgs, "Type"="Sex")
    if(!is.null(snps)) df_SNP <- data.frame("CpG" = SNP_cpgs, "Type"="SNP")

    if(sex & !is.null(snps)){
      x$excludedProbes$SexSNPs <- rbind(df_sex, df_SNP)
    } else{
      if(sex) x$excludedProbes$SexSNPs <- df_sex
      if(!is.null(snps)) x$excludedProbes$SexSNPs <- df_SNP
      if(!sex & is.null(snps)) x$excludedProbes$SexSNPs <- data.frame("CpG"=NULL, "Type"=NULL)
    }

    write.csv(x$excludedProbes$SexSNPs, getPath(x,"ExcludedProbes_SexSNPs", "output/ProbeSelection", ".csv"))
  }


  invisible(gc())
  invisible(x)
}

#' @title ConvertSet function
#' @description Converts the methylset into a beta or M set.
#'
#' @param x Epigenetic set object.
#' @param beta If TRUE, the methylset will be converted in a beta set.
#' @param M If TRUE, the methylset will be converted in a M set.
#' @param save If TRUE, the beta or M set will be saved to the "minfi_sets" sub-directory
#'
#' @return Updated Epigenetic set object.
#'
#' @importFrom minfi getBeta
#' @importFrom minfi getM
#' @importFrom utils write.csv
#' @export
#'
ConvertSet <- function(x, beta=TRUE, M=FALSE, save=FALSE){
  # betas
  if(beta){
    bset_temp <- getBeta(x$mset)
    bset <- bset_temp[!rowSums(!is.finite(bset_temp)),]
    if(!is.null(x$bset_raw)) x$bset_raw <- NULL
    x$bset <- bset
    if(save) save(bset, file=getPath(x, "bset"))

    # bset_na <- bset_temp[!(rownames(bset_temp) %in% rownames(bset)), as.logical(colSums(is.na(bset_temp) > 0))]
    # print(table(!(rownames(bset_temp) %in% rownames(bset))))
    bset_na <- bset_temp[!(rownames(bset_temp) %in% rownames(bset)),as.logical(colSums(is.na(bset_temp) > 0))]
    message("Number of CpGs that didn't pass the background check after conversion: ", nrow(bset_na))
    message("CpGs after Conversion: ", nrow(x$bset))
    if(!is.matrix(bset_na)){
      x$excludedProbes$AfterConversion <- data.frame("CpG"=rownames(bset_temp)[!(rownames(bset_temp) %in% rownames(bset))], "Type"="FailedBackgroundCheck_AfterConversion")
    } else {
      if(nrow(bset_na)>0){
        x$excludedProbes$AfterConversion <- data.frame("CpG"=rownames(bset_na), "Type"="FailedBackgroundCheck_AfterConversion")
      }else{
        x$excludedProbes$AfterConversion <- data.frame(row.names = c("CpG", "Type"))
      }
    }
    if(save) write.csv(x$excludedProbes$AfterConversion, getPath(x,"ExcludedProbes_Conversion", "output/ProbeSelection", ".csv"))
    rm(bset)
  }
  # M-values
  if(M){
    Mset <- getM(x$mset)
    x$Mset <- Mset
    if(save) save(Mset, file=getPath(x, "Mset"))
    rm(Mset)
  }
  invisible(gc())
  invisible(x)
}
