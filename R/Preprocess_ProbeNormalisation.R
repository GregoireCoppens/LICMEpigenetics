
library(minfi)
library(minfiDataEPIC)

# install.packages("tidyverse")
library(dplyr)
library(tidyr)
library(readxl)

# BiocManager::install("factoextra", update=FALSE)
library(factoextra)



#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
ProbeBackgroundCheck <- function(x){
  detP <- detectionP(x$rgset, type = "m+u")
  failed_01<-detP > 0.01
  # save(failed_01, file = detectionP_filename)

  df_rs <- tibble(names = rownames(failed_01), n_nonsig = rowSums(1*failed_01), n_cols = ncol(detP), perc = n_nonsig/n_cols*100) %>% arrange(perc)
  DetectionP_failed_s50 <- df_rs %>% filter(n_nonsig>ncol(x$rgset)/2)

  x$excludedProbes$FailedBackgroundProbes <- DetectionP_failed_s50
  invisible(x)
}

# ERROR with openblast, solution: https://support.bioconductor.org/p/122925/
#' Title
#'
#' @param x
#' @param quant
#' @param save
#'
#' @return
#' @export
#'
#' @examples
ProbeNormalisation <- function(x, quant=TRUE, save=FALSE){
  if(ncol(x$rgset)>=500) message("Warning: This set might be too large to fit in your RAM memory")
  if(!quant){
    x$mset <- preprocessFunnorm(x$rgset, ratioConvert = FALSE, verbose = 2)
    x$params$funnorm <- TRUE
  } else {
    x$mset <- preprocessQuantile(x$rgset)
    x$params$quant <- TRUE
  }
  if(save){
    mset <- x$mset
    save(mset, file=getPath(x, "mset"))
  }
  invisible(gc())
  invisible(x)
}


#' Title
#'
#' @param x
#' @param snps
#' @param sex
#' @param background
#' @param ExportSNPs
#'
#' @return
#' @export
#'
#' @examples
ProbeExclusion <- function(x, snps = c("CpG", "SBE"), sex = TRUE, background = TRUE, ExportSNPs = TRUE){
  # Removing CpGs
  # Source: http://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0194938.s001
  # Remove Sex
  message("CpGs before removal: ",nrow(x$mset))
  annotation <- getAnnotation(x$rgset)
  sex_cpgs <- c()
  SNP_cpgs <- c()
  if(sex){
    autosomes <-  annotation[!annotation$chr %in% c("chrX","chrY"), ]
    sex_cpgs <-  annotation[annotation$chr %in% c("chrX","chrY"), ] %>% rownames()
    x$mset <- x$mset[rownames(x$mset) %in% row.names(autosomes),]
    x$params$sexExcl = TRUE
    message("CpGs after sex removal: ",nrow(x$mset))
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
  }
  if(x$params$SnpExcl & x$params$sexExcl){
    x$params$sexExcl <- NULL
    x$params$SnpExcl <- NULL
    x$params$SnpSexExcl <- TRUE
  }
  if(ExportSNPs){
    df_sex <- data.frame("CpG" = sex_cpgs, "Type"="Sex")
    df_SNP <- data.frame("CpG" = SNP_cpgs, "Type"="SNP")
    x$excludedProbes$SexSNPs <- rbind(df_sex, df_SNP)

    write.csv(x$excludedProbes$SexSNPs, getPath(x,"ExcludedProbes_SexSNPs", "output/ProbeSelection", ".csv"))
  }

  if(background){
    x <- ProbeBackgroundCheck(x)
    badCpGs <- x$excludedProbes$FailedBackgroundProbes$names

    x$excludedProbes$backgroundCheck <- data.frame("CpG"=badCpGs, "Type"="FailedBackgroundCheck")

    write.csv(x$excludedProbes$backgroundCheck, getPath(x,"ExcludedProbes_BackgroundCheck", "output/ProbeSelection", ".csv"))

    goodCpGs <- rownames(x$mset)[-grep(paste0(badCpGs, collapse = "|"), rownames(x$mset))]
    x$mset <- x$mset[goodCpGs,]
    # message("Number of CpGs that didn't pass the background check: ",length(badCpGs))
    message("CpGs after background check: ",nrow(x$mset))
  }
  invisible(gc())
  invisible(x)
}

#' Title
#'
#' @param x
#' @param beta
#' @param M
#' @param save
#'
#' @return
#' @export
#'
#' @examples
ConvertSet <- function(x, beta=TRUE, M=TRUE, save=TRUE){
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
