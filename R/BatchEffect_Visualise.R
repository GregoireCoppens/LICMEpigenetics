#' @title  MergeMetadata function
#' @description Function to merge the metadata with the principal component data.
#'
#' @param x Epigenetic set object.
#' @param df Metadata
#' @param by Name of the column to merge existing data with
#'
#' @return Updated Epigenetic set object.
#'
#' @importFrom dplyr left_join
#' @export
#'
MergeMetadata <- function(x, df, by){
  x$metadata$df <- left_join(x$prc$prc_df, df, by=by)
  invisible(x)
}


#' @title BatchEffect_PointPlot Function
#' @description Make pointplots from the first 3 principical components from the beta set to assess wether batch effect is present. Every Plot should only contain maximum 6 groups to maintain an overview.
#'
#' @param x Epigenetic set object.
#' @param batch The non-biological factor which we want to investigate for batch effect.
#' @param n Number of samples per plot. It is advised to not exceed 6 groups per plot.
#' @param xlim range of the x-axis
#' @param ylim range of the y-axis
#'
#' @return Updated Epigenetic set object.
#' @importFrom rlang .data
#' @export
#'
#'
BatchEffect_PointPlot <- function(x, batch, n = 6, xlim = c(-1000, 2500), ylim = c(-1000, 1000)){
  if(n>6) warning("[Warning] Too many groups per plot. It is advised to not exceed 6 groups per plot (n=6).")

  x$metadata$df[,batch] <- as.factor(x$metadata$df[,batch])
  path <- getPath(x, paste0("BatchEffect_",batch,"_Pointplot"), "output/BatchEffect", ".pdf")
  pdf(path, onefile = TRUE)

  # Save the path to the output in the object
  if(is.null(x$prc$BatchEffect_pointPlot_paths)) x$prc$BatchEffect_pointPlot_paths <- c()
  x$prc$BatchEffect_pointPlot_paths <- unique(c(x$prc$BatchEffect_pointPlot_paths, path))

  # Print all in one plot
  print(ggplot(data = x$metadata$df, aes(x=.data$PC1, y=.data$PC2)) + geom_point() + coord_cartesian(xlim = xlim, ylim = ylim)) # Show all

  # Group per 6 chips
  batch_unique <- unique(x$metadata$df[,batch]) %>% sample()

  if(length(batch_unique)<n) n <- length(batch_unique)

  for(i in 1:floor(length(batch_unique)/n)){
    df_filter <- x$metadata$df %>% filter(get(batch) %in% batch_unique[(i*n-n+1):(i*n)])
    g <- ggplot(data = df_filter, aes(x=.data$PC1, y=.data$PC2)) + geom_point(aes_string(shape=batch)) + coord_cartesian(xlim = xlim, ylim = ylim)
    print(g)
    g2 <- ggplot(data = df_filter, aes(x=.data$PC1, y=.data$PC3)) + geom_point(aes_string(shape=batch)) + coord_cartesian(xlim = xlim, ylim = ylim)
    print(g2)
    g3 <- ggplot(data = df_filter, aes(x=.data$PC2, y=.data$PC3)) + geom_point(aes_string(shape=batch)) + coord_cartesian(xlim = xlim, ylim = ylim)
    print(g3)
  }

  # Remainder
  df_filter <- x$metadata$df %>% filter(get(batch) %in% batch_unique[c((length(batch_unique)-length(batch_unique)%%n)+1:length(batch_unique))])
  if(nrow(df_filter)>0){
    print(ggplot(data = df_filter, aes(x=.data$PC1, y=.data$PC2)) + geom_point(aes_string(shape=batch)) + coord_cartesian(xlim = xlim, ylim = ylim))
    print(ggplot(data = df_filter, aes(x=.data$PC1, y=.data$PC3)) + geom_point(aes_string(shape=batch)) + coord_cartesian(xlim = xlim, ylim = ylim))
    print(ggplot(data = df_filter, aes(x=.data$PC2, y=.data$PC3)) + geom_point(aes_string(shape=batch)) + coord_cartesian(xlim = xlim, ylim = ylim))
  }
  dev.off()
  invisible(x)
}


#' @title BatchEffect_PCAHeatmap function
#' @description Function that creates a p-value heatmap from categorical and continuous variables using F-tests and correlation tests (Pearson's product) to evaluate the presence of batch effect.
#'
#' @param x Epigenetic set object.
#' @param Cat_vars Categorical variables for the heatmap
#' @param Cont_vars Continues variables for the heatmap
#'
#' @return Updated Epigenetic set object.
#'
#' @import dplyr
#' @import ggplot2
#'
#' @importFrom rlang :=
#' @importFrom stats aov
#'
#' @importFrom rlang .data
#' @importFrom tidyr gather
#'
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#'
#' @export
#'

BatchEffect_PCAHeatmap <- function(x, Cat_vars = c(), Cont_vars=c()){
  if(is.null(Cat_vars) & is.null(Cont_vars)) stop("No categorical or continuous variables provided. Please fill in at least one categorical or continuous variable.")
  PC_count <- x$prc$n_pcs
  variance_percent <- x$prc$pca_eig$variance.percent[1:PC_count] %>% round(2)
  PCA_names <- paste0("PC",1:PC_count)

  # Categorical values
  df_heatmap_cat <- lapply(Cat_vars, function(i){
    df_temp <- x$metadata$df %>% select(starts_with("PC"), i)
    lapply(PCA_names, function(j){
      summarise(df_temp, !!i := summary(aov(get(j)~get(i), data=df_temp))[[1]][["Pr(>F)"]][[1]])
    }) %>% bind_rows()
  }) %>% bind_cols()

  # Continues variables
  df_heatmap_cont <- lapply(Cont_vars, function(i){
    df_temp <- x$metadata$df %>% select(starts_with("PC"), i)
    lapply(PCA_names, function(j){
      summarise(df_temp, !!i := cor.test(x = get(j), y = get(i), data=df_temp)$p.value)
    }) %>% bind_rows()
  }) %>% bind_cols()

  x$heatmap$df_heatmap <- bind_cols(df_heatmap_cat, df_heatmap_cont) %>%
    mutate(PC = paste0(PCA_names," (",variance_percent,"%)")) %>%
    gather("vars","pvalue", -.data$PC) %>%
    mutate(pgroup=factor(findInterval(.data$pvalue, c(0, 1e-10, 1e-5, 0.001, 0.05)),
      levels=c(1,2,3,4,5),
      labels = c("<1e-10", "<1e-5", "<0.001", "<0.05", "<1"))) %>%
    arrange(.data$vars, .data$PC)

  pdf(getPath(x, "PCAHeatmap", "output/BatchEffect", ".pdf"), onefile = TRUE)
  x$heatmap$hm <- ggplot(data= x$heatmap$df_heatmap, aes(x = .data$PC, y = .data$vars, fill=.data$pgroup)) +
    geom_tile(color="white", size=0.5) +
    scale_fill_manual(name="P-values",
                      values = c("<1e-10"="dark red", "<1e-5"="red", "<0.001"="orange", "<0.05"="pink", "<1"="white")) +
    theme(axis.ticks=element_blank()) +
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
    coord_equal()+
    ggtitle("P-Value heatmap") + xlab("Principal Components (Variance Explained)") + ylab("")
  print(x$heatmap$hm)
  dev.off()
  invisible(x)
}

#' @title BatchEffect_Visualise function
#' @description Function to visualise batch effect using point plots and heatmaps. This function combines the MergeMetadata, BatchEffect_PointPlot and BatchEffect_PCAHeatmap function.
#'
#' @param x Epigenetic set object.
#' @param df Metadata
#' @param by Name of the column to merge existing data with
#' @param batch The non-biological factor that we want to investigate for batch effect.
#' @param Cat_vars Categorical variables for the heatmap
#' @param Cont_vars Continues variables for the heatmap
#' @param n number of samples per plot
#' @param xlim range of the x-axis
#' @param ylim range of the y-axis
#'
#' @return Updated Epigenetic set object.
#'
#' @import dplyr
#' @import ggplot2
#'
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#'
#' @importFrom rlang :=
#'
#' @importFrom stats aov
#' @importFrom stats cor.test
#'
#' @export
#'
BatchEffect_Visualise <- function(x, df, by, batch, Cat_vars = c(), Cont_vars=c(), n = 6, xlim = c(-1000, 2500), ylim = c(-1000, 1000)){
  MergeMetadata(x, df, by) %>%
    BatchEffect_PointPlot(batch, n, xlim, ylim) %>%
    BatchEffect_PCAHeatmap(Cat_vars, Cont_vars) %>%
    invisible()
}
