#' Title
#'
#' @param x
#' @param df
#' @param by
#'
#' @return
#' @export
#'
#' @examples
MergeMetadata <- function(x, df, by){
  x$metadata$df <- left_join(x$prc$prc_df, df, by=by)
  invisible(x)
}


#' Title
#'
#' @param x
#' @param batch
#' @param n
#' @param xlim
#' @param ylim
#'
#' @return
#' @export
#'
#' @examples
BatchEffect_PointPlot <- function(x, batch, n = 6, xlim = c(-1000, 2500), ylim = c(-1000, 1000)){
  if(n>6) warning("[Warning] Too many groups per plot. It is advised to not exceed 6 groups per plot (n=6).")

  x$metadata$df[,batch] <- as.factor(x$metadata$df[,batch])
  path <- getPath(x, paste0("BatchEffect_",batch,"_Pointplot"), "output/BatchEffect", ".pdf")
  pdf(path, onefile = TRUE)

  # Save the path to the output in the object
  if(is.null(x$prc$BatchEffect_pointPlot_paths)) x$prc$BatchEffect_pointPlot_paths <- c()
  x$prc$BatchEffect_pointPlot_paths <- unique(c(x$prc$BatchEffect_pointPlot_paths, path))

  # Print all in one plot
  print(ggplot(data = x$metadata$df, aes(x=PC1, y=PC2)) + geom_point() + coord_cartesian(xlim = xlim, ylim = ylim)) # Show all

  # Group per 6 chips
  batch_unique <- unique(x$metadata$df[,batch]) %>% sample()

  if(length(batch_unique)<n) n <- length(batch_unique)

  for(i in 1:floor(length(batch_unique)/n)){
    df_filter <- x$metadata$df %>% filter(get(batch) %in% batch_unique[(i*n-n+1):(i*n)])
    g <- ggplot(data = df_filter, aes(x=PC1, y=PC2)) + geom_point(aes_string(shape=batch)) + coord_cartesian(xlim = xlim, ylim = ylim)
    print(g)
    g2 <- ggplot(data = df_filter, aes(x=PC1, y=PC3)) + geom_point(aes_string(shape=batch)) + coord_cartesian(xlim = xlim, ylim = ylim)
    print(g2)
    g3 <- ggplot(data = df_filter, aes(x=PC2, y=PC3)) + geom_point(aes_string(shape=batch)) + coord_cartesian(xlim = xlim, ylim = ylim)
    print(g3)
  }

  # Remainder
  df_filter <- x$metadata$df %>% filter(get(batch) %in% batch_unique[c((length(batch_unique)-length(batch_unique)%%n)+1:length(batch_unique))])
  if(nrow(df_filter)>0){
    print(ggplot(data = df_filter, aes(x=PC1, y=PC2)) + geom_point(aes_string(shape=batch)) + coord_cartesian(xlim = xlim, ylim = ylim))
    print(ggplot(data = df_filter, aes(x=PC1, y=PC3)) + geom_point(aes_string(shape=batch)) + coord_cartesian(xlim = xlim, ylim = ylim))
    print(ggplot(data = df_filter, aes(x=PC2, y=PC3)) + geom_point(aes_string(shape=batch)) + coord_cartesian(xlim = xlim, ylim = ylim))
  }
  dev.off()
  invisible(x)
}


#' Title
#'
#' @param x
#' @param Cat_vars
#' @param Cont_vars
#'
#' @return
#' @export
#'
#' @examples
BatchEffect_PCAHeatmap <- function(x, Cat_vars = c(), Cont_vars=c()){
  if(is.null(Cat_vars) & is.null(Cont_vars)) stop("No categorical or continuous variables provided. Please fill in at least one categorical or continuous variable.")
  PCA_count <- x$prc$n_pcas
  variance_percent <- x$prc$pca_eig$variance.percent[1:PCA_count] %>% round(2)
  PCA_names <- paste0("PC",1:PCA_count)

  # Categorical values
  df_heatmap_cat <- lapply(Cat_vars, function(i){
    x$metadata$df %>%
      select(starts_with("PC"), i) %>%
      {lapply(PCA_names, function(j){
        summarise(., !!i := summary(aov(get(j)~get(i), data=.))[[1]][["Pr(>F)"]][[1]])
      }) %>% bind_rows()
      }
  }) %>% bind_cols()

  # Continues variables
  df_heatmap_cont <- lapply(Cont_vars, function(i){
    x$metadata$df %>%
      select(starts_with("PC"), i) %>%
      {lapply(PCA_names, function(j){
        summarise(., !!i := cor.test(x = get(j), y = get(i), data=.)$p.value)
      }) %>% bind_rows()
      }
  }) %>% bind_cols()

  x$heatmap$df_heatmap <- bind_cols(df_heatmap_cat, df_heatmap_cont) %>%
    mutate(PC = paste0(PCA_names," (",variance_percent,"%)")) %>%
    gather("vars","pvalue", -PC) %>%
    mutate(pgroup=factor(
      findInterval(pvalue, c(0, 1e-10, 1e-5, 0.001, 0.05)),
      levels=c(1,2,3,4,5),
      labels = c("<1e-10", "<1e-5", "<0.001", "<0.05", "<1"))) %>%
    arrange(vars, PC)

  pdf(getPath(x, "PCAHeatmap", "output/BatchEffect", ".pdf"), onefile = TRUE)
  x$heatmap$hm <- ggplot(data= x$heatmap$df_heatmap, aes(x = PC, y = vars, fill=pgroup)) +
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

#' Title
#'
#' @param x
#' @param df
#' @param by
#' @param batch
#' @param Cat_vars
#' @param Cont_vars
#' @param n
#' @param xlim
#' @param ylim
#'
#' @return
#' @export
#'
#' @examples
BatchEffect_Visualise <- function(x, df, by, batch, Cat_vars = c(), Cont_vars=c(), n = 6, xlim = c(-1000, 2500), ylim = c(-1000, 1000)){
  MergeMetadata(x, df, by) %>%
    BatchEffect_PointPlot(batch, n, xlim, ylim) %>%
    BatchEffect_PCAHeatmap(Cat_vars, Cont_vars) %>%
    invisible()
}
