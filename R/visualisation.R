#' @export

draw_pca <- function(data_frame,
                     log = T,
                     group_by = NULL,
                     x_axis = "PC1",
                     y_axis = "PC2",
                     return_data = F) {

  if (log == T) {
    pca_dat <- prcomp(t(log2(data_frame)))
  } else {
    pca_dat <- prcomp(t(data_frame))
  }

  pca_df <- data.frame(pc1 = pca_dat$x[, x_axis],
                       pc2 = pca_dat$x[, y_axis],
                       groups = group_by)

  if (return_data == T) {
    return(pca_df)
  }

  stdev_val <- pca_dat$sdev^2
  stdev_val <- round(stdev_val/sum(stdev_val)*100, 1)

  ggplot(pca_df, aes(pc1, pc2)) +
    geom_point(cex = 3.3, shape = 21, aes(fill = groups), alpha = 0.8, stroke = 0.3) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          aspect.ratio = 1,
          legend.key = element_blank()) +
    labs(x = paste0(x_axis, ": ", stdev_val[as.numeric(str_extract(x_axis, "[0-9]"))], "%"),
         y = paste0(y_axis, ": ", stdev_val[as.numeric(str_extract(y_axis, "[0-9]"))], "%"))

}


#' @export

draw_volcano <- function(data_frame,
                         fdr_cutoff = 0.05,
                         logfc_cutoff = 1) {

  highlight <- subset(data_frame, FDR < fdr_cutoff & abs(logFC) > logfc_cutoff)

  ggplot(data_frame, aes(logFC, -log10(FDR))) +
    geom_point(shape = 21, cex = 2.4, alpha = 0.7, fill = "gray80") +
    geom_point(data = highlight, shape = 21, cex = 2.5, fill = "tomato") +
    ggrepel::geom_text_repel(data = highlight, label = rownames(highlight))

}








