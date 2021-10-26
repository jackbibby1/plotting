#' Draw PCA plot from data frame
#'
#' @param data_frame data frame with samples in columns
#' @param log should the data be log2 transformed?
#' @param group_by groups present in the samples supplied as character vector
#' @param x_axis principle components to be used on the x-axis
#' @param y_axis principle components to be used on the y-axis
#' @param return_data should the PCA data frame be returned
#'
#' @examples \dontrun{
#' draw_pca(expression_matrix, group_by = groups, return_data = F)
#' }
#'
#' @export

draw_pca <- function(data_frame,
                     log = T,
                     group_by = NULL,
                     x_axis = "PC1",
                     y_axis = "PC2",
                     return_data = T) {

  if (log == T) {
    pca_dat <- stats::prcomp(t(log2(data_frame)))
  } else {
    pca_dat <- stats::prcomp(t(data_frame))
  }

  pca_df <- data.frame(pc1 = pca_dat$x[, x_axis],
                       pc2 = pca_dat$x[, y_axis],
                       groups = group_by)

  if (return_data == T) {
    return(pca_df)
  }

  stdev_val <- pca_dat$sdev^2
  stdev_val <- round(stdev_val/sum(stdev_val)*100, 1)

  ggplot2::ggplot(pca_df, ggplot2::aes(pc1, pc2)) +
    ggplot2::geom_point(cex = 3.3, shape = 21, ggplot2::aes(fill = groups), alpha = 0.8, stroke = 0.3) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(fill = NA),
          aspect.ratio = 1,
          legend.key = ggplot2::element_blank()) +
    ggplot2::labs(x = paste0(x_axis, ": ", stdev_val[as.numeric(stringr::str_extract(x_axis, "[0-9]"))], "%"),
                  y = paste0(y_axis, ": ", stdev_val[as.numeric(stringr::str_extract(y_axis, "[0-9]"))], "%"))

}

#' Draw a volcano plot with FDR and logFC values
#'
#' @param data_frame Data frame with logFC, FDR, and gene columns at minimum
#' @param fdr_cutoff FDR cutoff for significant points
#' @param logfc_cutoff LogFC cutoff for significant points
#'
#' @examples \dontrun{
#' draw_volcano(degs, fdr_cutoff = 0.01, logfc_cutoff = 0.7)
#' }
#'
#' @export

draw_volcano <- function(data_frame,
                         fdr_cutoff = 0.05,
                         logfc_cutoff = 1) {

  highlight <- subset(data_frame, FDR < fdr_cutoff & abs(logFC) > logfc_cutoff)

  ggplot2::ggplot(data_frame, ggplot2::aes(logFC, -log10(FDR))) +
    ggplot2::geom_point(shape = 21, cex = 2.4, alpha = 0.7, fill = "gray80") +
    ggplot2::geom_point(data = highlight, shape = 21, cex = 2.5, fill = "tomato") +
    ggrepel::geom_text_repel(data = highlight, label = rownames(highlight))

}









