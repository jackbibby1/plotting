#' Draw PCA plot from data frame
#'
#' @param df data frame with samples in columns
#' @param log should the data be log2 transformed?
#' @param grouping groups present in the samples supplied as character vector
#' @param x_axis principle components to be used on the x-axis
#' @param y_axis principle components to be used on the y-axis
#' @param return_data should the PCA data frame be returned
#'
#' @examples \dontrun{
#' draw_pca(expression_matrix, grouping = groups, return_data = F)
#' }
#'
#' @export

draw_pca <- function(df,
                     log = T,
                     grouping = NULL,
                     x_axis = "PC1",
                     y_axis = "PC2",
                     return_data = F) {

  if (log == T) {
    pca_dat <- stats::prcomp(t(log2(df)))
  } else {
    pca_dat <- stats::prcomp(t(df))
  }

  pca_df <- data.frame(pc1 = pca_dat$x[, x_axis],
                       pc2 = pca_dat$x[, y_axis],
                       groups = grouping)

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
#' @param df Data frame with logFC, FDR, and gene columns at minimum
#' @param fdr_cutoff FDR cutoff for significant points
#' @param logfc_cutoff LogFC cutoff for significant points
#'
#' @examples \dontrun{
#' draw_volcano(degs, fdr_cutoff = 0.01, logfc_cutoff = 0.7)
#' }
#'
#' @export

draw_volcano <- function(df,
                         fdr_cutoff = 0.05,
                         logfc_cutoff = 1) {

  highlight <- subset(df, FDR < fdr_cutoff & abs(logFC) > logfc_cutoff)

  ggplot2::ggplot(df, ggplot2::aes(logFC, -log10(FDR))) +
    ggplot2::geom_point(shape = 21, cex = 2.4, alpha = 0.7, fill = "gray80") +
    ggplot2::geom_point(data = highlight, shape = 21, cex = 2.5, fill = "tomato") +
    ggrepel::geom_text_repel(data = highlight, label = rownames(highlight))

}


#' Boxplot of gene from expression file, by group
#'
#' @param df data frame of expression matrix with genes as rownames
#' @param gene gene of interest
#' @param grouping groups for condition, can be as a factor
#'
#' @importFrom magrittr "%>%"
#'
#' @export


draw_boxplot <- function(df, gene, grouping) {

  use_gene <- df %>%
    tibble::rownames_to_column("gene_symbol") %>%
    dplyr::filter(gene_symbol == gene) %>%
    dplyr::pull(gene_symbol)

  samples <- colnames(df)

  df[use_gene, ] %>%
    tidyr::pivot_longer(cols = samples, names_to = "sample", values_to = "expression") %>%
    dplyr::mutate(groups = grouping) %>%
    ggplot2::ggplot(ggplot2::aes(groups, expression)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(shape = 21, cex = 2.5, color = 'black', width = 0.2, height = 0, ggplot2::aes(fill = groups)) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::labs(y = "Norm intensity") +
    ggplot2::ggtitle(use_gene)
}



