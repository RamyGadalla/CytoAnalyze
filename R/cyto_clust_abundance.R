#' Differential abundance test for cytometry data using edgeR
#'
#'
#' @description  Performs differential abundance test for cell clusters of single-cell cytometry data using edgeR test.
#'
#' @param SE                  \code{SummarizedExperiment} object or path to \code{SummarizedExperiment} object.
#' @param groups              Experimental groups to generate the plots for. It has to be a variable in the data i.e. in \code{rowdata()}. Levels of this variable will be compared.
#' @param contrast            A character vector of length 2, specify two levels of "groups" variable to build contrast matrix. \code{(contrast[1]-constrast[2])}.
#' @param exp_subject         Default is NULL. If the experimental design is paired, this variable will be used in the formula to build paired design matrix.
#' @param paired              Logical (default is FALSE). If TRUE a paired experimental design is performed. The pairing is done on the exp_subject variable.
#' @param output_folder       Path to folder to receive pdf file of the plots and tables.
#'
#'
#'
#' @importFrom stats formula
#' @importFrom stats model.matrix
#' @importFrom S4Vectors metadata
#' @importFrom limma makeContrasts
#' @importFrom stats reorder
#' @import ggplot2
#'
#' @return
#' A table (csv file) with log fold change and statistical values for the comparison between the two contrasting groups using edgeR. Also, the function outputs a column plot of the log-fold changes of clusters between groups.
#'
#'
#' @seealso \href{https://github.com/lmweber/diffcyt}{diffcyt}, \href{https://bioconductor.org/packages/release/bioc/html/edgeR.html}{edgeR}
#'
#' @export
#'
#'
#'
#'
cyto_clust_abundance <- function(SE,
                                 groups = "groups",
                                 contrast,
                                 exp_subject = NULL,
                                 paired = FALSE,
                                 output_folder = ".") {


  message("Differential abndance edgeR...")

  if (paired) {
    if (is.null(exp_subject)) {
    stop("the exp_subject argument must be specified since paried is TRUE.")
    }
    message("Paired design...")
    design_matrix <- model.matrix(formula(paste0("~0+", groups, "+", exp_subject)), metadata(SE)$experiment_info)
    colnames(design_matrix) <- c(levels(as.factor(metadata(SE)$experiment_info[,groups])), levels(as.factor(metadata(SE)$experiment_info[,exp_subject]))[-1])

  } else {
    design_matrix <- model.matrix(formula(paste0("~0+", groups)), metadata(SE)$experiment_info)
    colnames(design_matrix) <- levels(as.factor(metadata(SE)$experiment_info[,groups]))
  }


  # Generate a contrast matrix
  myargs = list(contrasts = paste0(contrast[1],"-", contrast[2]), levels = design_matrix)
  contrast_matrix <- do.call(makeContrasts, myargs)



  d_counts_SE <- diffcyt::calcCounts(SE)

  res_DA <- diffcyt::testDA_edgeR(d_counts_SE, design_matrix, contrast_matrix, min_cells = 1, min_samples = 1)
  res_DAedgeR <- as.data.frame(rowData(res_DA))


  dir.create(file.path(output_folder, "stats", "Cluster abundance"), recursive = TRUE, showWarnings = FALSE)
  write.csv(res_DAedgeR, row.names = FALSE, file.path(output_folder, "stats", "Cluster abundance", "Differential abundace.csv"))


  # Generate columns plot

  p <- ggplot(res_DAedgeR, aes(x=reorder(.data$cluster_id, - .data$logFC), y = .data$logFC)) +
    geom_col (aes(fill = .data$logFC),color = "black", width = 1) +
    scale_fill_gradient2(low="#4062d9", high="#ed0707", midpoint = 0, mid = "white") +
    geom_hline(yintercept = 0, linetype="solid") +
    geom_hline(yintercept = 1:-1, linetype="dashed") +
    labs(title = "Clusters Differential Abundance - log fold change", subtitle = paste(contrast[1], "vs", contrast[2]), x = "cluster id") +
    theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank(), plot.background = element_rect(color = "white"),
          panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),
          axis.text.x = element_text(size=10, angle = 45, hjust=1),
          plot.subtitle = element_text(hjust = 0.5), plot.margin = margin(1,2,1,2, "cm")) +
    geom_text(aes(label = ifelse(.data$p_adj < 0.05, "*", "")), size=8)


  plot(p)

  ggsave("column_plot_cluster_logFC.pdf", plot = p, device = "pdf",
         path = file.path(output_folder, "stats", "Cluster abundance"), width = 15, height=8)

  message("Done differential abundace.\nPlot is saved at ", file.path(output_folder, "stats", "Cluster abundance"))
}

