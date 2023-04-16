#' Differential expression of biomarkers for clustered single-cell cytometry data using limma
#'
#'
#' @description  Performs differential expression analysis to clustered cytometry data using limma.
#'
#'
#' @param SE                 \code{SummarizedExperiment} object or path to \code{SummarizedExperiment} object.
#' @param contrast           A character vector of length 2, specify two levels of "groups" variable to build contrast matrix. \code{(contrast[1]-constrast[2])}.
#' @param groups             Experimental groups to generate the plots for. It has to be a variable in the data i.e. in \code{rowdata()}. Levels of this variable will be compared.
#' @param exp_subject        Default is NULL. If the experimental design is paired, this variable will be used in the formula to build paired design matrix.
#' @param paired             Logical (default is FALSE). If TRUE a paried experimental design is performed. The pairing is done on the exp_subject variable.
#' @param output_folder      Path to folder to receive limma output table.
#'
#'
#' @importFrom limma makeContrasts
#' @importFrom stats model.matrix
#' @importFrom stats formula
#' @importFrom S4Vectors metadata
#'
#'
#'
#' @return
#' A table (csv file) with log fold change and statistical values for the comparison between the two contrasting groups using limma.
#'
#'
#' @seealso   \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma}.\href{https://github.com/lmweber/diffcyt}{diffcyt}
#'
#'
#' @export
#'
#'
cyto_clust_limma <- function (SE,
                              contrast,
                              groups = "groups",
                              exp_subject = NULL,
                              paired = FALSE,
                              output_folder = "."
                             ) {

message("Differential expression limma...")
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

rownames(colData(SE)) <- colData(SE)$marker_name

d_medians_SE <- diffcyt::calcMedians(SE)
marker_flag <-  metadata(d_medians_SE)$id_state_markers | metadata(d_medians_SE)$id_type_markers

res_DS <- diffcyt::testDS_limma(d_counts_SE, d_medians_SE,
                       design_matrix, contrast_matrix, plot=FALSE, markers_to_test = marker_flag,
                       min_cells = 1, min_samples = 1)
as.data.frame(rowData(res_DS)) -> res_DS_limma


dir.create(file.path(output_folder, "stats", "Differential expression"), recursive = TRUE, showWarnings = FALSE)
write.csv(res_DS_limma, row.names = FALSE, file.path(output_folder, "stats", "Differential expression", "Differential_expression_limma.csv"))

message("Done limma.")

}
