#' Merge clusters of single-cell cytometry data
#'
#' @param SE                 SummarizedExperiment object or path to SummarizedExperiment object.
#' @param clusters           Vector of the cluster IDs that needs to be merged. All th clusters will be merged in the first cluster in the vector.
#'                           c(5,8,23), cluster 8 and 23 will be merged in cluster 5.
#' @param reassign           Logical. If TRUE, the cluster IDs are reassigned new IDs starting from 1.
#' @param output_folder      Path of output directory to receive the subsetted SummarizedExperiment.
#' @param SE_name            Name of the new subsetted SummarizedExperiment.
#'
#' @return
#' SummarizedExperiment object with the chosen clusters merged.
#'
#' @details
#' Occasionally, analyst might find that some clusters are over-clustered and wish to merge these clusters togther.
#'
#' @export
#'
#'
#'
#'
cluster_merge <- function (SE,
                           clusters,
                           reassign = TRUE,
                           output_folder = NULL,
                           SE_name = NULL) {


  if (is.character(SE)){
    readRDS(SE)
  }

  if (!is(SE, "SummarizedExperiment")){
    stop("SE must be a summarizedExperiment object")
  }

  if(!"cluster_id" %in% colnames(rowData(SE))) {
    stop("rowData(SE) does not contain column 'cluster_id'.")
  }


  for (i in 1:length(clusters)) {

    rowData(SE)$cluster_id[rowData(SE)$cluster_id == clusters[i]] <- clusters[1]

  }

  rowData(SE)$cluster_id <- droplevels(rowData(SE)$cluster_id)

  # reassign
  if (reassign) {

    levels(rowData(SE)$cluster_id) <- as.factor(1:length(levels(rowData(SE)$cluster_id)))
  }

  # Save
  if(!is.null(output_folder)) {
    if(is.null(SE_name)){
      stop("output_folder is provided but not SE_name.")
    }
    saveRDS(SE, file.path(output_folder, paste0(SE_name, "_SE.rds")))
  }

  return(SE)
}




