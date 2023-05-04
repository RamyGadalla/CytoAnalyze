#' Drop clusters of single cell cytometry data
#'
#' @description Remove specific clusters from SummarizedExperiment object
#'
#' @param SE               SummarizedExperiment object or path to SummarizedExperiment object.
#' @param clusters         Vector of the cluster IDs that needs to be removed.
#' @param reassign         Logical. If TRUE, the cluster IDs are reassigned new IDs starting from 1.
#' @param output_folder    Path of output directory to receive the subsetted SummarizedExperiment.
#' @param SE_name          Name of the new subsetted SummarizedExperiment.
#'
#' @return
#' @export
#'
#' @examples
cluster_drop <- function(SE,
                          clusters,
                          reassign = FALSE,
                          output_folder = NULL,
                          SE_name = NULL) {

if (is.character(SE)) {
  SE <- readRDS(SE)
}

if (!is(SE, "SummarizedExperiment")) {
  stop('SE needs to be a SummarizedExperiment object or path to SummarizedExperiment object')
}

for (cluster in clusters) {
SE <- SE[rowData(SE)$cluster_id != cluster, ]
}

rowData(SE)$cluster_id <- droplevels(rowData(SE)$cluster_id)

# Reassign cluster_ids
if (reassign){
  levels(rowData(SE)$cluster_id) <- c(1:length(levels(rowData(SE)$cluster_id)))
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

