# drop clusters

drop_clusters <- function(SE,
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

