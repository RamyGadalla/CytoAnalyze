#' Insert metadata to SummarizedExperiment object
#'
#' @param SE                   \code{SummarizedExperiment} object or path to \code{SummarizedExperiment} object.
#' @param meta                 Path to csv file that contain the data to insert.
#' @param add_column           Character vector of any length. It represents the column names of the new columns to be inserted.
#' @param meta_column          Character vector of the same length as \code{add_column} of the names of the columns in \code{meta} that will be inserted.
#' @param matching_column      Character vector of length 1. It represents an existing matching column in \code{meta} and \code{rowData()}
#' @param output_folder        Path of output directory to recieve the subsetted SummarizedExperiment.
#' @param SE_name              Name of the new subsetted SummarizedExperiment.
#'
#'
#' @import SummarizedExperiment
#'
#' @return
#' SummarizedExperiment with added data to \code{rowData()}
#'
#' @export
#'
#'
insert_metadata <- function (SE,
                      meta,
                      add_column,
                      meta_column,
                      matching_column,
                      output_folder = NULL,
                      SE_name = NULL) {

if (is.character(SE)) {
  SE <- readRDS(SE)
}

if (!is(SE, "SummarizedExperiment")) {
  stop('SE needs to be a SummarizedExperiment object or path to SummarizedExperiment object')
}


meta <- read.csv(meta)

if (!matching_column %in% colnames(meta) | !matching_column %in% colnames(rowData(SE))) {
  stop("matching_column needs to be in both meta and rowData().")
}


#rowData

for (i in 1:length(add_column)){

rowData(SE)[[add_column[[i]]]] <- as.factor(meta[[meta_column[[i]]]][match(rowData(SE)[[matching_column]], meta[[matching_column]])])

}


# metadata
# experiment_info
for (i in 1:length(add_column)){

  metadata(SE)$experiment_info[[add_column[[i]]]] <- as.factor(meta[[meta_column[[i]]]][match(metadata(SE)$experiment_info[[matching_column]], meta[[matching_column]])])

}

if(!is.null(output_folder)) {
  if(is.null(SE_name)){
    stop("output_folder is provided but not SE_name.")
  }
  saveRDS(SE_Final, file.path(output_folder, paste0(SE_name, cat(add_column, sep = "_"), "_SE.rds")))
}

return(SE)

}


