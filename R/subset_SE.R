#' Subsetting SummarizedExperiment
#'
#' @description Subsets SummarizedExperiment object for CytoAnalyze analysis pipeline.
#'
#' @param SE                \code{SummarizedExperiment} object or path to \code{SummarizedExperiment} object.
#' @param group             Character. Column name in \code{rowData()} to subset from.
#' @param no_keep           Character vector of the levels in \code{group} that need to be removed.
#' @param keep              Character vector of the levels in \code{group} that need to be kept.
#' @param output_folder     Path of output directory to recieve the subsetted SummarizedExperiment.
#' @param SE_name           Name of the new subsetted SummarizedExperiment.
#'
#' @return
#' SummarizedExperiment
#'
#'
#'
#' @export
#'
#'
subset_SE <- function (SE,
                       group,
                       no_keep,
                       keep = NULL,
                       output_folder,
                       SE_name) {

if (is.character(SE)) {
  SE <- readRDS(SE)
}

if (!is(SE, "SummarizedExperiment")) {
  stop('SE needs to be a SummarizedExperiment object or path to SummarizedExperiment object')
}


  if(!group %in% colnames(rowData(SE))) {
    stop('Provided "groups" is not in the data.')
  }

if(!no_keep %in% levels(rowData(SE)[,group])) {
  stop(paste(no_keep, "is not a level in", group))
}

if(!is.factor(rowData(SE)[,group])) {
  stop(paste(group, "must be a factor type."))
}


if (is.null(keep) & is.null(no_keep)) {
  stop("No levels are provided to subset.")
}



if (!is.null(keep) & !is.null(no_keep)) {

  for (level in no_keep) {

    SE_temp <- SE[rowData(SE)[group] != level,]

  }

  for (level in keep) {

    SE_temp <- SE[rowData(SE)[group] == level,]

  }

}


if (!is.null(no_keep) & is.null(keep)) {

  for (level in no_keep) {

    SE_temp <- SE[rowData(SE)[,group] != level,]

  }
}


if (!is.null(keep) & is.null(no_keep)) {

  for (level in keep) {

    SE_temp <- SE[rowData(SE)[group] == level,]

  }
}


#drop_levels
# have to loop here over all columns in rowdata that are factors

factor_cols <- colnames(rowData(SE)[,sapply(rowData(SE), is.factor)])

for (cols in factor_cols) {

  rowData(SE_temp)[,cols] <- droplevels(rowData(SE_temp)[,cols])

  levels(rowData(SE_temp)[,cols])
}



#metadata

experiment_info <- metadata(SE_temp)$experiment_info[metadata(SE_temp)$experiment_info$sample_id %in% levels(rowData(SE_temp)$sample_id),]

num_cells <- metadata(SE_temp)$n_cells[names(metadata(SE_temp)$n_cells) %in% levels(rowData(SE_temp)$sample_id)]


# Reconstruct SE


SE <- SummarizedExperiment(assay = assay(SE_temp),
                            rowData = rowData(SE_temp),
                            colData = colData(SE_temp),
                            metadata = list(
                              experiment_info = experiment_info,
                              n_cells = num_cells))

names(assays(SE)) <- "exprs"

if(!is.null(output_folder)) {
  saveRDS(SE, file.path(output_folder, paste0(SE_name, "_SE.rds")))
}

return(SE)


}

