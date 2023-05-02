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
#' @param rm_na             Remove NA values.
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
                       no_keep = NULL,
                       keep = NULL,
                       output_folder,
                       SE_name,
                       rm_na=TRUE) {


  if (is.character(SE)) {
    SE <- readRDS(SE)
  }

  if (!is(SE, "SummarizedExperiment")) {
    stop('SE needs to be a SummarizedExperiment object or path to SummarizedExperiment object')
  }


  if(!group %in% colnames(rowData(SE))) {
    stop('Provided "groups" is not in the data.')
  }

  # NULL can be logical. remove thAT PART OR MAKE IT INCLUDE NULL
  #if(!no_keep %in% levels(rowData(SE)[,group])) {
  #  stop(paste(no_keep, "is not a level in", group))
  #}

  if(!is.factor(rowData(SE)[,group])) {
    stop(paste(group, "must be a factor type."))
  }


  if(is.null(no_keep) & is.null(keep) & !rm_na){
    stop("no parameters provided to subset.")
  }


  if (rm_na) {
    SE <- SE[!is.na(rowData(SE)[,group]),]
  }


  # build empty SE with slot names as input SE

  temp_assay <- matrix(0, nrow = 0, ncol = ncol(assay(SE)))
  colnames(temp_assay) <- colnames (assay(SE))

  temp_rowdata <-  data.frame(matrix(nrow=0, ncol = ncol(rowData(SE))))
  colnames(temp_rowdata) <- colnames(rowData(SE))


  SE_Final <- SummarizedExperiment(assay=temp_assay,
                                   rowData = temp_rowdata,
                                   metadata = list())

  names(assays(SE_Final)) <- "exprs"


  if (!is.null(keep) & !is.null(no_keep)) {

    for (level in no_keep) {

      SE_temp <- SE[rowData(SE)[,group] != level,]
      SE_Final <- rbind(SE_Final,SE_temp)

    }

    for (level in keep) {

      SE_temp <- SE[rowData(SE)[,group] == level,]
      SE_Final<- rbind(SE_Final,SE_temp)

    }

  }


  if (!is.null(no_keep) & is.null(keep)) {

    for (level in no_keep) {

      SE_temp <- SE[rowData(SE)[,group] != level,]
      SE_Final <- rbind(SE_final,SE_temp)

    }
  }


  if (!is.null(keep) & is.null(no_keep)) {

    for (level in keep) {

      SE_temp <- SE[rowData(SE)[,group] == level,]
      SE_Final<- rbind(SE_Final,SE_temp)

    }
  }


  #drop unused levels

  factor_cols <- colnames(rowData(SE_Final)[,sapply(rowData(SE_Final), is.factor)])

  for (cols in factor_cols) {

    rowData(SE_Final)[,cols] <- droplevels(rowData(SE_Final)[,cols])

    levels(rowData(SE_Final)[,cols])
  }



  #metadata

  experiment_info <- metadata(SE)$experiment_info[metadata(SE)$experiment_info$sample_id %in% levels(rowData(SE_Final)$sample_id),]

  num_cells <- metadata(SE)$n_cells[names(metadata(SE)$n_cells) %in% levels(rowData(SE_Final)$sample_id)]


  # Reconstruct SE


  SE_Final <- SummarizedExperiment(assay = assay(SE_Final),
                                   rowData = rowData(SE_Final),
                                   colData = colData(SE_Final),
                                   metadata = list(
                                     experiment_info = experiment_info,
                                     n_cells = num_cells))


  names(assays(SE_Final)) <- "exprs"
  rownames(metadata(SE_Final)$experiment_info) <- as.character(c(1:nrow(metadata(SE_Final)$experiment_info)))

  if(!is.null(output_folder)) {
    saveRDS(SE_Final, file.path(output_folder, paste0(SE_name, "_SE.rds")))
  }

  return(SE_Final)

}

