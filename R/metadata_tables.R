#' Extract metadata tables from raw fcs files
#'
#' @description Extract two metadata tables from fcs files to be completed by the analyst for downstream analysis. First table is the experimental metadata (experiment_info).
#'              Second table is the markers used metadata (marker_info).
#'
#'
#' @param input_folder    Path to folder containing raw fcs files.
#' @param output_folder   Path to folder to receive the tables.
#' @param clean_name      If TRUE, the returned marker_info table will have the biological marker names cleaned and readable in column "clean_name". If FALSE (default)
#'                        the column will be empty and the analyst will be expected to fill this column manually. The name cleaning is done assuming the conventional and the most common marker naming
#'                        during the data acquisition. If not sure, leave it to the default value.
#' @import flowCore stringi
#' @importFrom utils read.csv write.csv
#'
#' @return                csv files exported to output_folder.
#'                        Analyst are required to fill these tables. For marker_info table, column $marker_class should have only one of the three values c("none", "type", "state), and column
#'                        $channel_type should only have one of the values c("biological', "non_biological"). Failure to do that will cause errors downsteam.
#'
#'
#'
#' @note clean_name = TRUE will remove the common used characters in naming channels to make the biomarkers more easily identified. But it
#' could result in some unclear names occasionally. So, the analyst's attention and edits is recommended after table is exported. Generally, clean-names
#' should not contain any special characters and should not start with a number.
#'
#' @export

metadata_tables <- function (input_folder,
                             output_folder,
                             clean_name = FALSE) {

  dir.create(file.path(output_folder, "metadata_tables"), showWarnings = FALSE)

  fcs_files <- list.files(path = input_folder, pattern = ".fcs$", full.names=TRUE)

  # write out experiment_info table that will be filled with experimental details.
  #read the whole data set to make sure all files have same panel.
  all_files_set <- flowCore::read.flowSet(fcs_files, transformation = FALSE, truncate_max_range = FALSE)
  experiment_info <- flowCore::pData(all_files_set)
  colnames(experiment_info) <- "sample_id"
  experiment_info$groups
  write.csv(experiment_info, row.names = FALSE, file.path(output_folder, "metadata_tables", "experiment_info.csv"))


  # write out marker_info table that will contain information about the markers panels used.

  marker_info <-
    flowCore::pData(flowCore::read.FCS(fcs_files[[1]], transformation=FALSE, truncate_max_range=FALSE)@parameters)[,c("name", "desc")]


  marker_info$clean_name <- ""

  if (clean_name) {
         marker_info$clean_name <- stringi::stri_replace_all_regex (marker_info$desc,
                                                                    pattern =c("^(.*?)_", "EQ", ":", "::", "-"),
                                                                    replacement = c(""),
                                                                    vectorize_all = FALSE
                                                                    )

  }


    colnames(marker_info)[2] <- "marker_name"
    marker_info$marker_class <- ""
    marker_info$channel_type <- ""

    write.csv(marker_info, row.names = FALSE, file.path(output_folder, "metadata_tables", "marker_info.csv"))

    message("metadata tables exported to the specified output folder.")

}


