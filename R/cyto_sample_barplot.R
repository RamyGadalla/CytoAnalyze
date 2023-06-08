#' Single cell cytometry data visualization
#'
#'
#' @description Plot a stacked columns plot for the proportion of clusters in each sample in clustered cytometry data.
#'
#' @param SE               \code{SummarizedExperiment} object or path to \code{SummarizedExperiment} object.
#' @param samples          Character. It must be in \code{colnames} \code{(rowData())}, indicates the sample name/ID.
#' @param reorder          A variable in the dataframe to use to rearrange levels of "samples" for visualization.
#' @param output_folder    Path to folder to receive pdf file of the columns plot.
#'
#' @return
#' stacked columns plot for cluster proportions within all samples in the dataset.
#'
#'
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#'
#' @return Barplot for cluster proportions in each sample in the data, and csv file with this data.
#'
#'
#'
#' @export
#'
cyto_samples_barplot <- function (SE,
                                 samples = "sample_id",
                                 reorder = NULL,
                                 output_folder = ".") {



  if (is.character(SE)) {
    SE <- readRDS(SE)
  }

  if (!is(SE, "SummarizedExperiment")) {
    stop('SE needs to be a SummarizedExperiment object or path to SummarizedExperiment object')
  }


  dir.create(file.path(output_folder, "Tables"), showWarnings = FALSE)
  df_for_plotting <- as.data.frame(rowData(SE))


  if (!is.null(reorder)) {
    df_for_plotting$patient_id <- factor(df_for_plotting[,samples], levels = unique(df_for_plotting$reorder[order(df_for_plotting[,reorder])]))

  }

  # counting cells per group - cluster_id combination

  df <- df_for_plotting %>%
    select(!!sym(samples), !!sym("cluster_id")) %>%
    group_by(!!sym(samples), !!sym("cluster_id")) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    group_by(!!sym(samples)) %>%
    mutate("total_count_group" = sum(count)) %>%
    mutate("cluster_percent" = count/!!sym("total_count_group") *100)

  write.csv(subset(df, select=-total_count_group), row.names = FALSE, file.path(output_folder, "Tables", "clusters_proportions_per_sample.csv"))

  # Colour palatte
  cluster_number <- length(levels(rowData(SE)$cluster_id))
  set.seed(723451)
  col_pal <- unname(Polychrome::createPalette(cluster_number,  c("#9c0505", "#04c204", "#0202bf")))
  set.seed(NULL)

  message("Generating samples-cluster proportions barplot...")



  p <- ggplot(df, aes(y=.data$cluster_percent, x=!!sym(samples))) +
    geom_col(aes(fill=.data$cluster_id), position = "fill") +
    ylab("proportion") + xlab("samples") +
    theme(axis.line=element_line(color="black"), axis.text.x = element_text(angle = 90, size=10, face="bold"),
          axis.title.x = element_text(size= 12, face="bold"), panel.background = element_blank(),
          axis.title.y = element_text(size=12, face="bold")) +
    scale_fill_manual(values = col_pal)

  ggsave("sample_cluster_col_plot.pdf", plot = p, device="pdf", path=file.path(output_folder), width=8, height=12)

}
