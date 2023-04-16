
#' Calculate and plot correlation of cluster abundance.


#' @description Calculates and plots correlation of cluster abundance in clustered single-cell cytometry data.
#'
#' @param SE                \code{SummarizedExperiment} object or path to \code{SummarizedExperiment} object.
#' @param assay_name        A character vector of length 1 indicates the assay to be used in plotting the UMAP. By default, it is \code{"exprs"}.
#' @param groups            Experimental groups to generate the plots for.
#' @param output_folder     Path to folder to receive pdf file of the plots.
#' @param ...               Extra arguments passed to \code{\link[corrplot]{corrplot}}
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr count
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom dplyr n
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom stats cor
#'
#' @importFrom tidyr spread
#' @importFrom Hmisc rcorr
#'
#'
#' @return  correlation plots between cell-clusters for clustered cytometry single cell data.
#' @export
#'
#'
#'
#'
#'
cyto_clust_corr <- function (SE,
                      assay_name = "exprs",
                      groups = NULL,
                      output_folder = ".",
                      ...
                      ) {

  if (is.character(SE)) {
    SE <- readRDS(SE)
  }

  if (!is(SE, "SummarizedExperiment")) {
    stop('SE needs to be a SummarizedExperiment object or path to SummarizedExperiment object')
  }


  dir.create(file.path(output_folder, "Tables"), showWarnings = FALSE)

  df_for_plotting <- as.data.frame(rowData(SE))

  message("calculating clusters correlation in samples...")
  df <- df_for_plotting %>%
    select(!!sym("sample_id"), !!sym("cluster_id")) %>%
    group_by(!!sym("sample_id") ,!!sym("cluster_id")) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    group_by(!!sym("sample_id")) %>%
    mutate("total_count_group" = sum(count)) %>%
    mutate("cluster_percent" = count/!!sym("total_count_group") *100) %>%
    ungroup() %>%
    select("sample_id", "cluster_id", "cluster_percent")

  df <- tidyr::spread(df, "cluster_id", "cluster_percent", fill = 0)

  correlation_Matrix <- cor(as.matrix(df[,-1]) ,method = "pearson", use ="complete.obs")

  correlation_Matrix_Significance <-  Hmisc::rcorr(as.matrix(df[,-1]))

  pvals <- correlation_Matrix_Significance$P
  pvals[is.na(pvals)] <- 1


  write.csv(correlation_Matrix, row.names = FALSE, file.path(output_folder, "Tables", "samples correlation matrix.csv"))
  write.csv(as.matrix(correlation_Matrix_Significance$P), row.names = FALSE, file.path(output_folder, "Tables", "samples correlation matrix pvals.csv" ))


  pdf(file = file.path(output_folder, "cluster abundaces correlation in samples.pdf"))

  corrplot::corrplot(correlation_Matrix,
                     type="upper",
                     order="hclust",
                     tl.col="black",
                     tl.srt=90,
                     tl.cex=0.7,
                     p.mat = pvals,
                     pch=3,
                     pch.cex = 1,
                     insig = "label_sig",
                     bg=  "#edece8",
                     mar=c(0,0,2,0),
                     title = "Clusters correlation in all samples")

  dev.off()

  if (!is.null(groups)) {
    # Correlation within groups

    message("Calculating clusters correlation in groups...")
    # split main dataframe into sub dataframe based on groups level and store them in a list.
    groups_level <- levels(df_for_plotting[,groups])
    levels_list <- list()

    for (group in groups_level) {
      df_group <- df_for_plotting[df_for_plotting[groups]==group,]
      levels_list[[group]] <- df_group
    }

    for (group in names(levels_list)) {

      if(length(unique(levels_list[[group]]$sample_id)) <= 4) {
        message("\033[31mCan't calculate correlation matrix for group ", group, " as it has only 4 or less samples\033[0m")
        next
      }

      df <- levels_list[[group]] %>%
        select(!!sym("sample_id"), !!sym("cluster_id")) %>%
        group_by(!!sym("sample_id") ,!!sym("cluster_id")) %>%
        summarise(count = n()) %>%
        ungroup() %>%
        group_by(!!sym("sample_id")) %>%
        mutate("total_count_group" = sum(count)) %>%
        mutate("cluster_percent" = count/!!sym("total_count_group") *100) %>%
        ungroup() %>%
        select("sample_id", "cluster_id", "cluster_percent")

      df <- tidyr::spread(df, "cluster_id", "cluster_percent", fill = 0)

      correlation_Matrix <- cor(as.matrix(df[,-1]) ,method="pearson",use="complete.obs")

      correlation_Matrix_Significance <-  Hmisc::rcorr(as.matrix(df[,-1]))

      pvals <- correlation_Matrix_Significance$P
      pvals[is.na(pvals)] <- 1

      write.csv(correlation_Matrix, row.names = FALSE, file.path(output_folder, "Tables", paste0(group, " correlation matrix.csv" )))
      write.csv(as.matrix(correlation_Matrix_Significance$P), row.names = FALSE, file.path(output_folder, "Tables", paste0(group, " correlation matrix pvals.csv" )))

      pdf(file = file.path(output_folder, paste0("cluster abundaces correlation in group ", group, ".pdf")))

      corrplot::corrplot(correlation_Matrix,
                         type = "upper",
                         order = "hclust",
                         tl.col = "black",
                         tl.srt= 90,
                         tl.cex= 0.7,
                         p.mat = pvals,
                         pch = 3,
                         pch.cex = 1,
                         insig = "label_sig",
                         bg =  "#edece8",
                         mar = c(0,0,2,0),
                         title = paste0("Clusters correlation in group: ", group))



      dev.off()

    }
  }

  message("Plot(s) saved at ", output_folder)
}
