
#' Stacked columns and pie chart for single cell cytometry data
#'
#'
#' @description  Generates stacked column and pie charts plots for cell-cluster abundance in clustered single-cell cytometry data.
#'
#' @param SE               \code{SummarizedExperiment} object or path to \code{SummarizedExperiment} object.
#' @param groups           Experimental groups to generate the plots for.
#' @param output_folder    Path to folder to receive pdf file of the plots.
#'
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr count
#' @importFrom dplyr ungroup
#' @importFrom dplyr n
#' @importFrom dplyr summarise
#' @import ggplot2
#'
#' @return Saved pdf file containing the plots.
#' @export
#'
#'
#'
#'
#'
cyto_clust_colsPie <- function (SE,
                                groups = "groups",
                                output_folder = "."
) {


  dir.create(file.path(output_folder, "Tables"), showWarnings = FALSE)
  df_for_plotting <- as.data.frame(rowData(SE))

  # counting cells per group - cluster_id combination

  df <- df_for_plotting %>%
    select(!!sym(groups), !!sym("cluster_id")) %>%
    group_by(!!sym(groups), !!sym("cluster_id")) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    group_by(!!sym(groups)) %>%
    mutate("total_count_group" = sum(count)) %>%
    mutate("cluster_percent" = count/!!sym("total_count_group") *100)

  df <- as.data.frame(df)

  write.csv(df, row.names = FALSE, file.path(output_folder, "Tables", "clusters_proportions_per_group.csv"))

  # Colour palatte
  cluster_number <- length(levels(rowData(SE)$cluster_id))
  set.seed(723451)
  col_pal <- unname(Polychrome::createPalette(cluster_number,  c("#9c0505", "#04c204", "#0202bf")))
  set.seed(NULL)

  message("Generating stacked column plot for groups...")
  p_stack <- ggplot(df, aes(x= "", y= .data$cluster_percent, fill= .data$cluster_id)) +
    geom_bar(width = 0.4, position=position_stack(), stat = "identity") +
    scale_fill_manual(values = col_pal) +
    xlab("Groups") + ylab("Cluster %") +
    facet_wrap(sym(groups)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          strip.text = element_text(size=20, face = "bold")) +
    guides(fill=guide_legend(title = "Cluster ID"))


  ggsave("stacked column plot.pdf", plot = p_stack, device = "pdf", path= output_folder)

  message("Generating pie chart for groups...")
  p_pie <- p_stack + coord_polar(theta="y") + xlab("") + ylab("")

  ggsave("pie chart.pdf", plot = p_pie, device = "pdf", path= output_folder)


}

