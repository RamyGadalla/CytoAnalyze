#'
#' Violin plots for single-cell cytometry data
#'
#' @description Generates violin plots of marker expression intensity on cellular clustered data.
#'
#' @param SE             \code{SummarizedExperiment} object or path to \code{SummarizedExperiment} object.
#' @param assay_name     Name of \code{assay()} slot to be used in generating the violin plots.
#' @param marker_list    Character vector containing the name of the markers to generate violin plots for. It must be included
#'                       in \code{clean_name} column of \code{colData(SE)}
#' @param groups         Experimental groups to generate the plots for.
#' @param output_folder  Path to folder to receive pdf file of the violin plots.
#'
#'
#' @importFrom Polychrome createPalette
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr all_of
#'
#' @import ggplot2
#'
#' @return Saved pdf file containing violin plots.
#' @export
#'
#'

cyto_clust_violin <- function (SE,
                               assay_name = "exprs",
                               marker_list = "all",
                               groups = NULL,
                               output_folder = ".") {

  message("Violin plots...")
  rownames(colData(SE)) <- colData(SE)$clean_name

  if (length(marker_list) == 1 & all(marker_list == "all")) {

    marker_list <- colData(SE)[colData(SE)$channel_type == "biological","clean_name"]
    df_for_plotting <- data.frame(rowData(SE), assay(SE, assay_name)[,as.character(marker_list)], check.names = FALSE)

  } else if (length(marker_list) == 1) {

    df_for_plotting <- data.frame(rowData(SE), assay(SE, assay_name)[,as.character(marker_list)], check.names = FALSE, fix.empty.names = FALSE)
    colnames(df_for_plotting)[ncol(df_for_plotting)] <- marker_list

  } else if (length(marker_list) >= 1) {

    marker_list <- as.character(marker_list)
    df_for_plotting <- data.frame(rowData(SE), assay(SE, assay_name)[,as.character(marker_list)], check.names = FALSE)

  }


  # Colour palatte
  cluster_number <- length(levels(rowData(SE)$cluster_id))
  set.seed(723451)
  col_pal <- unname(Polychrome::createPalette(cluster_number,  c("#9c0505", "#04c204", "#0202bf")))
  set.seed(NULL)

  dir.create(file.path(output_folder, "violins"), showWarnings = FALSE)


    df <- df_for_plotting %>% select("cluster_id", all_of(marker_list))

    dir.create(file.path(output_folder, "violins", "Overall markers intensity"), showWarnings = FALSE)

    for (marker in marker_list) {

      message("Generating violin plots...", marker)
      p_violin <- suppressWarnings(ggplot(df, aes_string(x= "cluster_id", y= df[,marker])) +
        geom_violin(width=3, aes(fill = .data$cluster_id)) +
        scale_fill_manual(values = col_pal) +
        theme(panel.grid= element_blank(), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
              axis.text.x=element_text(size=15, color = "black"),
              axis.text.y=element_text(size=10, color = "black"),
              axis.line = element_line(color = "black"),
              strip.text.x = element_text(size = 17, face = "bold"),
              strip.background = element_rect(color = "black", fill = "#edece8")) +
        xlab("cluster_id") + ylab(marker) +
        geom_boxplot(width=0.05, outlier.shape = NA) +
        guides(fill=guide_legend(title = "Cluster ID")))

      ggsave(paste0(marker,"_violin.pdf"), plot = p_violin, device = "pdf", path= file.path(output_folder, "violins", "Overall markers intensity"), width = 25, height = 6)

    }


  if (!is.null(groups)) {

    df <- df_for_plotting %>% select("cluster_id", groups, all_of(marker_list))

    dir.create(file.path(output_folder, "violins", "Markers intensity per group"), showWarnings = FALSE)

    for (marker in marker_list) {

      message("Generating violin plots...", marker, " for ", groups, " levels.")
      p_violin_groups <- suppressWarnings(ggplot(df, aes_string(x= groups, y= df[,marker], fill=groups)) +
        geom_violin(width =1) +
        facet_grid(. ~ .data$cluster_id) +
        scale_fill_brewer(palette = "Dark2") +
        theme(panel.grid= element_blank(), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
              axis.text.x=element_text(size=15, color = "black", angle = 90),
              axis.text.y=element_text(size=10, color = "black"),
              axis.line = element_line(color = "black"),
              strip.text.x = element_text(size = 17, face = "bold"),
              strip.background = element_rect(color = "black", fill = "#edece8")) +
        xlab(groups) + ylab(marker) +
        geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white"))



      ggsave(paste0(marker,"_groups_violin.pdf"), plot = p_violin_groups, device = "pdf", path= file.path(output_folder, "violins", "Markers intensity per group"), width = 35, height = 6)

    }
  }

}

