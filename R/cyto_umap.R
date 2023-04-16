#' Visualize phenograph clusters of cytometry single-cells data on UMAP axes.
#'
#'
#' @description Plot phenograph clusters of cytometry single-cells data on UMAP axes.
#'
#' @param SE                   \code{SummarizedExperiment} object or path to  \code{SummarizedExperiment} object.
#' @param output_folder        Path to folder to receive the generated plots.
#' @param groups               Experimental groups to generate plots for.
#' @param marker_int           Logical. If TRUE, UMAPs for markers intensity are generated.
#' @param marker_list          Character vector containing the name of the proteins of interest to generate a UMAP for. It must be included
#'                             in \code{clean_name} column of \code{colData(SE)}
#' @param assay_name           A single character vector indicates the assay to be used in plotting the UMAP. By default, it is \code{"exprs"}.
#' @param plotwrap             Logical. If TRUE, a universal plot of the UMAPs of the proteins listed in \code{marker_list} is generated and saved.
#' @param DS                   Numeric. Number of cells to downsample for visualization only. If \code{groups} is NULL, this value will be used to downsample the whole dataset.
#'                             If groups is not NULL, this value will be used to downsample from each group provided.
#' @param label_size           Numeric. Controls the cluster label size.
#' @param dot_size             Numeric. Controls the plot dot size.
#' @param width                The plot width in inches.
#' @param height               The plot height in inches.
#' @param ...                  Extra arguments to \code{ggplot()}
#'
#'
#
#' @importFrom stats median
#' @importFrom ggrastr geom_point_rast
#' @importFrom rlang sym
#' @importFrom viridis scale_colour_viridis
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr slice_sample
#' @importFrom dplyr select
#' @importFrom dplyr summarize_all
#' @importFrom patchwork wrap_plots
#'
#'
#'
#'
#' @export
#'
#'
#'
cyto_umap <- function (SE,
                       output_folder,
                       groups = NULL,
                       marker_int = TRUE,
                       marker_list ="all",
                       assay_name ="exprs",
                       plotwrap = FALSE,
                       DS = 5000,
                       label_size = 5,
                       dot_size = 0.1,
                       width = 8,
                       height = 8,
                       ...) {

  if (is.character(SE)) {
    SE <- readRDS(SE)
  }

  if (!is(SE, "SummarizedExperiment")) {
    stop('SE needs to be a SummarizedExperiment object or path to SummarizedExperiment object')
  }


  if (is.character(groups)) {
    if(!groups %in% colnames(rowData(SE))) {
      stop('Provided "groups" is not in the data.')
    }
  }

  if (!all(isUnique(colData(SE)$clean_name))) {

    stop("clean_name in colData(SE) must be unique")

  }

  df_for_plotting <- as.data.frame(rowData(SE))

  if (is.null(groups)) {

    dwnsampl_df <- df_for_plotting[sample(nrow(df_for_plotting), DS),]

  } else {

    dwnsampl_df <- df_for_plotting %>% group_by(!!sym(groups)) %>% slice_sample(n = DS)

  }


  # Center the text label on the middle of the cluster.
  centres <- dwnsampl_df %>%
    group_by(!!sym("cluster_id")) %>%
    select(!!sym("UMAP1"),!!sym("UMAP2")) %>%
    summarize_all(median)

  # Colour palatte
  cluster_number <- length(levels(rowData(SE)$cluster_id))
  set.seed(723451)
  col_pal <- unname(Polychrome::createPalette(cluster_number,  c("#9c0505", "#04c204", "#0202bf")))
  set.seed(NULL)

  # UMAP - clusters
  message("Generating universal UMAP plot...")
  p_umap <- suppressWarnings(
    ggplot(dwnsampl_df, aes_string(x="UMAP1", y="UMAP2", colour="cluster_id")) +
      geom_point_rast(size = dot_size) +
      scale_color_manual(values = col_pal) +
      labs(x="UMAP1", y="UMAP2") +
      theme(panel.background = element_blank(), legend.key = element_rect(fill = "transparent", colour = "black"),
            axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"), axis.title = element_text(size=15, face = "bold"),
            strip.background = element_rect(fill = "transparent", colour = "black"),  strip.text = element_text(size=20, face = "bold")) +
      geom_text(data=centres, aes(label = !!sym("cluster_id")),  size = label_size, color = "black") +
      guides(colour = guide_legend(override.aes = list(size=3)))
  )

  #Export the plot
  dir.create(file.path(output_folder, "UMAPs"), showWarnings = FALSE)
  ggsave("universal_umap.pdf", plot = p_umap, device = "pdf", path= file.path(output_folder, "UMAPs"), width=width, height = height)

  if (!is.null(groups)) {
    message("Generating UMAP for groups...")
    p_umap_groups <- p_umap + facet_wrap(groups) + theme(strip.text = element_text(size = 10))
    ggsave("groups_umap.pdf", plot = p_umap_groups, device = "pdf", path= file.path(output_folder, "UMAPs"), width=width, height = height)
  }

  # UMAPs - markers intensities

  if(marker_int == TRUE) {

    rownames(colData(SE)) <- colData(SE)$clean_name
    dir.create(file.path(output_folder, "UMAPs", "Overall markers intensity"), showWarnings = FALSE)

    if (!is.null(groups)) {
      dir.create(file.path(output_folder, "UMAPs", "Markers intensity per group"), showWarnings = FALSE)
    }

    if (length(marker_list) == 1 & all(marker_list == "all")) {

      plot_list_overall <- list ()
      plot_list_groups <- list()

      for (marker in colData(SE)[colData(SE)$marker_class != "none","clean_name"]) {

        message("Generating ", marker, " UMAP...")
        SE_temp <- SE[,colData(SE)$clean_name == marker]
        df <- data.frame(c(rowData(SE_temp), assay(SE_temp, assay_name)), check.names = FALSE)

        p_umap_marker <- suppressWarnings(
          ggplot(df, aes_string(x="UMAP1", y="UMAP2", colour= df[,marker])) +
            geom_point_rast(size = dot_size) +
            labs(x="UMAP1", y="UMAP2") +
            ggtitle(paste0(marker," intensity")) +
            theme(panel.background = element_blank(), legend.key = element_rect(fill = "transparent", colour = "black"),
                  axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"), axis.title = element_text(size=15, face = "bold"),
                  strip.background = element_rect(fill = "transparent", colour = "black"), strip.text = element_text(size=5, face = "bold")) +
            scale_colour_viridis(option="plasma")
        )

        ggsave(paste0(marker,"_umap.pdf"), plot = p_umap_marker, device = "pdf", path= file.path(output_folder, "UMAPs", "Overall markers intensity"), width=width, height = height)

        p_umap_marker -> plot_list_overall[[marker]]

        if (!is.null(groups)) {

          p_umap_marker_group <- p_umap_marker + facet_wrap(groups) + theme(strip.text = element_text(size = 10))
          ggsave(paste0(marker,"_per_group_umap.pdf"), plot = p_umap_marker_group, device = "pdf", path= file.path(output_folder, "UMAPs", "Markers intensity per group"), width=width, height = height)
          p_umap_marker_group -> plot_list_groups[[marker]]
        }
      }

      if (plotwrap) {
        p_wrap <- patchwork::wrap_plots(plot_list_overall, width=2, height = 1)
        ggsave("plot_wrap.pdf", plot = p_wrap, device = "pdf", path= file.path(output_folder, "UMAPs", "Overall markers intensity"))
      }

      if (plotwrap & !is.null(groups)) {
        p_wrap <- patchwork::wrap_plots(plot_list_groups, width=2, height = 1)
        ggsave("plot_wrap.pdf", plot = p_wrap, device = "pdf", path = file.path(output_folder, "UMAPs", "Markers intensity per group"))
      }

    } else if (length(marker_list) >= 1 ) {

      plot_list_overall <- list ()
      plot_list_groups <- list()

      for (marker in marker_list) {

        message("Generating ", marker, " UMAP...")
        SE_temp <- SE[,colData(SE)$clean_name == marker]
        df <- data.frame(c(rowData(SE_temp), assay(SE_temp, assay_name)), check.names = FALSE)

        p_umap_marker <- suppressWarnings(
          ggplot(df, aes_string(x="UMAP1", y="UMAP2", colour= df[,marker])) +
            geom_point_rast(size = dot_size) +
            labs(x="UMAP1", y="UMAP2") +
            ggtitle(paste0(marker," intensity")) +
            theme(panel.background = element_blank(), legend.key = element_rect(fill = "transparent", colour = "black"),
                  axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"), axis.title = element_text(size=15, face = "bold"),
                  strip.background = element_rect(fill = "transparent", colour = "black"), strip.text = element_text(size=5, face = "bold")) +
            scale_colour_viridis(option="plasma")
        )

        ggsave(paste0(marker,"_umap.pdf"), plot = p_umap_marker, device = "pdf", path= file.path(output_folder, "UMAPs", "Overall markers intensity"), width=width, height = height)

        p_umap_marker ->  plot_list_overall[[marker]]

        if (!is.null(groups)) {
          p_umap_marker_group <- p_umap_marker + facet_wrap(groups) + theme(strip.text = element_text(size = 10))
          ggsave(paste0(marker,"_per_group_umap.pdf"), plot = p_umap_marker_group, device = "pdf", path= file.path(output_folder, "UMAPs", "Markers intensity per group"), width=width, height = height)
          p_umap_marker_group -> plot_list_groups[[marker]]
        }
      }

      if (plotwrap) {
        p_wrap <- patchwork::wrap_plots(plot_list_overall)
        ggsave("plot_wrap.pdf", plot = p_wrap, device = "pdf", path= file.path(output_folder, "UMAPs", "Overall markers intensity"))
      }


      if (plotwrap & !is.null(groups)) {
        p_wrap <- patchwork::wrap_plots(plot_list_groups)
        ggsave("plot_wrap.pdf", plot = p_wrap, device = "pdf", path= file.path(output_folder, "UMAPs", "Markers intensity per group"))
      }

    }
  }
}
