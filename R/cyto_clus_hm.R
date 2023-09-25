#' Heatmap for clustered single-cell cytometry data
#'
#'
#' @description A wrapper around \code{\link{pheatmap}} package to facilitate  generation of heatmap for clustered single-cell cytoemetry data.
#'
#' @param SE              \code{SummarizedExperiment} object or path to \code{SummarizedExperiment} object.
#' @param output_folder   Path to folder to recieve pdf file of the heatmap.
#' @param assay_name      Name of \code{assay()} slot to be used in generating the heatmap.
#' @param markers         A numeric vector specifying the indices of markers to be plotted. By default, it plots all the markers that are not labelled as "none" in "marker_class" column in \code{colData(SE)}.
#' @param ...             Extra arguments for \code{pheatmap()}
#'
#' @return Heatmap clusters as rows and median protein expression as column.
#'
#' @note The protein expression is scaled across the clusters.
#'
#' @importFrom Polychrome createPalette
#' @importFrom pheatmap pheatmap
#' @importFrom MatrixGenerics colMedians
#' @importFrom RColorBrewer brewer.pal
#'
#'
#' @export
#'
#'
#'
cyto_clust_hm <- function (SE,
                          output_folder = ".",
                          assay_name = "exprs",
                          markers = NA,
                          ...) {

  if (is.character(SE)) {
    SE <- readRDS(SE)
  }

  if (!is(SE, "SummarizedExperiment")) {
    stop('SE needs to be a SummarizedExperiment object or path to SummarizedExperiment object')
  }


  dir.create(file.path(output_folder, "Tables"), showWarnings = FALSE)

  #building data frame

  if (length(markers) == 1 & all(is.na(markers))) {
    marker_flag <- which(colData(SE)$marker_class != "none")
  } else if (length(markers) >= 1 & all(!is.na(markers))) {
    marker_flag <- markers
  }


  cluster_names <- levels(rowData(SE)$cluster_id)
  df <- data.frame(matrix(ncol=length(marker_flag), nrow=length(cluster_names)))
  colnames(df) <- colnames(SE)[marker_flag]
  rownames(df) <- cluster_names


  #calculating Median expression for each cluster
  i=1
  for(cluster in cluster_names)
  {
    cluster_rows <- which(rowData(SE)$cluster_id == cluster)
    cluster_assays <- assay(SE[cluster_rows,], assay_name)
    cluster_assays <- cluster_assays[,marker_flag]

    df[i,] <- colMedians(data.matrix(cluster_assays))
    i = i+1
  }


  write.csv(df, row.names = FALSE, file.path(output_folder, "Tables", "heatmap_markers_median.csv"))

  #Annotation dataframe

  clusters <- unique(rowData(SE)$cluster_id)
  ordered <- order(clusters)

  annot_df_col <- data.frame(row.names =rownames(df), Cluster_id = as.factor(ordered[clusters]))

  cluster_number <- length(levels(rowData(SE)$cluster_id))
  set.seed(723451)
  col_pal <- unname(Polychrome::createPalette(cluster_number,  c("#9c0505", "#04c204", "#0202bf")))
  set.seed(NULL)

  #Colors annotation
  names(col_pal) <- annot_df_col$Cluster_id



  message("Generating heatmap plot...")

  pheatmap(t(df),
                   scale = "column",
                   color =rev(brewer.pal(n = 7, name ="RdBu")),
                   border_color = "black",
                   cluster_rows=TRUE,
                   cluster_cols = TRUE,
                   annotation_col = annot_df_col,
                   annotation_colors = list(Cluster_id = col_pal),
                   main = "Marker Expression Intensity",
                   angle_col = 45,
                   height = 12,
                   width = 12,
                   filename = file.path(output_folder, "clusters_heatmap.pdf")
  )

}
