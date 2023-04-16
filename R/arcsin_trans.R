#' Arcsin data transformation for single-cell cytometry data
#'
#'
#' @description Performs data transformation to cytometry single-cell data.
#'
#' @param input_folder     Path to folder containing fcs files.
#' @param cofactor         Numeric. The value to be transformed is divided by this cofactor before arcsin transformation.
#'                         This value controls the extent of the linear region around 0. Default is 5.
#' @param output_folder    Path to folder to save the transformed fcs files.
#' @param channel.plot     Logical (default is FALSE). Plot each channel before and after transformation for a randam fcs file in the dataset for visual inspection of transformation. Plots are exported to \code{output_folder}.
#'
#' @import patchwork
#' @import flowCore
#' @import ggplot2
#'
#' @return Arcsin transformed fcs files. And, histograms of all channels of one randomly chosen file.
#'
#' @export
#'
#'
#'
#'
arcsin_trans <- function(input_folder,
                         cofactor = 5,
                         output_folder = ".",
                         channel.plot = FALSE) {

  fcs_files <- list.files(path = input_folder, pattern = ".fcs$", full.names=TRUE, recursive = TRUE)


  dir.create(file.path(output_folder,"ArcsinTrans"), showWarnings = FALSE)

  for (file in fcs_files) {

    message("Transforming....", file)

    temp <- flowCore::read.FCS(file, transformation = FALSE, truncate_max_range = FALSE)

    temp@exprs <- asinh(temp@exprs/cofactor)
    #translist <- flowCore::transformList(colnames(temp), arcsinhTransform(a=0, b=1/cofactor))

    #trans_fcs <- flowCore::transform(temp, translist)

    flowCore::write.FCS(temp, file.path(output_folder, "ArcsinTrans", paste0(tools::file_path_sans_ext(basename(file)),"_arcsin_trans.fcs")))


  }

  message("Done transformation \n.\n.")

  if (channel.plot) {

    message("Plotting histograms")

    dir.create(file.path(output_folder,"ArcsinTrans", "Before and after histograms"), showWarnings = FALSE)

    trans_files <- list.files(path = file.path(output_folder, "ArcsinTrans"),
                              pattern = ".fcs$", full.names=TRUE)

    random_file <- sample (c(1:length(fcs_files)), 1)

    untrans <- read.FCS(fcs_files[random_file], transformation = FALSE, truncate_max_range = FALSE)
    trans <- read.FCS(trans_files[random_file], transformation = FALSE, truncate_max_range = FALSE)

    message(fcs_files[random_file], " was chosen randomly for plotting \n.\n.")

    colnames(untrans@exprs) <- make.names(untrans@parameters@data$desc, unique = TRUE)
    colnames(trans@exprs) <- make.names(trans@parameters@data$desc, unique = TRUE)

    channels <- colnames(untrans@exprs)
    channels <- channels[!grepl("NA", channels) & !channels == ""]

    for (channel in channels) {

      message("Plotting....", channel)

      p1 <- ggplot() +
        geom_density(as.data.frame(untrans@exprs), mapping = aes(x = !!sym(channel)), color = "red", alpha = 0.2) +
        ggtitle(paste0(channel, " Before transformation"), subtitle = tools::file_path_sans_ext(basename(fcs_files[random_file]))) +
        theme_classic()

      p2 <- ggplot() +
        geom_density(as.data.frame(trans@exprs), mapping = aes(x = !!sym(channel)), color = "black", alpha = 0.2) +
        ggtitle(paste0(channel, " After transformation"), subtitle = tools::file_path_sans_ext(basename(trans_files[random_file]))) +
        theme_classic()

      p3 <- p1 + p2

      ggsave(paste0(channel, "_histogram.png"), plot = p3, device = "png",
             path = file.path(output_folder, "ArcsinTrans", "Before and after histograms"), width = 10, height=5)

    }

  }

}
