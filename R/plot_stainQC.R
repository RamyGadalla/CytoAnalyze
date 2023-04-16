#' Plot histograms for visual inspection of batch controls for single-cell cytometry data
#'
#'
#' @description Plot histograms of data channels in fcs files for visual inspection and comparison between batch controls.
#'
#' @param input_folder       Path to folder containing fcs files.
#' @param pattern            A common string to look for among fcs files in the \code{input_folder} e.g. "control" , "healthy_donor".etc.
#'                           It could also be a regular expression.
#' @param output_folder      Path to folder to save the transformed fcs files.
#'
#' @return
#' Histogram plots of all channels in fcs files in all control files.
#'
#' @importFrom ggridges geom_density_ridges2
#' @import ggplot2
#' @importFrom patchwork wrap_plots
#'
#' @importFrom stringi stri_replace_all_regex
#' @importFrom Biobase sampleNames
#'
#' @export


plot_stainQC <- function (input_folder,
                          pattern,
                          output_folder) {

  fcs_files <- list.files(path = input_folder, pattern = paste0(".*", pattern, ".*fcs$.*"), full.names=TRUE, recursive = TRUE)

  print(fcs_files)

  all_controls <- flowCore::read.flowSet(fcs_files, transformation = FALSE, truncate_max_range = FALSE)



# Making sure all marker names are identical in fcs files.

  flag <- logical()

  for (i in 2:length(all_controls)) {
    flag[1] <- c(FALSE)

    if (!identical(all_controls[[1]]@parameters@data$desc, all_controls[[i]]@parameters@data$desc)) {
      cat("\033[31mParamter names in",sampleNames(all_controls)[i], "do not match those in", sampleNames(all_controls)[1], "\n\033[0m")
      flag[i] <- TRUE
    } else {
      flag[i] <- FALSE
    }
  }

  if (any(flag)) {
    stop("can't move forward with paramter names not matching.")
  }

  # clean marker names
  for (i in 1:length(all_controls)) {

    all_controls[[i]]@parameters$desc <-  stringi::stri_replace_all_regex (all_controls[[i]]@parameters$desc,
                                                                           pattern =c("^(.*?)_", "EQ", ":", "::", "-"),
                                                                           replacement = c(""),
                                                                           vectorize_all = FALSE
    )

  }


  channels <- all_controls[[1]]@parameters$desc
  channels <- channels[!is.na(channels) & !channels == ""]
  channels <- make.names(channels, unique = TRUE)

  dir.create(file.path(output_folder, "QC plots"), showWarnings = FALSE)


  all_plots <- list()

  message("\n.\n.Plotting histograms")

  for (channel in channels) {

    df <- data.frame()

    for (i in 1:length(all_controls)) {

      colnames(all_controls[[i]]@exprs) <- make.names(all_controls[[i]]@parameters@data$desc, unique=TRUE)
      df_temp <- cbind(as.data.frame(all_controls[[i]]@exprs[,channel]), rep(sampleNames(all_controls)[i], length(all_controls[[i]]@exprs[,channel])))
      colnames(df_temp) <- c(channel, "sample_name")

      df <-  rbind(df, df_temp)

    }

    message("Plotting...", channel)

    p <- ggplot(df, aes(x = !!sym(channel), y= .data$sample_name, fill= .data$sample_name)) +
      geom_density_ridges2(scale = 4, alpha = 0.5) +
      ggtitle(channel) +
      theme_classic() + theme(axis.text.y=element_text(size=12,face="bold"),
                              axis.title.y = element_text(size = 14, face="bold"))

    ggsave(paste0(channel, "staining_QC", ".pdf"), plot = p, device= "pdf", path= file.path(output_folder, "QC plots"), width= 12)

    all_plots[[channel]] <- p
  }

  #wrap plot
  p_all <- wrap_plots(all_plots, ncol = 1)
  #ggsave("All_staining_QC.pdf", plot = p_all, device= "pdf", path= file.path(output_folder, "QC plots"))

}

