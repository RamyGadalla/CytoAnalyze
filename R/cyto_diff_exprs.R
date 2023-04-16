#' Differential expression analysis of biomarkers in clustered single-cell cytometry data
#'
#'
#'
#' @description  Performs GLM modeling to clustered single-cell cytometry data. Treats selected biomarkers as explanatory variables and the experimental conditions as response variable.
#'
#'
#'
#' @param SE                \code{SummarizedExperiment} object or path to \code{SummarizedExperiment} object.
#' @param groups            Experimental groups to generate the plots for. It has to be a variable in the data i.e. in \code{rowdata()}. Levels of this variable will be compared.
#' @param assay_name        A character vector of length 1 indicates the assay to be used in plotting the UMAP. By default, it is \code{"exprs"}.
#' @param exp_subject       A character vector of length 1 indicates the experimental subject e.g patient_id, sample_id..etc.
#' @param count_threshold   An Integer (default is 100). Clusters with cell count that are equal or less than this value will be excluded from the analysis.
#' @param output_folder     Path to folder to receive pdf file of the plots and tables.
#' @param num_boot          Numeric. Indicates the number of bootstrapping runs for the bootstrapped generalized linear model. Only used for the unpaired experimental design.
#' @param paired            logical (default is FALSE). Indicate whether the experiment design has a pairing variable. If TRUE, the exp_subject will be used as the pairing variable.
#'
#'
#'
#' @importFrom CytoGLMM cytoglm
#' @importFrom CytoGLMM cytoglmm
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr left_join
#' @importFrom dplyr join_by
#' @import ggplot2
#' @import SummarizedExperiment
#' @importFrom magrittr %<>%
#' @importFrom limma makeContrasts
#' @importFrom stats model.matrix
#' @importFrom tibble as_tibble
#' @importFrom tibble rownames_to_column
#' @importFrom stats median
#' @importFrom stats formula
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @importFrom patchwork wrap_plots
#'
#'
#' @seealso \href{https://github.com/ChristofSeiler/CytoGLMM}{CytoGLMM}.
#'
#' @return
#' Forest plot for marker of interest on each cluster available in the data.
#' Table of coefficients and significance values.
#'
#'
#' @export
#'
#'
#'
cyto_diff_exprs <- function(SE,
                            groups = "groups",
                            assay_name = "exprs",
                            exp_subject = "sample_id",
                            count_threshold = 100,
                            output_folder = ".",
                            num_boot = 1000,
                            paired = FALSE) {

  if (!dir.exists(output_folder)) {
    stop("output_folder does not exist.")
  }


  if (is.character(SE)) {
    SE <- readRDS(SE)
  }

  if (!is(SE, "SummarizedExperiment")) {
    stop('SE needs to be a SummarizedExperiment object or path to SummarizedExperiment object')
  }

  if (is.character(groups)) {
    if(!groups %in% colnames(rowData(SE))) {
      stop('Choosen "groups" is not in the data.')
    }
  }

  if (is.character(exp_subject)) {
    if(!exp_subject %in% colnames(rowData(SE))) {
      stop('Choosen "exp_subject" is not in the data.')
    }
  }

  #CytoGLMM

  message("Differential expression CytoGLMM...")


  if(is.null(metadata(SE))) {
    stop("metadata(SE) can not be null")
  }

  rownames(colData(SE)) <- colData(SE)$clean_name

  df <- as.data.frame(cbind(rowData(SE), assay(SE, assay_name)), optional = FALSE)

  if(is.null(metadata(SE))) {
    stop("metadata(SE) can not be null")
  }

  if (!paired) {

    if(any(!isUnique(metadata(SE)$experiment_info[,exp_subject]))) {
      stop("\033[31mSome of the exp_subject values provided are paired, but the argument paired is set to FALSE.\nRerun cyto_diff_exprs() with paired = FALSE or provide a paired exp_subject.\033[0m")
    }

    i <- 0
    plot_list <- list()
    coeff_table <- data.frame()

    for (cluster in levels(df$cluster_id)) {

      if (nrow(df[df$cluster_id == cluster,]) <= count_threshold) {

        i = i+ 1

      } else if (nrow(df[df$cluster_id == cluster,]) > count_threshold) {

        glm_fit <- CytoGLMM::cytoglm(df[df$cluster_id == cluster,],
                                     protein_names = make.names(colData(SE)[colData(SE)$marker_class =="state","clean_name"]),
                                     condition = groups,
                                     group = exp_subject,
                                     num_cores = parallel::detectCores()-1,
                                     num_boot = num_boot)

        coeff_pvals <- summary(glm_fit)

        if (any(is.na(coeff_pvals))) {
          message("\033[31mWarning: ", paste(coeff_pvals$protein_name[which(is.na(coeff_pvals$pvalues_unadj))], collapse=" "), " on cluster ", cluster, " is producing NA value. This cluster will be skipped. Rerun cyto_diff_exprs() without this marker(s).\033[0m")
          next
        }

        df_for_plotting <- glm_fit$tb_coef %>%
          group_by(.data$protein_name) %>%
          summarize(
            coeff_median = median(.data$coeff),
            min_median = quantile(.data$coeff, probs = 0.05),
            max_median = quantile(.data$coeff, probs = 0.95)
          ) %>%
          left_join(coeff_pvals, by = join_by("protein_name"))

        df_for_plotting$cluster_id <- as.factor(cluster)

        coeff_table <- rbind(coeff_table, df_for_plotting)

        p <- ggplot(df_for_plotting,
                    aes(x = .data$coeff_median, y = .data$protein_name)) +
          geom_vline(xintercept = 0, color = "red") +
          geom_point(size = 2) +
          geom_errorbarh(aes(xmin = .data$min_median, xmax = .data$max_median), color="#055C9D") +
          ggtitle(paste("Cluster", cluster, "(count:", nrow(df[df$cluster_id == cluster,]), ")")) +
          xlab("Coefficient median") + ylab("Biomarker") +
          theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                             panel.background = element_rect(fill = "#efefef"),
                             axis.title = element_text(size = 16, face = "bold"),
                             axis.text=element_text(size=12, face="bold"),
                             plot.title = element_text(size = 15, face = "bold")) +
          geom_text(aes(label = ifelse(.data$pvalues_adj < 0.05, "*", "")), size=7, nudge_x = 0.1, nudge_y = 0.1)

        dir.create(file.path(output_folder,"stats", "Differential expression", "forest_plots_unpaired"), recursive = TRUE, showWarnings = FALSE)

        ggsave(paste0("cluster_", cluster, "_forestplot_glm.pdf"), plot = p, device = "pdf",
               path = file.path(output_folder, "stats", "Differential expression", "forest_plots_unpaired"), height = 8, width = 12)

        plot_list[[cluster]] <- p

      }
    }

    message("The GLM formula used is:\n", glm_fit$formula_str)

    message("\033[31mNote: ", i, " clusters were lower than the count threshold and excluded out of the analysis.\033[0m")

    write.csv(coeff_table, row.names = FALSE, file.path(output_folder,"stats", "Differential expression", "coeff_table_unpaired.csv"))

    p_all <- patchwork::wrap_plots(plot_list)
    ggsave(paste0("All_clusters_forestplot_glm.pdf"), plot =  p_all, device = "pdf",
           path = file.path(output_folder, "stats", "Differential expression", "forest_plots_unpaired"), height = 25, width = 30)

  }


  if (paired) {

    # If there is no pairing in the variable provided.
    if(!any(!isUnique(metadata(SE)$experiment_info[,exp_subject]))){
      stop("\033[31mThe values in exp_subject provided all are unique values. There is no pairing, but the argument paired is set to TRUE.\nRerun cyto_diff_exprs() with paired = FALSE or provide a paired exp_subject.\033[0m")
    }

    message("Paired design adjusting for random effect: ", exp_subject, " ...")

    i <- 0
    plot_list <- list()
    coeff_table <- data.frame()

    for (cluster in levels(df$cluster_id)) {

      if (nrow(df[df$cluster_id == cluster,]) <= count_threshold) {

        i = i+ 1

      } else if (nrow(df[df$cluster_id == cluster,]) > count_threshold) {

        glmm_fit <- CytoGLMM::cytoglmm(df[df$cluster_id == cluster,],
                                       protein_names = make.names(colData(SE)[colData(SE)$marker_class =="state","clean_name"]),
                                       condition = groups ,
                                       group = exp_subject,
                                       num_cores = parallel::detectCores()-1)

        ## this part copied from plot.cytoglmm() for visualization
        alpha <- 0.05
        ci <- qnorm(1-alpha/2)

        summ <- summary(glmm_fit$glmmfit)
        coeff_pvals <- summary(glmm_fit)

        if (any(is.na(coeff_pvals))) {
          stop("\n\033[31mone or more cluster has a very low cell count. Try increasing the count threshold.\033[0m")
        }

        df_for_plotting <- summ$coefficient
        df_for_plotting<- df_for_plotting[-1,]
        df_for_plotting %<>%
          as.data.frame %>%
          rownames_to_column(var = "protein_name") %>%
          as_tibble
        df_for_plotting %<>%
          mutate(high = df_for_plotting$Estimate+ci*df_for_plotting$`Std. Error`,
                 low = df_for_plotting$Estimate-ci*df_for_plotting$`Std. Error`)
        df_for_plotting <- df_for_plotting[df_for_plotting$protein_name %in% glmm_fit$protein_names,]
        df_for_plotting <- left_join(df_for_plotting, coeff_pvals, by = join_by("protein_name"))


        df_for_plotting$cluster_id <- as.factor(cluster)
        coeff_table <- rbind(coeff_table, df_for_plotting)

        p <- ggplot( df_for_plotting, aes(x = .data$Estimate, y = .data$protein_name)) +
          geom_vline(xintercept = 0,color = "red") +
          geom_point(size = 2) +
          geom_errorbarh(aes(xmin = .data$low, xmax = .data$high), color="#055C9D") +
          ggtitle(paste("Cluster", cluster, "(count:", nrow(df[df$cluster_id == cluster,]), ")")) +
          xlab("Coefficient estimate") + ylab("Biomarker") +
          theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                             panel.background = element_rect(fill = "#efefef"),
                             axis.title = element_text(size = 16, face = "bold"),
                             axis.text=element_text(size=12, face="bold"),
                             plot.title = element_text(size = 15, face = "bold")) +
          geom_text(aes(label = ifelse(.data$pvalues_adj < 0.05, "*", "")), size=7, nudge_x = 0.1, nudge_y = 0.1)



        dir.create(file.path(output_folder,"stats", "Differential expression", "forest_plots_paired"), recursive = TRUE, showWarnings = FALSE)

        ggsave(paste0("cluster_", cluster, "_forestplot_glm.pdf"), plot = p, device = "pdf",
               path = file.path(output_folder, "stats", "Differential expression", "forest_plots_paired"), height = 8, width = 12)

        plot_list[[cluster]] <- p

      }
    }

    message("\033[31mNote: ", i, " clusters were lower than the count threshold and excluded out of the analysis.\033[0m")

    write.csv(coeff_table, row.names = FALSE, file.path(output_folder,"stats", "Differential expression", "coeff_table_paired.csv"))

    p_all <- patchwork::wrap_plots(plot_list)
    ggsave(paste0("All_clusters_forestplot_glm.pdf"), plot =  p_all, device = "pdf",
           path = file.path(output_folder, "stats", "Differential expression", "forest_plots_paired"), height = 25, width = 30)

  }

  message("Done CytoGLMM.")
}



