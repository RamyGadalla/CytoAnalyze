

#' @title Phenograh clustering k parameter sweeps
#'
#' @description   Test multiple values of k on clustering a given data set. Other clustering parameters are defaulted to simulate Phenograph algorithm.
#'
#'
#'
#'
#'
#' @param SE                Path to SummarizedExperiment object created by \code{\link{build_SE}}
#' @param k                 Integer vector represent k values to be tested.
#' @param export            Logical. If TRUE (default), plot is saved in output_folder.
#' @param output_folder     Path to folder to receive the output.
#' @param assay_name        Single character vector of the name of the \code{assay()} slot to be used in clustering.
#' @param show.plot         Logical. If TRUE final plot is  saved and shown. Default (FALSE) will save the plot in working directory without showing it.
#'
#' @seealso  For more information and more metrics to assess clustering parameters; see \href{https://bioconductor.org/packages/devel/bioc/vignettes/bluster/inst/doc/diagnostics.html}{Assorted clustering diagnostics}
#' @import bluster BiocParallel ggplot2 diffcyt SummarizedExperiment
#' @return Silhouette plot for the k values selected.
#'
#' @note The performance of \code{k_sweep()} is affected significantly by the number of cores and memory available for computation, and also the size of the dataset.
#'
#' @export

k_sweep <- function (SE,
                     k,
                     export= TRUE,
                     output_folder =".",
                     assay_name = "exprs",
                     show.plot=FALSE) {


  if (is.character(SE)) {
    SE <- readRDS(SE)
  }


  if (is(SE, "SummarizedExperiment")) {

    message("k_sweep() progress bar....")
    sweep <- bluster::clusterSweep(assay(SE, assay_name)[,colData(SE)$marker_class == "type"],
                                   BLUSPARAM=bluster::NNGraphParam(),
                                   k=as.integer(k),
                                   type = "jaccard",
                                   cluster.fun = "louvain",
                                   BPPARAM = BiocParallel::MulticoreParam(progressbar = TRUE))

    colnames(sweep$clusters) <-  as.character(k)


    ## assess clustering via silhouette
    # In principle, we could choose the clustering with the greatest separation for further analysis; however, this tends to be disappointing as it often favors overly broad clusters

    sil <- vapply(as.list(sweep$clusters),
                  function(x) mean(bluster::approxSilhouette(SummarizedExperiment::assay(SE, assay_name),x)$width),0)

    #Plot silhouette values

    df <- data.frame(k= factor(names(sil),levels=k), sil)

    message("Generating the averge silhouette plot")
    p <- ggplot(df,
                aes(x=k , sil, group=1)) +
      geom_line(linetype="dashed") +
      geom_point(shape=19, color="red") +
      theme_classic(base_size = 15) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      xlab("k") +
      ylab("Average silhouette width")


    if (export) {
      suppressMessages(ggsave(file.path(output_folder, "k_sweep.pdf"), plot = p))
      message("Plot is saved at ", output_folder)
    }


    if (show.plot) {
      plot(p)
    }

  }

}

