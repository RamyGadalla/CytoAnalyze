% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cyto_clust_abundance.R
\name{cyto_clust_abundance}
\alias{cyto_clust_abundance}
\title{Differential abundance test for cytometry data using edgeR}
\usage{
cyto_clust_abundance(
  SE,
  groups = "groups",
  contrast,
  exp_subject = NULL,
  paired = FALSE,
  output_folder = "."
)
}
\arguments{
\item{SE}{\code{SummarizedExperiment} object or path to \code{SummarizedExperiment} object.}

\item{groups}{Experimental groups to generate the plots for. It has to be a variable in the data i.e. in \code{rowdata()}. Levels of this variable will be compared.}

\item{contrast}{A character vector of length 2, specify two levels of "groups" variable to build contrast matrix. \code{(contrast[1]-constrast[2])}.}

\item{exp_subject}{Default is NULL. If the experimental design is paired, this variable will be used in the formula to build paired design matrix.}

\item{paired}{Logical (default is FALSE). If TRUE a paired experimental design is performed. The pairing is done on the exp_subject variable.}

\item{output_folder}{Path to folder to receive pdf file of the plots and tables.}
}
\value{
A table (csv file) with log fold change and statistical values for the comparison between the two contrasting groups using edgeR. Also, the function outputs a column plot of the log-fold changes of clusters between groups.
}
\description{
Performs differential abundance test for cell clusters of single-cell cytometry data using edgeR test.
}
\seealso{
\href{https://github.com/lmweber/diffcyt}{diffcyt}, \href{https://bioconductor.org/packages/release/bioc/html/edgeR.html}{edgeR}
}
