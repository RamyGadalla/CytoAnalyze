% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset_SE.R
\name{subset_SE}
\alias{subset_SE}
\title{Subsetting SummarizedExperiment for CytoAnalyze}
\usage{
subset_SE(
  SE,
  group,
  no_keep = NULL,
  keep = NULL,
  output_folder = NULL,
  SE_name = NULL,
  rm_na = TRUE
)
}
\arguments{
\item{SE}{\code{SummarizedExperiment} object or path to \code{SummarizedExperiment} object.}

\item{group}{Character. Column name in \code{rowData()} to subset from.}

\item{no_keep}{Character vector of the levels in \code{group} that need to be removed.}

\item{keep}{Character vector of the levels in \code{group} that need to be kept.}

\item{output_folder}{Path of output directory to recieve the subsetted SummarizedExperiment.}

\item{SE_name}{Name of the new subsetted SummarizedExperiment.}

\item{rm_na}{Remove NA values.}
}
\value{
SummarizedExperiment
}
\description{
Subsets SummarizedExperiment object for CytoAnalyze analysis pipeline.
}
