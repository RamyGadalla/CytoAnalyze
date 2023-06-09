% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_SE.R
\name{build_SE}
\alias{build_SE}
\title{Build SummarizedExperiment object}
\usage{
build_SE(
  input_folder,
  experiment_info,
  marker_info,
  output_folder = ".",
  export = FALSE,
  arcsinh_cofactor = NULL
)
}
\arguments{
\item{input_folder}{Path to folder containing raw fcs files.}

\item{experiment_info}{Path to experiment_info table; one of the csv files generated by \code{\link{metadata_tables}}.}

\item{marker_info}{Path to marker_info table; one of the csv files generated by \code{\link{metadata_tables}}.}

\item{output_folder}{Path to folder to recieve the SummarizedExperiment object.}

\item{export}{If TRUE, the SummarizedExperiment object gets saved in the outpput_folder location.}

\item{arcsinh_cofactor}{Cofactor value for the data arcsinh transformation. Default value is 5, which is the typically recommended data for mass cytometry.
Recommended value for flow cytometry is 150.}
}
\value{
SummarizedExperiment object
}
\description{
Converts a set of fcs files into SummarizedExperiment object for further analysis.
}
\seealso{
\code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}} \code{\link[diffcyt]{prepareData}} \code{\link[flowCore]{read.flowSet}}
}
