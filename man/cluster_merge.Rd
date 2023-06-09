% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_merge.R
\name{cluster_merge}
\alias{cluster_merge}
\title{Merge clusters of single-cell cytometry data}
\usage{
cluster_merge(
  SE,
  clusters,
  reassign = TRUE,
  output_folder = NULL,
  SE_name = NULL
)
}
\arguments{
\item{SE}{SummarizedExperiment object or path to SummarizedExperiment object.}

\item{clusters}{Vector of the cluster IDs that needs to be merged. All th clusters will be merged in the first cluster in the vector.
c(5,8,23), cluster 8 and 23 will be merged in cluster 5.}

\item{reassign}{Logical. If TRUE, the cluster IDs are reassigned new IDs starting from 1.}

\item{output_folder}{Path of output directory to receive the subsetted SummarizedExperiment.}

\item{SE_name}{Name of the new subsetted SummarizedExperiment.}
}
\value{
SummarizedExperiment object with the chosen clusters merged.
}
\description{
Merge clusters of single-cell cytometry data
}
\details{
Occasionally, analyst might find that some clusters are over-clustered and wish to merge these clusters togther.
}
