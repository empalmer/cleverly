% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algorithm1.R
\name{get_clusters}
\alias{get_clusters}
\title{Get Cluster Membership from Pairwise Differences}
\usage{
get_clusters(v, K, P)
}
\arguments{
\item{v}{A matrix of final pairwise differences from Algorithm 1, used to infer which subjects belong to the same cluster.}

\item{K}{The number of response variables.}

\item{P}{The number of B-spline coefficients (equal to \code{order + nknots}).}
}
\value{
A list containing:
\describe{
\item{clusters}{A vector of cluster membership assignments.}
\item{indices}{A list of indices corresponding to each cluster.}
\item{n_clusters}{The total number of unique clusters identified.}
}
}
\description{
Derives cluster assignments based on the final pairwise difference estimates (\code{v}) from Algorithm 1.
}
