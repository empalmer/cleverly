% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{get_Y_ij0}
\alias{get_Y_ij0}
\title{Get \eqn{Y_{ij0}}}
\usage{
get_Y_ij0(i, j, Y, i_index)
}
\arguments{
\item{i}{subject index}

\item{j}{time index}

\item{Y}{Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.}

\item{mi_vec}{vector of the number of timepoints for each sample. Of length n}
}
\value{
Scalar of the total sum constraint for a given i, j
}
\description{
Get \eqn{Y_{ij0}}
}
