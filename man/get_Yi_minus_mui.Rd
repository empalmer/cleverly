% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{get_Yi_minus_mui}
\alias{get_Yi_minus_mui}
\title{Get the Yi minus mui term of the gradient}
\usage{
get_Yi_minus_mui(i, Y, mi_vec, i_index, beta, Z, B, K)
}
\arguments{
\item{i}{subject index}

\item{Y}{Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.}

\item{mi_vec}{vector of the number of timepoints for each sample. Of length n}

\item{beta}{matrix of beta (or beta hat) of dimension (P*K) x L}

\item{Z}{Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.}

\item{B}{B spline basis matrix of dimension (N x P)}

\item{K}{Number of responses}
}
\value{
Vector of length Kmi
}
\description{
Get the Yi minus mui term of the gradient
}
