% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{get_alpha_ijk}
\alias{get_alpha_ijk}
\title{Calculate alpha for a single \eqn{i = 1, j = 1, k = 1} from \eqn{\beta}}
\usage{
get_alpha_ijk(i, j, k, beta_ks, Z_ij, B_ij, i_index, P, L)
}
\arguments{
\item{i}{subject index}

\item{j}{time index}

\item{k}{response index}

\item{beta_ks}{list of beta for each k}

\item{Z_ij}{vector of length l.}

\item{B_ij}{B spline basis for i and j}

\item{i_index}{starting index of the ith subject in the data}

\item{P}{Number of B-spline coefficients (order + nknots)}

\item{L}{Number of external variables}
}
\value{
Numeric (1x1)
}
\description{
Calculate alpha for a single \eqn{i = 1, j = 1, k = 1} from \eqn{\beta}
}
