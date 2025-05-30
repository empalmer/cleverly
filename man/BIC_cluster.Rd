% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BIC.R
\name{BIC_cluster}
\alias{BIC_cluster}
\title{Compute BIC for Clustered Functional Data}
\usage{
BIC_cluster(y_ra_df, K, L, n_clusters, mi_vec, nknots, order)
}
\arguments{
\item{y_ra_df}{A data frame containing columns \code{yhat} (fitted values) and \code{y} (observed values), both assumed to be in relative abundance form.}

\item{K}{Number of responses}

\item{L}{Number of external variables}

\item{n_clusters}{Number identified clusters}

\item{mi_vec}{Vector of number of timepoints for each subject}

\item{nknots}{The number of internal knots used in the B-spline basis.}

\item{order}{The order of the B-spline basis.}
}
\value{
A single numeric value representing the modified BIC for the fitted model.
}
\description{
Calculates a modified Bayesian Information Criterion (BIC) for a fitted clustering model on functional data.
The BIC is computed based on the log-transformed relative abundance data and includes a penalty term proportional to the number of clusters and B-spline coefficients.
}
