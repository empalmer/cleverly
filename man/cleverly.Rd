% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleverly.R
\name{cleverly}
\alias{cleverly}
\title{cleverly}
\usage{
cleverly(
  Y,
  Z,
  subject_ids,
  time,
  lp,
  gamma,
  d = 3,
  nknots = 3,
  order = 3,
  tol = 1e+06,
  smax = 100
)
}
\arguments{
\item{Y}{either a data frame or matrix of numeric response variables. Each response should be a separate column. Each row should be a separate subject/time combination. There should be M total rows. Must be ordered in time.}

\item{Z}{Matrix or data frame containing a column for each external variable. There should be M rows and L columns.}

\item{subject_ids}{either a vector of length(Y) or a column reference if Y is a data frame}

\item{time}{either a vector of length(Y) or a column reference if Y is a data frame. Must be numeric}

\item{lp}{either a numeric index of which external variable to cluster on, or the name of the column of Z that contains the clustering variable. Specify numeric 0 to cluster via baseline.}

\item{gamma}{Vector of length (L + 1) of penalization hyper parameters}

\item{d}{Order for the difference matrix}

\item{nknots}{Number of knots for the B-spline basis}

\item{order}{Order of the B-spline basis}

\item{tol}{Tolerance for convergence}

\item{smax}{Maximum number of iterations}
}
\description{
cleverly
}
