% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algorithm1.R
\name{algorithm1}
\alias{algorithm1}
\title{Algorithm 1}
\usage{
algorithm1(
  Y,
  Z,
  time,
  mi_vec,
  lp,
  gammas,
  psi,
  phi,
  tau,
  theta,
  C,
  d,
  nknots,
  order,
  tol,
  max_outer_iter,
  max_admm_iter
)
}
\arguments{
\item{Y}{Matrix of counts Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.}

\item{Z}{matrix that starts with a column of 1s. Of dimension M x L + 1}

\item{mi_vec}{count vector of length n for the number of time points for each i}

\item{lp}{clustering index}

\item{d}{Order for the difference matrix}

\item{nknots}{Number of knots for the B-spline basis}

\item{order}{Order of the B-spline basis}

\item{tol}{Tolerance for convergence}

\item{max_admm_iter}{}
}
\description{
Algorithm 1
}
