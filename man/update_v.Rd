% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algorithm3.R
\name{update_v}
\alias{update_v}
\title{Update v}
\usage{
update_v(beta, lp, lambda, Kappa, P, gammas, tau, theta, psi)
}
\arguments{
\item{beta}{matrix of beta (or beta hat) of dimension (P*K) x L}

\item{lp}{clustering index (integer between 0 and L)}

\item{lambda}{ADMM vector of dimension (P*|Kappa|) x1}

\item{Kappa}{keeps track of each possible response pair}

\item{P}{Number of B-spline coefficients (order + nknots)}

\item{gammas}{Vector of dimension L + 1 for penalizing the D matrix}

\item{tau}{MCP hyper parameter.}

\item{theta}{ADMM hyper parameter.}

\item{psi}{Hyperparameter for clustering penalty (larger drives pairwise differences to zero)}
}
\value{
Vector of updated vs for ADMM step
}
\description{
Update v step of ADMM algorithm. V is the pairwise differences vector of dimension |Kappa|P
}
