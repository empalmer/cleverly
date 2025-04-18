% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algorithm3.R
\name{update_lambda}
\alias{update_lambda}
\title{Update lambda step of the admm algorithm}
\usage{
update_lambda(beta_lp, v, lambda, A, theta)
}
\arguments{
\item{beta_lp}{beta_lp is the lp column of beta}

\item{v}{differences for each pair of betas of dimension (P*|Kappa|) x1}

\item{lambda}{value of lambda from the previous iteration}

\item{A}{Matrix keeping track of pairs of responses to calculate pairwise differences}

\item{theta}{ADMM hyper parameter.}
}
\value{
vector of updated lambdas for the ADMM step
}
\description{
Update lambda step of the admm algorithm
}
