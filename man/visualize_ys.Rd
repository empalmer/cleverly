% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization_helpers.R
\name{visualize_ys}
\alias{visualize_ys}
\title{Visualize just the Ys.}
\usage{
visualize_ys(Y, is, time, beta)
}
\arguments{
\item{Y}{Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.}

\item{is}{numeric of length M of individual ids}

\item{time}{vector of time values for each subject/time}

\item{beta}{matrix of beta (or beta hat) of dimension (P*K) x L}
}
\value{
ggplot object
}
\description{
Visualize just the Ys.
}
