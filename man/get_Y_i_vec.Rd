% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{get_Y_i_vec}
\alias{get_Y_i_vec}
\title{Get Yi as a vector}
\usage{
get_Y_i_vec(i, mi_vec, Y, Y_mat)
}
\arguments{
\item{i}{subject index}

\item{mi_vec}{vector of the number of timepoints for each sample. Of length n}

\item{Y}{Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.}

\item{Y_mat}{(optional) pre-calculated Yi-matrix from Y}
}
\value{
Vector of length K * mi
}
\description{
different k responses are listed after each other. so it goes
y(i,j=1,k=1), ... yi(j = mi, k = 1), yi(j = 1, k = 2)...
NOPE! Other way around...
}
