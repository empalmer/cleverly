#' cleverly:
#'
#' cleverly is an R package implementing a novel clustering method for longitudinal microbiome data that incorporates external variables into the modeling framework. Unlike traditional clustering approaches, which often ignore time structure and external influences, cleverly uses a B-spline smoothing approach combined with Dirichlet-multinomial generalized estimating equations (GEEs) to account for compositional and temporal correlation in sequencing data. This enables clustering of microbial taxa based on either baseline behavior or their response to external variables, such as environmental or physiological conditions.
#'
#' @section Mypackage functions:
#' The mypackage functions ...
#'
"_PACKAGE"
#' @name cleverly
#' @useDynLib cleverly, .registration = TRUE
## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
