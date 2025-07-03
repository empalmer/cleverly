
# Calculate cluster BIC ---------------------------------------------------

#' Compute BIC for Clustered Functional Data
#'
#' Calculates a modified Bayesian Information Criterion (BIC) for a fitted clustering model on functional data.
#' The BIC is computed based on the log-transformed relative abundance data and includes a penalty term proportional to the number of clusters and B-spline coefficients.
#'
#' @param y_ra_df A data frame containing columns \code{yhat} (fitted values) and \code{y} (observed values), both assumed to be in relative abundance form.
#' @param K Number of responses
#' @param n_clusters Number identified clusters
#' @param mi_vec Vector of number of timepoints for each subject
#' @param nknots The number of internal knots used in the B-spline basis.
#' @param order The order of the B-spline basis.
#' @param L Number of external variables
#'
#' @returns A single numeric value representing the modified BIC for the fitted model.

BIC_cluster <- function(y_ra_df,
                        K,
                        L,
                        n_clusters,
                        mi_vec,
                        nknots,
                        order){

  N <- sum(mi_vec)
  yhat_y <- log(y_ra_df$yhat + .01) - log(y_ra_df$y + .01)
  first_term <- sum(yhat_y^2)/(N*K)
  second_term <- log(N * K) * n_clusters * (L + 1) * (order + nknots)/(N*K)

  BIC <- log(first_term) + second_term
  return(list(BIC = BIC,
              first_term = log(first_term),
              second_term = second_term))
}


#' BIC_cluster_group
#'
#' @param y_hat_counts Estimated y hat in counts
#' @param y_counts True Y counts
#' @param beta_group Re-fit group betas
#' @param clusters Clusters
#' @param M Total time points x samples
#' @param K Number of responses
#' @param nknots Number of knots
#' @param order Order of the B-spline basis
#'
#' @returns List of BIC, first term, and second term

BIC_cluster_group <- function(y_hat_counts,
                              y_counts,
                              beta_group,
                              clusters,
                              M,
                              K,
                              nknots,
                              order) {

  n_clusters <- clusters$no

  y_log <- log(y_counts + 1)
  y_hat_log <- log(y_hat_counts + 1)

  N <- M * K

  first_term <- log((1/N) * sum((y_log - y_hat_log)^2)) * N
  second_term <- log(N) * (order + nknots) * n_clusters

  BIC_group <- first_term + second_term

  return(list(BIC = BIC_group,
              first_term = first_term,
              second_term = second_term))
}


