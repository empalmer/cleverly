# Nuisance parameters ----------------------------------------------------

#' Phi
#'
#' @param K Number of responses
#' @param M Number of samples times timepoints for each sample
#' @param pearson_residuals list of pearson residuals for each i
#' @param L Number of external variables
#' @param P Number of bspline basis parameters (order + nknots)
#'
#' @returns scalar phi
#' @export
get_phi <- function(pearson_residuals, K, M, L, P){

  phi <- sum(unlist(pearson_residuals)^2) / (K * M - ((L + 1) * P))
  if (phi > 1e8) {
    warning("Phi is too large; setting to 1e6.")
    print("Phi is too large; setting to 1e6.")
    phi <- 1e8
  }

  return(phi)
}


#' Get pearson residuals
#'
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param beta beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of responses
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param Y0 Vector of total count for each sample
#' @param alpha list of alpha that can be subsetted by i and j
#' @param i_index starting index of the ith subject in the data
#' @param M Number of samples times timepoints for each sample
#'
#' @returns List of Pearson residuals for all i, j, k, indexed by i
#' @export
get_pearson_residuals <- function(Y,
                                  Y0,
                                  beta,
                                  alpha,
                                  Z,
                                  B,
                                  K,
                                  mi_vec,
                                  i_index,
                                  M){
  n <- length(mi_vec)
  r_list <- c()
  for (i in 1:n) {
    r_list[[i]] <- get_pearson_residual_i(Y = Y,
                                          Y0 = Y0,
                                          i = i,
                                          beta = beta,
                                          alpha = alpha,
                                          Z = Z,
                                          B = B,
                                          K = K,
                                          mi_vec = mi_vec,
                                          i_index = i_index)
  }
  return(r_list)
}

#' Get the pearson residual for a given i
#'
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param i subject index
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of responses
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param Y0 Vector of total count for each sample
#' @param alpha list of alpha that can be subsetted by i and j
#' @param i_index starting index of the ith subject in the data
#'
#' @returns pearson residual vector for each i
#' @export
get_pearson_residual_i <- function(Y,
                                   Y0,
                                   i,
                                   beta,
                                   alpha,
                                   Z,
                                   B,
                                   K,
                                   mi_vec,
                                   i_index){

  Yi_minus_mui <- get_Yi_minus_mui(i = i,
                                   Y = Y,
                                   Y0 = Y0,
                                   mi_vec = mi_vec,
                                   i_index = i_index,
                                   alpha = alpha,
                                   K = K)
  ri <- numeric(mi_vec[i]*K)
  for (j in 1:mi_vec[i]) {
    Yij_minus_muij <- Yi_minus_mui[((j - 1)*K + 1):(K*j)]
    alpha_ij <- alpha[[i]][[j]]
    Y_ij0 <- get_Y_ij0(i = i,
                       j = j,
                       Y0 = Y0,
                       i_index = i_index)
    alpha_ij0 <- sum(alpha_ij)
    # Should be of length K.
    denom <- sqrt( Y_ij0 *
                     (Y_ij0 + alpha_ij0) / (1 + alpha_ij0) *
                     alpha_ij / alpha_ij0 * (1 - alpha_ij / alpha_ij0) )

    too_small <- denom < 1e-5
    if (any(too_small, na.rm = TRUE)) {
      warning("Some denominator values in the pearson residuals are too small; replacing with 1e-3.")
      print("Some denominator values in the pearson residuals are too small; replacing with 1e-6.")
      denom[too_small] <- 1e-5
    }

    rij <- Yij_minus_muij / denom

    ri[((j - 1)*K + 1):(K*j)] <- rij

  }

  return(ri)
}
