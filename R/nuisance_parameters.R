# Nuisance parameters ----------------------------------------------------

#' Phi
#'
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param beta beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of responses
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns scalar phi
#' @export
get_phi <- function(Y, beta, Z, B, K, mi_vec, i_index){

  r <- get_pearson_residuals(Y = Y,
                             beta = beta,
                             Z = Z,
                             B = B,
                             K = K,
                             mi_vec = mi_vec,
                             i_index = i_index)
  M <- sum(mi_vec)
  phi <- sum(r^2) / (K*M - 1)

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
#'
#' @returns Pearson residuals vector for all i, j, k
#' @export
get_pearson_residuals <- function(Y, beta, Z, B, K, mi_vec, i_index){
  r <- c()
  for (i in 1:length(mi_vec)) {
    ri <- get_pearson_residual_i(Y = Y,
                                 i = i,
                                 beta = beta,
                                 Z = Z,
                                 B = B,
                                 K = K,
                                 mi_vec = mi_vec,
                                 i_index = i_index)
    r <- c(r, ri)
  }
  return(r)
}

#' Get the pearson residual for a given i
#'
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param phi Current value of overdispersion parameter
#' @param i subject index
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of responses
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns pearson residual vector for each i
#' @export
get_pearson_residual_i <- function(Y,
                                   i,
                                   beta,
                                   Z,
                                   B,
                                   K,
                                   mi_vec,
                                   i_index){

  Yi_minus_mui <- get_Yi_minus_mui(i = i,
                                   Y = Y,
                                   mi_vec = mi_vec,
                                   i_index = i_index,
                                   beta = beta,
                                   Z = Z,
                                   B = B,
                                   K = K)
  ri <- numeric(mi_vec[i]*K)
  for (j in 1:mi_vec[i]) {
    Yij_minus_muij <- Yi_minus_mui[((j - 1)*K + 1):(K*j)]
    alpha_ij <- get_alpha_ij(beta = beta,
                             i = i,
                             j = j,
                             Z = Z,
                             B = B,
                             K = K,
                             i_index = i_index)
    Y_ij0 <- get_Y_ij0(i = i,
                       j = j,
                       Y = Y,
                       i_index = i_index)
    alpha_ij0 <- sum(alpha_ij)
    # Should be of length K.
    rij <- Yij_minus_muij /
      sqrt( Y_ij0 *
              (Y_ij0 + alpha_ij0) / (1 + alpha_ij0) *
              alpha_ij/alpha_ij0 * (1 - alpha_ij/alpha_ij0))
    ri[((j - 1)*K + 1):(K*j)] <- rij
  }

  return(ri)
}
