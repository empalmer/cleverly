# Nuisance parameters ----------------------------------------------------

#' Phi
#'
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param beta beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of responses
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param Y0
#' @param alpha
#' @param i_index
#' @param M
#'
#' @returns scalar phi
#' @export
get_phi <- function(Y, Y0, beta, alpha, Z, B, K, mi_vec, i_index, M){

  r <- get_pearson_residuals(Y = Y,
                             Y0 = Y0,
                             beta = beta,
                             alpha = alpha,
                             Z = Z,
                             B = B,
                             K = K,
                             mi_vec = mi_vec,
                             i_index = i_index,
                             M = M)
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
#' @param Y0
#' @param alpha
#' @param i_index
#' @param M
#'
#' @returns Pearson residuals vector for all i, j, k
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
  r <- c()
  for (i in 1:n) {
    ri <- get_pearson_residual_i(Y = Y,
                                 Y0 = Y0,
                                 i = i,
                                 beta = beta,
                                 alpha = alpha,
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
#' @param i subject index
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of responses
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param Y0
#' @param alpha
#' @param i_index
#'
#' @returns pearson residual vector for each i
#' @export
get_pearson_residual_i <- function(Y,
                                   Y0,
                                   i,
                                   beta,
                                   alpha = alpha,
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
    rij <- Yij_minus_muij /
      sqrt( Y_ij0 *
              (Y_ij0 + alpha_ij0) / (1 + alpha_ij0) *
              alpha_ij/alpha_ij0 * (1 - alpha_ij/alpha_ij0))
    ri[((j - 1)*K + 1):(K*j)] <- rij
  }

  return(ri)
}
