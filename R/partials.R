# Partials -----------------------------------------------------

#' Get partials matrix for i, j, l
#'
#' @param i subject index
#' @param j time index
#' @param l external variable index
#' @param K Number of responses
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns Matrix of dimension PK x K
#' @export
#'
get_partials_ijl <- function(i, j, l, mi_vec, i_index, K, Y, Z, B, beta){
  Y_ij0 <- get_Y_ij0(i = i,
                     j = j,
                     Y = Y,
                     i_index = i_index)
  U_ij <- get_U_ij(
    i = i,
    j = j,
    beta = beta,
    Z = Z,
    B = B,
    K = K,
    i_index = i_index
  )
  Z_ijl <- get_Z_ijl(
    i = i,
    j = j,
    l = l,
    Z = Z,
    i_index = i_index
  )
  B_ij <- get_B_ij(i = i,
                   j = j,
                   B = B,
                   i_index = i_index)

  partials_ijl <- Y_ij0 * Z_ijl * kronecker(U_ij, B_ij)
  return(partials_ijl)
}



#' Get full partials matrix for a given i, l
#'
#' @param i subject index
#' @param l external variable index
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#'
#' @returns Matrix of dimension KP x Kmi
#' @export
#'
get_partials_il <- function(i, l, Y, Z, B, beta, mi_vec, i_index){
  P <- ncol(B)
  K <- ncol(Y)
  mi <- mi_vec[i]
  partials_il <- matrix(nrow = K*P, ncol = K*mi)
  for (j in 1:mi) {
    partials_il[, ((j - 1)*K + 1):(j*K)] <-
      get_partials_ijl(i = i,
                       j = j,
                       l = l,
                       i_index = i_index,
                       mi_vec = mi_vec,
                       K = K,
                       Y = Y,
                       Z = Z,
                       B = B,
                       beta = beta)
  }
  return(partials_il)
}



#' Partials list (of matrices) with respect to given l
#'
#' @param Y
#' @param l
#' @param mi_vec
#' @param i_index
#' @param beta
#' @param Z
#' @param B
#'
#' @returns
#' @export
get_partials_l_list <- function(Y, l, mi_vec, i_index, beta, Z, B){
  M <- nrow(Y)
  L <- ncol(Z) - 1
  K <- ncol(Y)
  n <- length(mi_vec)
  partials_list <- list()

  for (i in 1:n) {
    partials_list[[i]] <- get_partials_il(i = i,
                                          l = l,
                                          Y = Y,
                                          Z = Z,
                                          B = B,
                                          beta = beta,
                                          mi_vec = mi_vec,
                                          i_index = i_index)
  }

  return(partials_list)
}
