# Partials -----------------------------------------------------

#' Get partials matrix for i, j, l
#'
#' @param i subject index
#' @param j time index
#' @param l external variable index
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param i_index
#' @param Y0
#' @param alpha_ij
#'
#' @returns Matrix of dimension PK x K
#' @export
#'
get_partials_ijl <- function(i,
                             j,
                             l,
                             mi_vec,
                             i_index,
                             Y0,
                             Z,
                             B,
                             beta,
                             alpha_ij){
  Y_ij0 <- get_Y_ij0(i = i,
                     j = j,
                     Y0 = Y0,
                     i_index = i_index)
  U_ij <- get_U_ij(alpha_ij = alpha_ij)
  Z_ijl <- get_Z_ijl(i = i,
                     j = j,
                     l = l,
                     Z = Z,
                     i_index = i_index)
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
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Y0
#' @param alpha_i
#' @param i_index
#' @param P
#' @param K
#'
#' @returns Matrix of dimension KP x Kmi
#' @export
#'
get_partials_il <- function(i,
                            l,
                            Y0,
                            Z,
                            B,
                            beta,
                            alpha_i,
                            mi_vec,
                            i_index,
                            P,
                            K){
  mi <- mi_vec[i]
  partials_il <- matrix(nrow = K*P, ncol = K*mi)

  for (j in 1:mi) {
    alpha_ij <- alpha_i[[j]]
    partials_il[, ((j - 1)*K + 1):(j*K)] <- get_partials_ijl(i = i,
                                                             j = j,
                                                             l = l,
                                                             i_index = i_index,
                                                             mi_vec = mi_vec,
                                                             Y0 = Y0,
                                                             Z = Z,
                                                             B = B,
                                                             beta = beta,
                                                             alpha_ij = alpha_ij)
  }
  return(partials_il)
}



#' Partials list (of matrices) with respect to given l
#'
#' @param l
#' @param mi_vec
#' @param i_index
#' @param beta
#' @param Z
#' @param B
#' @param Y0
#' @param alpha
#' @param P
#' @param K
#'
#' @returns
#' @export
get_partials_l_list <- function(Y0,
                                l,
                                mi_vec,
                                i_index,
                                beta,
                                alpha,
                                Z,
                                B,
                                P,
                                K){
  n <- length(mi_vec)
  partials_list <- list()


  for (i in 1:n) {
    partials_list[[i]] <- get_partials_il(i = i,
                                          l = l,
                                          Y0 = Y0,
                                          Z = Z,
                                          B = B,
                                          beta = beta,
                                          alpha_i = alpha[[i]],
                                          mi_vec = mi_vec,
                                          i_index = i_index,
                                          P = P,
                                          K = K)
  }

  return(partials_list)
}
