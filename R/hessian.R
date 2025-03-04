#' Get hessian elements
#'
#' Get for for the ith subject deriving with respect to the lth external variable
#'
#' @param i subject index
#' @param l external variable index
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param phi Current value of overdispersion parameter
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns Matrix of dimension KP times KP
#' @export
#'
get_Hessian_il <- function(i, l, Y, mi_vec, beta, Z, B, phi){
  # Third term:
  K <- ncol(Y)
  partials_il <- get_partials_il(i = i,
                                l = l,
                                Y = Y,
                                mi_vec = mi_vec,
                                beta = beta,
                                Z = Z,
                                B = B)

  Vi_inv <- get_Vi_inv(i = i,
                        Y = Y,
                        mi_vec = mi_vec,
                        phi = phi,
                        beta = beta,
                        Z = Z,
                        B = B,
                        K = K)


  hessian_il <- -partials_il %*% Vi_inv %*% t(partials_il)

  return(hessian_il)
}


#' Get hessian with respect to the lth external variable
#'
#' @param l external variable index
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param phi Current value of overdispersion parameter
#' @param C Constant for determining the hessian change.
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns Matrix of dimension KP times KP
#' @export
#'
get_Hessian_l <- function(l, Y, mi_vec, beta, Z, B, phi, C) {
  n <- length(mi_vec)


  hessian_l <- matrix(0, nrow = nrow(beta), ncol = nrow(beta))
  for (i in 1:n) {
    hessian_il <- get_Hessian_il(i = i,
                                 l = l,
                                 Y = Y,
                                 mi_vec = mi_vec,
                                 beta = beta,
                                 Z = Z,
                                 B = B,
                                 phi = phi)
    hessian_l <- hessian_l + hessian_il

  }


  # Diagonalize hessian:
  diag_hessian <- -pmax(diag(-hessian_l), C)

  C_hessian <- diag(diag_hessian)

  return(C_hessian)
  #return(hessian_l)
}
