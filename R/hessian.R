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
#' @param i_index
#' @param Vi_inv (optional) Vi inverse If supplied is faster
#' @param partials_il
#'
#' @returns Matrix of dimension KP times KP
#' @export
#'
get_dHessian_il <- function(i,
                            l,
                            Y,
                            mi_vec,
                            i_index,
                            beta,
                            Z,
                            B,
                            phi,
                            Vi_inv,
                            partials_il){
  # Third term:
  K <- ncol(Y)
  partials_il <- get_partials_il(i = i,
                                l = l,
                                Y = Y,
                                mi_vec = mi_vec,
                                i_index = i_index,
                                beta = beta,
                                Z = Z,
                                B = B)

  if (missing(Vi_inv)) {
    stop("Missing Vi inverse argument")
    # Vi_inv <- get_Vi_inv(i = i,
    #                      Y = Y,
    #                      mi_vec = mi_vec,
    #                      i_index = i_index,
    #                      phi = phi,
    #                      beta = beta,
    #                      Z = Z,
    #                      B = B,
    #                      K = K)
  }
  #slower
  #hessian_il <- -partials_il %*% Vi_inv %*% t(partials_il)

  #faster
  hessian_il <- -partials_il %*% tcrossprod(Vi_inv, partials_il)

  d_hessian_il <- diag(hessian_il)
  return(d_hessian_il)
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
#' @param i_index
#' @param V_inv (optional) List of V inverses for each i. If supplied is faster
#'
#' @returns Matrix of dimension KP times KP
#' @export
#'
get_Hessian_l <- function(l,
                          Y,
                          mi_vec,
                          i_index,
                          beta,
                          Z,
                          B,
                          phi,
                          C,
                          V_inv,
                          partials_l) {
  n <- length(mi_vec)
  #hessian_l <- matrix(0, nrow = nrow(beta), ncol = nrow(beta))
  d_hessian_l <- numeric(nrow(beta))
  # Hessian will be diagonalized, so we can treat it as a vector
  for (i in 1:n) {
    dhessian_il <- get_dHessian_il(i = i,
                                   l = l,
                                   Y = Y,
                                   mi_vec = mi_vec,
                                   i_index = i_index,
                                   beta = beta,
                                   Z = Z,
                                   B = B,
                                   phi = phi,
                                   Vi_inv = V_inv[[i]],
                                   partials_il = partials_l[[i]])
    d_hessian_l <- d_hessian_l + dhessian_il
  }

  diag_hessian <- -pmax(-d_hessian_l, C)
  C_hessian <- diag(diag_hessian)

  return(C_hessian)
}
