#' Get hessian elements
#'
#' Get for for the ith subject deriving with respect to the lth external variable
#'
#' @param i
#' @param l
#' @param Y
#' @param mis
#' @param beta
#' @param Z
#' @param B
#' @param phi
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
#' @param l
#' @param Y
#' @param mis
#' @param beta
#' @param Z
#' @param B
#' @param phi
#' @param C
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
