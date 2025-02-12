# Various -----------------------------------------------------------------
#' Get \eqn{Y_{ij0}}$
#'
#' @param i subject
#' @param j time index
#' @param Y response
#' @param mi number of time points for subject i
#'
#' @returns
#'
#' @examples
get_Yij0 <- function(i, j, Y, mi) {
  Yij0 <- sum(Y[((i - 1) * mi + j), ])
  return(Yij0)
}



#' Get \eqn{B_ij}$
#'
#' @param i subject
#' @param j time index
#' @param B B-spline basis matrix
#' @param mi number of time points for subject i
#'
#' @returns
#'
#' @examples
get_Bij <- function(i, j, B, mi) {
  Bij <- B[((i - 1) * mi + j), ]
  return(Bij)
}





#' Get \eqn{Z_{ijl}}$
#'
#' @param i subject
#' @param j time index
#' @param l external variable index
#' @param Z matrix of external variables
#'
#' @returns
#'
#' @examples
get_Z_ijl <- function(i, j, l, Z) {
  Z_ijl <- Z[i + (j - 1), (l + 1)]
  return(Z_ijl)
}
