#' Algorithm 1
#'
#' @param Y Matrix of counts Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param is vector of numeric ids for each subject
#' @param mis count vector of length n for the number of timepoints for each i
#' @param lp clustering index
#' @param Z matrix that starts with a column of 1s. Of dimension M x L + 1
#' @param gamma Vector of length (L + 1) of penalization hyper parameters
#' @param d Order for the difference matrix
#' @param nknots Number of knots for the B-spline basis
#' @param order Order of the B-spline basis
#' @param tol Tolerance for convergence
#' @param smax Maximum number of iterations
#'
#' @returns
#' @export
#'
#' @examples
algorithm1 <- function(Y,
                       Z,
                       is,
                       time,
                       mis,
                       lp,
                       gamma,
                       phi,
                       d,
                       nknots,
                       order,
                       tol,
                       smax) {

  P <- nknots + order
  L <- ncol(Z) - 1
  K <- ncol(Y)

  lp_minus <- setdiff(0:L , lp)

  beta <- initialize_beta(K = K, L = L, P = P)

  # Calculate B-spline basis based on time for each subject/time
  B <- get_B(time, order, nknots)


  diff <- 100
  s <- 0
  while ((diff > tol) & (s < smax)) {
    beta_lp_minus <- algorithm2(Y = Y,
                                Z = Z,
                                is = is,
                                mis = mis,
                                lp = lp,
                                B = B,
                                beta = beta,
                                D = D,
                                gamma = gamma,
                                phi)
    #beta_lp <- algorithm3(Y, Z, lp, B, beta_lp_minus, gamma) # DUMMY STAND IN


    #diff <- sum(abs(beta_lp - beta))
    #beta <- beta_lp
    beta <- beta_lp_minus
    #
    s <- s + 1
  }

  return(list(B = B,
              beta = beta,
              P = P,
              L = L,
              K = K))
}







