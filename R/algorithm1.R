#' Algorithm 1
#'
#' @param Y Matrix of counts Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param Z matrix that starts with a column of 1s. Of dimension M x L + 1
#' @param time
#' @param mi_vec count vector of length n for the number of time points for each i
#' @param lp clustering index
#' @param gammas
#' @param psi
#' @param phi
#' @param tau
#' @param theta
#' @param d Order for the difference matrix
#' @param nknots Number of knots for the B-spline basis
#' @param order Order of the B-spline basis
#' @param tol Tolerance for convergence
#' @param max_outer_iter
#' @param max_admm_iter

#'
#' @returns
#' @export
#'
#' @examples
algorithm1 <- function(Y,
                       Z,
                       time,
                       mi_vec,
                       lp,
                       gammas,
                       psi,
                       phi,
                       tau,
                       theta,
                       d,
                       nknots,
                       order,
                       tol,
                       max_outer_iter,
                       max_admm_iter) {

  # Get algorithm constants.
  P <- nknots + order
  L <- ncol(Z) - 1
  K <- ncol(Y)

  # Get indices for the non-cluster responses
  lp_minus <- setdiff(0:L , lp)

  # Calculate B-spline basis based on time for each subject/time
  B <- get_B(time = time,
             order = order,
             nknots = nknots)
  # Calculate the difference operator
  D <- get_D(K = K,
             d = d,
             order = order,
             nknots = nknots)
  # Create the matrix keeping track of pairwise difference indices
  Kappa <- t(combn(K, 2))
  A <- get_A(Kappa = Kappa,
             K = K,
             P = P)

  # Initialize beta vector to 0s.
  beta <- initialize_beta(K = K, L = L, P = P)
  # loop initialization
  loop_list_beta <- list()
  loop_list_diff <- list()
  diff <- 100
  s <- 1
  while ((diff > tol) & (s < max_outer_iter)) {
    # Solution for l_p minus, with l_p fixed
    beta_lp_minus <- algorithm2(Y = Y,
                                Z = Z,
                                mi_vec = mi_vec,
                                lp = lp,
                                B = B,
                                beta = beta,
                                D = D,
                                gammas = gammas,
                                phi = phi,
                                C = C)
    # Solution for l p (ADMM) with beta l_p minus fixed
    beta_lp <- algorithm3(Y = Y,
                          Z = Z,
                          lp = lp,
                          B = B,
                          beta = beta_lp_minus,
                          A = A,
                          P = P,
                          Kappa = Kappa,
                          gammas = gammas,
                          tau = tau,
                          theta = theta,
                          psi = psi,
                          max_admm_iter = max_admm_iter)


    # Difference in betas between this loop and the last
    diff <- sum(abs(beta_lp - beta)) # matrix difference
    beta <- beta_lp

    loop_list_beta[s] <- beta
    loop_list_diff[s] <- diff

    # Increment loop
    s <- s + 1
  }

  return(list(beta = beta))
}







