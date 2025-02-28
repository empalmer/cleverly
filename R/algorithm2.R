#' Algorithm 2: Updating the non clustering betas
#'
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param lp clustering index (integer between 0 and L)
#' @param B B spline basis matrix of dimension (N x P)
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param D Matrix of dth order weighted difference operator
#' @param gammas Vector of dimension L + 1 for penalizing the D matrix
#' @param phi Current value of overdispersion parameter
#' @param C Constant for determining the hessian change.
#'
#' @returns Vector of length PK x L
#' @export
algorithm2 <- function(Y, Z, mi_vec, lp, B, beta, D, gammas, phi, C){
  L <- ncol(Z) - 1
  lp_minus <- NULL

  # Loop only through non-clustering values
  l_loop <- setdiff(0:L , lp)


  for (l in l_loop) {
    gamma_l <- gammas[l + 1]
    beta_l <- beta[, l + 1]


    # function of beta_l_list not beta_l
    gradient_l <- get_gradient_l(Y = Y,
                                 mi_vec = mi_vec,
                                 l = l,
                                 phi = phi,
                                 beta = beta,
                                 Z = Z,
                                 B = B) # function of beta_l_list not beta_l
    Hessian_l <- get_Hessian_l(l = l,
                               Y = Y,
                               mi_vec = mi_vec,
                               beta = beta,
                               Z = Z,
                               B = B,
                               phi = phi,
                               C = C)
    # function of beta_l_list not beta_l

    first_term <- -Hessian_l + gamma_l*D
    second_term <- gradient_l - Hessian_l %*% beta_l

    beta_l_s <- MASS::ginv(first_term) %*% second_term

    # Update:
    beta[,l + 1] <- beta_l_s
  }

  return(beta)
}

