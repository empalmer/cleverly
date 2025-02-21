#' Algorithm 2: Updating the non clustering betas
#'
#' @param Y
#' @param Z
#' @param mi_vec
#' @param lp
#' @param B
#' @param beta
#' @param D
#' @param gammas
#' @param phi
#' @param C
#'
#' @returns Vector of length PK x L
#' @export
#'
#' @examples
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

