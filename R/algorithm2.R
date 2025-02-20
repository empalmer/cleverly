

algorithm2 <- function(Y, Z, is, mis, lp, B, beta, D, gammas, phi){
  L <- ncol(Z) - 1
  lp_minus <- NULL

  # Loop only through non-clustering values
  l_loop <- setdiff(0:L , lp)


  for (l in l_loop) {
    gamma_l <- gammas[l]
    beta_l <- beta[,l]

    # function of beta_l_list not beta_l
    gradient_l <- get_gradient_l(Y = Y,
                                 is = is,
                                 mis = mis,
                                 l = l,
                                 phi = phi,
                                 beta = beta,
                                 Z = Z,
                                 B = B) # function of beta_l_list not beta_l
    Hessian_l <- NULL # function of beta_l_list not beta_l

    #first_term <- Hessian_l - gamma_l*D
    #second_term <- -gradient_l + Hessian_l %*% beta_l
    #beta_l_s <- solve(first_term) %*% second_term

    # Update:
    #beta[,l + 1] <- beta_l_s
  }

  return(beta)
}

