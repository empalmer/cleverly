

algorithm2 <- function(Y, Z, lp, B, beta, D, gammas){
  L <- ncol(Z) - 1
  lp_minus <- NULL

  # Loop only through non-clustering values
  l_loop <- setdiff(0:L , lp)


  for (l in l_loop) {
    gamma_l <- gammas[l]
    beta_l <- beta[,l]

    gradient_l <- NULL # function of beta_l_list not beta_l
    Hessian_l <- NULL # function of beta_l_list not beta_l

    first_term <- Hessian_l - gamma_l*D
    second_term <- -gradient_l + Hessian_l %*% beta_l
    beta_l_s <- solve(first_term) %*% second_term

    # Update:
    beta[,l + 1] <- beta_l_s
  }

  return(beta)
}

