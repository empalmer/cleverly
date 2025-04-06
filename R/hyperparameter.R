


# Get hyperparameters -----------------------------------------------------



get_hyperparameters <- function(Y, Z, mi_vec, lp, B, D, K, P){
  L <- ncol(Z) - 1
  beta_init <- matrix(0, nrow = K * P, ncol = L + 1)
  gammas <- get_gammas(Y, Z, mi_vec, lp, B, beta = beta_init, D, K, P)
  return(list(
    gammas = gammas
  ))
}


get_gammas <- function(Y, Z, mi_vec, lp, B, beta, D, K, P){

  gamma_grid <- seq(1, 200, length.out = 100)

  bic <- numeric(length(gamma_grid))

  for (g in 1:length(gamma_grid)) {
    gammas <- rep(gamma_grid[g], ncol(Z))
    fit <- algorithm2(Y = Y,
                      Z = Z,
                      mi_vec = mi_vec,
                      lp = -1,
                      B = B,
                      beta = beta,
                      D = D,
                      gammas = gammas,
                      phi = 1,
                      C = 1)
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



    y_hat <- estimate_y(beta = beta,
                        B = B,
                        Z = Z,
                        K = K,
                        Y = Y,
                        time = time)

    bic[g] <- fit
  }

  return(gammas_selected)


}



# Calculate cluster BIC ---------------------------------------------------

BIC_cluster <- function(y_ra_df,
                        K, n_clusters,
                        mi_vec,
                        nknots, order){

  N <- sum(mi_vec)
  yhat_y <- log(y_ra_df$yhat + .1) - log(y_ra_df$y + .1)
  first_term <- sum(yhat_y^2)/(N*K)
  second_term <- log(N * K) * n_clusters * (order + nknots)/(N*K)

  BIC <- log(first_term) + second_term
  return(BIC)
}
