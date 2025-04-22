
#' cleverly_bestpsi
#'
#' Run cleverly for a range of hyperparameter psi values.
#'
#' @param psi_min minimum psi
#' @param psi_max maximum psi
#' @param npsi Number of psi values to test
#' @param parralel Run in parallel? T/F, if T, use all available cores
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param time vector of time values for each subject/time
#' @param lp clustering index (integer between 0 and L)
#' @param response_type Counts or continuous response
#' @param cor_str Type of correlation structure (IND, CON, AR1)
#' @param gammas Vector of dimension L + 1 for penalizing the D matrix
#' @param tau MCP hyper parameter.
#' @param theta ADMM hyper parameter.
#' @param C Constant for determining the hessian change.
#' @param d Order for the difference matrix
#' @param nknots Number of knots for the B-spline basis
#' @param order Order of the B-spline basis
#' @param epsilon_b Tolerance for alg 1 convergence
#' @param epsilon_r Tolerance for ADMM convergence
#' @param epsilon_d Tolerance for ADMM convergence
#' @param max_outer_iter Max number of iterations for the outer loop (Algorithm 1)
#' @param max_admm_iter Max number of iterations for the ADMM loop
#' @param max_2_iter Maximum number of iterations for algorithm 2 to run each loop
#' @param epsilon_2 Tolerance for convergence of algorithm 2
#' @param run_min
#'
#' @returns list of
#' @export
cleverly_bestpsi <- function(psi_min,
                             psi_max,
                             npsi,
                             parralel = FALSE,
                             Y,
                             Z,
                             time,
                             lp = 0,
                             response_type = "counts",
                             cor_str = "IND",
                             gammas,
                             tau = 8/100,
                             theta = 300,
                             C = 10,
                             d = 2,
                             run_min = 3,
                             nknots = 3,
                             order = 3,
                             epsilon_b = 1e-3,
                             epsilon_r = 1e-3,
                             epsilon_d = 1e-3,
                             max_outer_iter = 10,
                             max_admm_iter = 100,
                             max_2_iter = 100,
                             epsilon_2 = 1e-3){

  psis <- seq(psi_min, psi_max, length.out = npsi)


  if (parralel) {
    future::plan(future::multisession, workers = future::availableCores())
    res_list <- furrr::future_map(psis, ~cleverly(Y = Y,
                                                  Z = Z,
                                                  subject_ids = individual,
                                                  lp = 0,
                                                  time = time,
                                                  # Hyperparameters
                                                  gammas = gammas,
                                                  tau = tau,
                                                  theta = theta,
                                                  psi = ..1,
                                                  C = 100,
                                                  run_min = run_min,
                                                  # Iterations max
                                                  max_admm_iter = max_admm_iter,
                                                  max_outer_iter = max_outer_iter,
                                                  max_2_iter = max_2_iter,
                                                  # Convergence criteria
                                                  epsilon_r = .001,
                                                  epsilon_d = .05,
                                                  epsilon_b = .01,
                                                  epsilon_2 = .001,
                                                  cor_str = cor_str))
  } else {
    res_list <- list()
    for (p in 1:length(psis)) {
      psi <- psis[p]

      res_list[[p]] <- cleverly(Y = Y,
                                Z = Z,
                                subject_ids = individual,
                                time = time,
                                lp = lp,
                                response_type = response_type,
                                cor_str = cor_str,
                                gammas = gammas,
                                psi = psi,
                                tau = tau,
                                theta = theta,
                                C = C,
                                d = d,
                                run_min = run_min,
                                nknots = nknots,
                                order = order,
                                epsilon_b = epsilon_b,
                                epsilon_r = epsilon_r,
                                epsilon_d = epsilon_d,
                                max_outer_iter = max_outer_iter,
                                max_admm_iter = max_admm_iter,
                                max_2_iter = max_2_iter,
                                epsilon_2 = epsilon_2)
    }
  }

  best <- which.min(purrr::map_dbl(res_list, ~.x$BIC))
  res <- res_list[[best]]
  res$all_clusters_psi <- purrr::map(res_list, ~.x$clusters)


  print(paste0("all clusters: psi:", psis,", cluster:", purrr::map(res_list, ~.x$clusters)))
  print(paste0("chosen psi cluster", purrr::map_dbl(res_list, ~.x$clusters$no)[best]))
  print(paste0("chosen psi", psis[best]))


  sim_result <- list("chosen_cluster" = res$clusters,
                     "possible_cluster" = res$all_clusters_psi,
                     "chosen_psi" = psis[best],
                     "y_hat_init" = res$y_hat_init,
                     "y_hat_final" = res$y_hat)

  cluster <- res$clusters$membership
  true_cluster <- rep(1:3, each = 4)
  sim_result$cluster_result <- data.frame("rand" = fossil::rand.index(cluster, true_cluster),
                                          "adj.rand" = mclust::adjustedRandIndex(cluster, true_cluster),
                                          "jacc" = length(intersect(cluster, true_cluster)) /
                                            length(union(cluster, true_cluster)),
                                          "miss" = mclust::classError(classification = cluster,
                                                                               class = true_cluster)$errorRate,
                                          "nclust" = res$clusters$no
                                          )

  return(sim_result)


}



#' cleverly_bestparam
#'
#' @param psi_min minimum psi value
#' @param psi_max maximum psi value
#' @param npsi Number of psi values to test
#' @param parralel Run in parallel? T/F, if T, use all available cores
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param time vector of time values for each subject/time
#' @param lp clustering index (integer between 0 and L)
#' @param response_type Counts or continuous response
#' @param cor_str Type of correlation structure (IND, CON, AR1)
#' @param gammas Vector of dimension L + 1 for penalizing the D matrix
#' @param tau MCP hyper parameter.
#' @param theta ADMM hyper parameter.
#' @param C Constant for determining the hessian change.
#' @param d Order for the difference matrix
#' @param nknots Number of knots for the B-spline basis
#' @param order Order of the B-spline basis
#' @param epsilon_b Tolerance for alg 1 convergence
#' @param epsilon_r Tolerance for ADMM convergence
#' @param epsilon_d Tolerance for ADMM convergence
#' @param max_outer_iter Max number of iterations for the outer loop (Algorithm 1)
#' @param max_admm_iter Max number of iterations for the ADMM loop
#' @param max_2_iter Maximum number of iterations for algorithm 2 to run each loop
#' @param epsilon_2 Tolerance for convergence of algorithm 2
#'
#' @returns list
#' @export
cleverly_bestparam <- function(param_grid,
                               parralel = FALSE,
                               Y,
                               Z,
                               time,
                               lp = 0,
                               response_type = "counts",
                               cor_str = "IND",
                               gamma_length,
                               theta = 300,
                               C = 10,
                               d = 2,
                               nknots = 3,
                               order = 3,
                               epsilon_b = 1e-3,
                               epsilon_r = 1e-3,
                               epsilon_d = 1e-3,
                               max_outer_iter = 10,
                               max_admm_iter = 100,
                               max_2_iter = 100,
                               epsilon_2 = 1e-3){



  if (parralel) {
    future::plan(future::multisession, workers = future::availableCores())
    res_list <- furrr::future_pmap(param_grid, ~cleverly(Y = Y,
                                                  Z = Z,
                                                  subject_ids = individual,
                                                  lp = 0,
                                                  time = time,
                                                  # Hyperparameters
                                                  gammas = rep(..3, gamma_length),
                                                  tau = ..2,
                                                  theta = theta,
                                                  psi = ..1,
                                                  C = 100,
                                                  # Iterations max
                                                  max_admm_iter = max_admm_iter,
                                                  max_outer_iter = max_outer_iter,
                                                  max_2_iter = max_2_iter,
                                                  # Convergence criteria
                                                  epsilon_r = .001,
                                                  epsilon_d = .05,
                                                  epsilon_b = .01,
                                                  epsilon_2 = .001,
                                                  cor_str = cor_str))
  } else {
    res_list <- list()
    for (p in 1:length(param_grid)) {
      psi <- param_grid[p]

      res_list[[p]] <- cleverly(Y = Y,
                                Z = Z,
                                subject_ids = individual,
                                time = time,
                                lp = lp,
                                response_type = response_type,
                                cor_str = cor_str,
                                gammas = gammas,
                                psi = psi,
                                tau = tau,
                                theta = theta,
                                C = C,
                                d = d,
                                nknots = nknots,
                                order = order,
                                epsilon_b = epsilon_b,
                                epsilon_r = epsilon_r,
                                epsilon_d = epsilon_d,
                                max_outer_iter = max_outer_iter,
                                max_admm_iter = max_admm_iter,
                                max_2_iter = max_2_iter,
                                epsilon_2 = epsilon_2)
    }
  }


  best <- which.min(purrr::map_dbl(res_list, ~.x$BIC))
  res <- res_list[[best]]
  res$all_clusters_psi <- purrr::map(res_list, ~.x$clusters)


  sim_result <- list("chosen_cluster" = res$clusters,
                     "possible_cluster" = res$all_clusters_psi,
                     "chosen_params" = param_grid[best,],
                     "y_hat_init" = res$y_hat_init,
                     "y_hat_final" = res$y_hat)

  cluster <- res$clusters$membership
  true_cluster <- c(1, 1, 1, 1,
                    2, 2, 2, 2,
                    3, 3, 3, 3)

  cluster_table <- table(cluster, rep(1:3, each = 4))
  miss_rate <- sum(diag(cluster_table))/sum(cluster_table)

  sim_result$cluster_result <- data.frame("rand" = fossil::rand.index(cluster, true_cluster),
                                          "adj.rand" = mclust::adjustedRandIndex(cluster, true_cluster),
                                          "jacc" = length(intersect(cluster, true_cluster)) /
                                            length(union(cluster, true_cluster)),
                                          "miss" = miss_rate
  )

  return(sim_result)


}


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
