# Algorithm 3 -------------------------------------------------------------

#' ADMM algorithm for updating the clustering beta
#'
#' Iterates between updating v, beta, and lambda
#'
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param lp clustering index (integer between 0 and L)
#' @param B B spline basis matrix of dimension (N x P)
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param A Matrix keeping track of pairs of responses to calculate pairwise differences
#' @param tau MCP hyper parameter.
#' @param theta ADMM hyper parameter.
#' @param psi Hyperparameter for clustering penalty (larger drives pairwise differences to zero)
#' @param max_admm_iter Max number of iterations for the ADMM loop
#' @param P Number of B-spline coefficients (order + nknots)
#' @param C Constant for determining the hessian change.
#' @param D Matrix of dth order weighted difference operator
#' @param Kappa keeps track of each possible response pair
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param gammas Vector of dimension L + 1 for penalizing the D matrix
#' @param phi Current value of overdispersion parameter
#' @param epsilon_r Tolerance for ADMM convergence
#' @param epsilon_d Tolerance for ADMM convergence
#'
#' @returns List with the updated betas, updated vs, and lists tracking the betas, vs, and lambdas for each iteration
#' @export
algorithm3 <- function(Y,
                       Z,
                       lp,
                       B,
                       beta,
                       lambda,
                       A,
                       P,
                       C,
                       D,
                       Kappa,
                       mi_vec,
                       gammas,
                       tau,
                       theta,
                       psi,
                       phi,
                       max_admm_iter,
                       epsilon_r,
                       epsilon_d) {
  #print("Starting Algorithm 3")
  # Iterations for admm loop
  t <- 1
  K <- ncol(Y)
  # Initialize all return lists
  beta_admm_track <- list()
  v_admm_track <- list()
  lambda_admm_track <- list()
  cluster_list <- list()
  r_list <- list()
  d_list <- list()
  diff_admm <- c()
  phi_track <- c()

    pb_admm <- utils::txtProgressBar(min = 0,
                              max = max_admm_iter,
                              style = 3)

  # Initialize lambda to all 0
  #lambda <- numeric(nrow(Kappa)*P)
  # Initialize v to be zero just for checking the difference
  # The first v will be based on the current beta.
  v <- numeric(nrow(Kappa)*P)
  pb <- utils::txtProgressBar(min = 0,
                              max = max_admm_iter,
                              style = 3)
  while (t <= max_admm_iter) {
    # if (t == 1 | t %% 5 == 0) {
    #   print(paste("ADMM Iteration: ", t))
    # }

    v_new <- update_v(beta = beta,
                      lp = lp,
                      lambda = lambda,
                      Kappa = Kappa,
                      P = P,
                      gammas = gammas,
                      tau = tau,
                      theta = theta,
                      psi = psi)

    beta_new <- update_beta_admm(Y = Y,
                                 beta = beta,
                                 lp = lp,
                                 v = v_new,
                                 lambda = lambda,
                                 theta = theta,
                                 gammas = gammas,
                                 D = D,
                                 A = A,
                                 Z = Z,
                                 B = B,
                                 phi = phi,
                                 C = C,
                                 mi_vec = mi_vec)
    lambda_new <- update_lambda(beta_lp = beta_new[,lp + 1],
                                v = v_new,
                                lambda = lambda,
                                A = A,
                                theta = theta)

    # Update dispersion parameter
    phi <- get_phi(Y = Y,
                   beta = beta_new,
                   Z = Z,
                   B = B,
                   K = ncol(Y),
                   mi_vec = mi_vec)

    # Keep track of loop
    phi_track <- c(phi_track, phi)
    beta_admm_track[[t]] <- beta_new
    v_admm_track[[t]] <- v_new
    lambda_admm_track[[t]] <- lambda_new
    cluster_list[[t]] <- get_clusters(v_new, K, P)

    # Check for convergence
    diff <- sum(abs(beta_new - beta))
    diff_admm <- c(diff_admm, diff)
    r_norm <- get_r_norm(
      beta = beta_new,
      v = v_new,
      Kappa = Kappa,
      P = P,
      lp = lp)


    d_norm <- get_d_norm(v_new = v_new,
                         v = v,
                         theta = theta,
                         Kappa = Kappa,
                         K = ncol(beta),
                         P = P)


    r_list[[t]] <- r_norm
    d_list[[t]] <- d_norm

    if (r_norm <= epsilon_r & d_norm <= epsilon_d) {
      # exit the loop
      break
    }

    if (t == 20) {
      #browser()
    }
    # Prepare for next iteration
    v <- v_new
    beta <- beta_new
    lambda <- lambda_new
    t <- t + 1
    utils::setTxtProgressBar(pb_admm, t)
    # pb_admm$tick()


  }
  utils::setTxtProgressBar(pb_admm, max_admm_iter)
  #print(paste("Last ADMM Iteration: ", t - 1))
  return(list(beta = beta,
              lambda = lambda,
              v = v,
              beta_admm_track = beta_admm_track,
              lambda_admm_track = lambda_admm_track,
              v_admm_track = v_admm_track,
              cluster_list = cluster_list,
              phi_track = phi_track,
              diff_admm = diff_admm,
              r_list = r_list,
              d_list = d_list))
}



# Get v, beta, lambda -----------------------------------------------------

#' Update v
#'
#' Update v step of ADMM algorithm. V is the pairwise differences vector of dimension |Kappa|P
#'
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param lp clustering index (integer between 0 and L)
#' @param lambda ADMM vector of dimension (P*|Kappa|) x1
#' @param Kappa keeps track of each possible response pair
#' @param gammas Vector of dimension L + 1 for penalizing the D matrix
#' @param tau MCP hyper parameter.
#' @param theta ADMM hyper parameter.
#' @param psi Hyperparameter for clustering penalty (larger drives pairwise differences to zero)
#' @param P Number of B-spline coefficients (order + nknots)
#'
#' @returns Vector of updated vs for ADMM step
#' @export
update_v <- function(beta,
                     lp,
                     lambda,
                     Kappa,
                     P,
                     gammas,
                     tau,
                     theta,
                     psi){

  mcp <- tau * theta / (tau * theta - 1)
  sigma <- psi/theta
  v <- numeric(nrow(Kappa) * P)

  for (kappa in 1:nrow(Kappa)) {
    # Get the pairwise betas
    k1 <- Kappa[kappa, 1]
    k2 <- Kappa[kappa, 2]
    beta_k1 <- beta[(P * (k1 - 1) + 1):(P * k1), lp + 1]
    beta_k2 <- beta[(P * (k2 - 1) + 1):(P * k2), lp + 1]
    # Get corresponding lambdas
    lambda_kappa <- lambda[(P * (kappa - 1) + 1):(P * kappa)]
    # BIG CHANGE HERE! ITS PLUS NOT MINUS TO MATCH SIGN
    u <- beta_k1 - beta_k2 + lambda_kappa/theta

    norm_u <- sum(u^2) # or norm(u, type = "2")
    # Also check if the norm is 0 to avoid dividing by zero
    if (norm_u >= tau * psi | norm_u == 0) {
      v[((kappa - 1) * P + 1):(kappa * P)] <- u
    } else {
      v[((kappa - 1) * P + 1):(kappa * P)] <- mcp * max(0, ( 1 - sigma/norm_u)) * u
    }
  }

  if (any(is.nan(v))) {
    stop("v is NaN")
  }

  return(v)
}


#' Update beta_lp in the ADMM algorithm
#'
#' @param Y Matrix of counts Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param lp clustering index (integer between 0 and L)
#' @param v differences for each pair of betas of dimension (P*|Kappa|) x1
#' @param lambda ADMM vector of dimension (P*|Kappa|) x1
#' @param theta ADMM hyper parameter.
#' @param gammas Vector of dimension L + 1 for penalizing the D matrix
#' @param D Matrix of dth order weighted difference operator
#' @param A Matrix keeping track of pairs of responses to calculate pairwise differences
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param phi Current value of overdispersion parameter
#' @param C Constant for determining the hessian change.
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns matrix of updated betas
#' @export
update_beta_admm <- function(Y,
                             beta,
                             lp,
                             v,
                             lambda,
                             theta,
                             gammas,
                             D,
                             A,
                             Z,
                             B,
                             phi,
                             C,
                             mi_vec){


  gamma <- gammas[lp + 1]
  H <- get_Hessian_l(l = lp,
                     Y = Y,
                     mi_vec = mi_vec,
                     beta = beta,
                     Z = Z,
                     B = B,
                     phi = phi,
                     C = C)
  Q <- get_gradient_l(Y = Y,
                      mi_vec = mi_vec,
                      l = lp,
                      phi = phi,
                      beta = beta,
                      Z = Z,
                      B = B)

  v_tilde <- v - lambda/theta
  beta_lp <- beta[,lp + 1]

  first_term <- -H + gamma * D + theta * t(A) %*% A
  second_term <- Q - H %*% beta_lp + theta * t(A) %*% v_tilde
  beta_lp_new <- MASS::ginv(first_term) %*% second_term

  # All other beta elements are fixed, only the lp column is updated.
  beta[,lp + 1] <- beta_lp_new


  if (any(is.nan(beta_lp_new))) {
    stop("Beta ADMM is NaN")
  }

  return(beta)
}

#' Update lambda step of the admm algorithm
#'
#' @param beta_lp beta_lp is the lp column of beta
#' @param v differences for each pair of betas of dimension (P*|Kappa|) x1
#' @param lambda value of lambda from the previous iteration
#' @param A Matrix keeping track of pairs of responses to calculate pairwise differences
#' @param theta ADMM hyper parameter.

#'
#' @returns vector of updated lambdas for the ADMM step
#' @export
update_lambda <- function(beta_lp, v, lambda, A, theta){
  lambda_new <- lambda + theta * (A %*% beta_lp - v)

  if (any(is.nan(lambda_new))) {
    stop("Lambda is NaN")
  }
  return(lambda_new)
}




# Convergence criteria ----------------------------------------------------


get_r_norm <- function(beta, v, Kappa, P, lp){
  r_kappa <- c()
  for (kappa in 1:nrow(Kappa)) {
    k1 <- Kappa[kappa,1]
    k2 <- Kappa[kappa,2]

    beta_k1 <- beta[(P * (k1 - 1) + 1):(P * k1), lp + 1]
    beta_k2 <- beta[(P * (k2 - 1) + 1):(P * k2), lp + 1]
    v_kappa <- v[((kappa - 1)*P + 1):(kappa*P)]

    diff <- beta_k1 - beta_k2 - v_kappa

    r_kappa <- c(r_kappa, diff)
  }
  norm <- sum(r_kappa^2)

  return(norm)
}

get_d_norm <- function(v_new, v, theta, Kappa, K, P){
  v_diff <- v_new - v
  diff_mat <- matrix(v_diff, nrow = P)
  s <- 0

  # for (k in K) {
  #   id_k1 <- which(Kappa[,1] == k)
  #   id_k2 <- which(Kappa[,2] == k)
  # }
  for (i in 1:K) {
    v1 <- rep(0, P)
    v2 <- rep(0, P)
    for (ix in 1:nrow(Kappa)) {
      l1 <- Kappa[ix, 1]
      l2 <- Kappa[ix, 2]
      if (l1 == i) {
        v1 <- v1 + diff_mat[, ix]
      } else {
        v1 <- v1
      }
      if (l2 == i) {
        v2 <- v2 + diff_mat[, ix]
      } else {
        v2 <- v2
      }
    }
    v <- (v1 - v2)^2
    s <- s + sum(v)
  }

  residual <- theta * sqrt(s)
  return(residual)

  #return(norm)
}
