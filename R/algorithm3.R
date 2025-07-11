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
#' @param Y0 Vector of total count for each sample
#' @param lambda ADMM vector of dimension (P*|Kappa|) x1
#' @param AtA Pre-calculated A transpose times A for speed
#' @param i_index starting index of the ith subject in the data
#' @param s Outer iteration number
#' @param L Number of external variables
#' @param K Number of responses
#' @param M Number of samples times timepoints for each sample
#' @param cor_str specified correlation structure
#' @param cor_blocks pre-calculated correlation structure
#' @param j1_j2_list used for AR1 and AR1 d correlations
#'
#' @returns List with the updated betas, updated vs, and lists tracking the betas, vs, and lambdas for each iteration
algorithm3 <- function(Y,
                       Y0,
                       Z,
                       lp,
                       B,
                       beta,
                       lambda,
                       A,
                       AtA,
                       P,
                       C,
                       D,
                       Kappa,
                       mi_vec,
                       i_index,
                       gammas,
                       tau,
                       theta,
                       psi,
                       phi,
                       max_admm_iter,
                       epsilon_r,
                       epsilon_d,
                       s,
                       L,
                       K,
                       M,
                       cor_str,
                       cor_blocks,
                       j1_j2_list) {

  # Initialize all return lists
  beta_admm_track <- list()
  v_admm_track <- list()
  lambda_admm_track <- list()
  cluster_list <- list()
  r_list <- list()
  u_list <- list()
  d_list <- list()
  diff_admm <- numeric(max_admm_iter)
  phi_track <- numeric(max_admm_iter)

  # cat(paste0("ADMM iteration: ", s, "\n"))
  #
  # # Initialize progress bar
  # pb_admm <- utils::txtProgressBar(min = 0,
  #                                  max = max_admm_iter,
  #                                  style = 3)
  # Initialize v to be zero just for checking the difference
  # The first v will be based on the current beta.
  v <- numeric(nrow(Kappa)*P)


  for (t in 1:max_admm_iter) {

    alpha <- get_alpha_list(beta = beta,
                            Z = Z,
                            B = B,
                            K = K,
                            i_index = i_index,
                            mi_vec = mi_vec,
                            L = L,
                            P = P)

    # Calculate pearson residuals (used to calculate phi and rho)
    # Creates a list of length i
    pearson_residuals <- get_pearson_residuals(Y = Y,
                                               Y0 = Y0,
                                               beta = beta,
                                               alpha = alpha,
                                               Z = Z,
                                               B = B,
                                               K = K,
                                               mi_vec = mi_vec,
                                               i_index = i_index,
                                               M = M)
    # Update dispersion parameter
    phi <- get_phi(pearson_residuals = pearson_residuals,
                   K = K,
                   M = M,
                   L = L,
                   P = P)

    # Update correlation parameter
    rho_cor <- get_rho(pearson_residuals = pearson_residuals,
                       phi = phi,
                       K = K,
                       mi_vec = mi_vec,
                       M = M,
                       cor_str = cor_str,
                       cor_blocks = cor_blocks,
                       j1_j2_list = j1_j2_list)

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
                                 Y0 = Y0,
                                 beta = beta,
                                 alpha = alpha,
                                 lp = lp,
                                 v = v_new$v,
                                 lambda = lambda,
                                 theta = theta,
                                 gammas = gammas,
                                 D = D,
                                 A = A,
                                 AtA = AtA,
                                 Z = Z,
                                 B = B,
                                 phi = phi,
                                 C = C,
                                 mi_vec = mi_vec,
                                 i_index = i_index,
                                 K = K,
                                 P = P,
                                 L = L,
                                 cor_str = cor_str,
                                 rho_cor = rho_cor)

    lambda_new <- update_lambda(beta_lp = beta_new[,lp + 1],
                                v = v_new$v,
                                lambda = lambda,
                                A = A,
                                theta = theta)

    # Keep track of loop
    phi_track[t] <- phi
    beta_admm_track[[t]] <- beta_new
    v_admm_track[[t]] <- v_new
    lambda_admm_track[[t]] <- lambda_new
    cluster_list[[t]] <- get_clusters(v_new$v, K, P)
    u_list[[t]] <- v_new$u_list

    # Check for convergence
    diff_admm[t] <- sum(abs(beta_new - beta))

    r_norm <- get_r_norm(beta = beta_new,
                         v = v_new$v,
                         Kappa = Kappa,
                         P = P,
                         lp = lp)

    d_norm <- get_d_norm(v_new = v_new$v,
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

    # Prepare for next iteration
    v <- v_new$v
    beta <- beta_new
    lambda <- lambda_new
    # utils::setTxtProgressBar(pb_admm, t)

  }
  # utils::setTxtProgressBar(pb_admm, max_admm_iter)
  # close(pb_admm)


  #print(paste("Last ADMM Iteration: ", t - 1))
  return(list(beta = beta,
              lambda = lambda,
              v = v,
              cluster_list = cluster_list
              ))
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
  u_list <- list()

  for (kappa in 1:nrow(Kappa)) {

    # Get the pairwise betas
    k1 <- Kappa[kappa, 1]
    k2 <- Kappa[kappa, 2]
    beta_k1 <- beta[(P * (k1 - 1) + 1):(P * k1), lp + 1]
    beta_k2 <- beta[(P * (k2 - 1) + 1):(P * k2), lp + 1]
    # Get corresponding lambdas
    lambda_kappa <- lambda[(P * (kappa - 1) + 1):(P * kappa)]
    # BIG CHANGE HERE (compared to COMPARING paper)! ITS PLUS NOT MINUS TO MATCH SIGN
    u <- beta_k1 - beta_k2 + lambda_kappa/theta

    norm_u <- sum(u^2) # or norm(u, type = "2")

    u_list[[kappa]] <- u

    # Also check if the norm is 0 to avoid dividing by zero
    if (norm_u >= tau * psi | norm_u == 0) {
      v[((kappa - 1) * P + 1):(kappa * P)] <- u
    } else {
      v[((kappa - 1) * P + 1):(kappa * P)] <- mcp * max(0, ( 1 - sigma/norm_u)) * u
    }
  }


  if (anyNA(v)) {
    stop("v is NaN")
  }

  return(list(v = v,
              u_list = u_list))
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
#' @param Y0 Vector of total count for each sample
#' @param alpha list of alpha that can be subsetted by i and j
#' @param AtA Pre-calculated A transpose times A for speed
#' @param i_index starting index of the ith subject in the data
#' @param K Number of responses
#' @param L Number of external variables
#' @param P Number of B-spline coefficients (order + nknots)
#' @param cor_str correlation structure
#' @param rho_cor value of rho
#'
#' @returns matrix of updated betas
update_beta_admm <- function(Y,
                             Y0,
                             beta,
                             alpha,
                             lp,
                             v,
                             lambda,
                             theta,
                             gammas,
                             D,
                             A,
                             AtA,
                             Z,
                             B,
                             phi,
                             C,
                             mi_vec,
                             i_index,
                             K,
                             L,
                             P,
                             cor_str,
                             rho_cor){


  gamma <- gammas[lp + 1]
  K <- ncol(Y)

  # Things to be calculated once per loop/beta update
  # Which are: alpha and V inverse and partials_l
  # alpha <- get_alpha_list(beta = beta,
  #                         Z = Z,
  #                         B = B,
  #                         K = K,
  #                         i_index = i_index,
  #                         mi_vec = mi_vec)
  # Get V inverse for all is
  # compute it just once first so we don't have to calculate it for both H and Q.

  V_inv <- get_V_inv(Y = Y,
                     Y0 = Y0,
                     mi_vec = mi_vec,
                     i_index = i_index,
                     phi = phi,
                     beta = beta,
                     Z = Z,
                     B = B,
                     K = ncol(Y),
                     alpha = alpha,
                     cor_str = cor_str,
                     rho_cor = rho_cor)

  partials_l <- get_partials_l_list(Y0 = Y0,
                                    l = lp,
                                    mi_vec = mi_vec,
                                    i_index = i_index,
                                    beta = beta,
                                    alpha = alpha,
                                    Z = Z,
                                    B = B,
                                    K = K,
                                    P = P)

  # Things we need to calculate new beta:
  H <- get_Hessian_l(l = lp,
                     Y0 = Y0,
                     mi_vec = mi_vec,
                     i_index = i_index,
                     beta = beta,
                     Z = Z,
                     B = B,
                     phi = phi,
                     C = C,
                     V_inv = V_inv,
                     partials_l = partials_l,
                     alpha = alpha,
                     K = K,
                     P = P)

  Q <- get_gradient_l(Y = Y,
                      Y0 = Y0,
                      mi_vec = mi_vec,
                      i_index = i_index,
                      l = lp,
                      phi = phi,
                      beta = beta,
                      alpha = alpha,
                      Z = Z,
                      B = B,
                      V_inv = V_inv,
                      partials_l = partials_l,
                      K = K,
                      P = P)
  v_tilde <- v - lambda/theta
  beta_lp <- beta[,lp + 1, drop = F]


  #first_term <- -H + gamma * D + theta * AtA
  #second_term <- Q - H %*% beta_lp + theta * crossprod(A, v_tilde)
  #beta_lp_new <- MASS::ginv(first_term) %*% second_term
  #beta_lp_new <- fast_mat_mult2(MASS::ginv(first_term), second_term)

  # theta_test_term <- -H + gamma * D
  # condition_num <- kappa(theta_test_term)


  # if (condition_num < 10) {
  #   theta
  # } else {
  #   theta <- eigen(A)
  # }

  beta_lp_new <- calculate_beta_lp_new(H = H,
                                       gamma = gamma,
                                       D = D,
                                       theta = theta,
                                       AtA = AtA,
                                       Q = Q,
                                       beta_lp = beta_lp,
                                       A = A,
                                       v_tilde = matrix(v_tilde))


  # All other beta elements are fixed, only the lp column is updated.
  beta[,lp + 1] <- beta_lp_new

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
update_lambda <- function(beta_lp, v, lambda, A, theta){
  #lambda_new <- lambda + theta * (A %*% beta_lp - v)
  lambda_new <- lambda + theta * (fast_mat_mult2(A, matrix(beta_lp)) - v)

  # if (any(is.nan(lambda_new))) {
  #   stop("Lambda is NaN")
  # }
  return(lambda_new)
}




# Convergence criteria ----------------------------------------------------


#' get_r_norm
#'
#' Get the L2 norm of convergence criteria r
#'
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param v differences for each pair of betas of dimension (P*|Kappa|) x1
#' @param Kappa keeps track of each possible response pair
#' @param P Number of B-spline coefficients (order + nknots)
#' @param lp clustering index (integer between 0 and L)
#'
#' @returns L2 norm of the difference between beta_k1 and beta_k2
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

#' get_d_norm
#'
#' @param v_new updated differences for each pair of betas of dimension (P*|Kappa|) x1
#' @param v differences for each pair of betas of dimension (P*|Kappa|) x1
#' @param theta ADMM hyper parameter.
#' @param Kappa keeps track of each possible response pair
#' @param K Number of responses
#' @param P Number of B-spline coefficients (order + nknots)
#'
#' @returns difference between v_new and v
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
