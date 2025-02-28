#' Algorithm 1
#'
#' @param Y Matrix of counts Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param time vector of time values for each subject/time
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param lp clustering index (integer between 0 and L)
#' @param gammas Vector of dimension L + 1 for penalizing the D matrix
#' @param psi Hyperparameter for clustering penalty (larger drives pairwise differences to zero)
#' @param tau MCP hyper parameter.
#' @param theta ADMM hyper parameter.
#' @param d Order for the difference matrix
#' @param nknots Number of knots for the B-spline basis
#' @param order Order of the B-spline basis
#' @param epsilon_b Tolerance for alg 1 convergence
#' @param epsilon_r Tolerance for ADMM convergence
#' @param epsilon_d Tolerance for ADMM convergence
#' @param max_outer_iter Max number of iterations for the outer loop (Algorithm 1)
#' @param max_admm_iter Max number of iterations for the ADMM loop
#' @param C Constant for determining the hessian change.
#'
#' @returns List of beta, clusters, y_hat, v, admm_diffs, admm_beta_list, loop_list_beta, loop_list_diff, phis_list
#' @export
algorithm1 <- function(Y,
                       Z,
                       time,
                       mi_vec,
                       lp,
                       gammas,
                       psi,
                       tau,
                       theta,
                       C,
                       d,
                       nknots,
                       order,
                       epsilon_b,
                       epsilon_r,
                       epsilon_d,
                       max_outer_iter,
                       max_admm_iter) {

  # Get algorithm constants.
  P <- nknots + order
  L <- ncol(Z) - 1
  K <- ncol(Y)
  M <- nrow(Y)

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
  Kappa <- t(utils::combn(K, 2))
  A <- get_A(Kappa = Kappa,
             K = K,
             P = P)

  # Initialize beta vector to 0s.
  beta <- initialize_beta(K = K, L = L, P = P)
  # Initialize phi to be 1
  #Dirichlet Multinomial over dispersion parameter.
  phi <- 1


  # loop initialization
  loop_list_beta <- list()
  loop_list_diff <- list()
  admm_beta_list <- list()
  admm_diffs <- list()
  phis_list <- list()
  diff <- 100
  s <- 1

  while ((diff > epsilon_b) & (s <= max_outer_iter)) {

    # Go straight to ADMM code if there is no external variables
    if (L > 0) {
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
    } else {
      beta_lp_minus <- beta
    }

    # Solution for l p (ADMM) with beta l_p minus fixed
    alg3 <- algorithm3(Y = Y,
                       Z = Z,
                       lp = lp,
                       B = B,
                       beta = beta_lp_minus,
                       A = A,
                       P = P,
                       C = C,
                       D = D,
                       Kappa = Kappa,
                       mi_vec = mi_vec,
                       gammas = gammas,
                       tau = tau,
                       theta = theta,
                       psi = psi,
                       phi = phi,
                       epsilon_r = epsilon_r,
                       epsilon_d = epsilon_d,
                       max_admm_iter = max_admm_iter)
    # phi <- estimate_phi()
    beta_lp <- alg3$beta
    admm_beta_list[[s]] <- alg3$beta_admm_track
    v <- alg3$v

    # Difference in betas between this loop and the last
    diff <- sum(abs(beta_lp - beta)) # matrix difference
    beta <- beta_lp

    loop_list_beta[[s]] <- beta
    loop_list_diff[[s]] <- diff
    admm_diffs[[s]] <- alg3$diff_admm
    phis_list[[s]] <- alg3$phi_track


    # Increment loop
    s <- s + 1
  }

  # After loop:
  # Calculte estimated ys and clusters
  y_hat <- estimate_y(beta = beta,
                      B = B,
                      Z = Z,
                      K = K)
  clusters <- get_clusters(v = v,
                           K = K,
                           P = P)

  return(list(beta = beta,
              clusters = clusters,
              y_hat = y_hat,
              v = v,
              admm_diffs = admm_diffs,
              admm_beta_list = admm_beta_list,
              loop_list_beta = loop_list_beta,
              loop_list_diff = loop_list_diff,
              phis_list = phis_list))
}





# Helpers for Algorithm 1 -------------------------------------------------
#' Estimate y-hat given beta
#'
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param B B spline basis matrix of dimension (N x P)
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param K Number of responses
#'
#' @returns Vector of RA values for y
#' @export
estimate_y <- function(beta, B, Z, K){
  P <- ncol(B)
  L <- ncol(Z) - 1
  M <- nrow(Z)
  yhat <- matrix(nrow = M, ncol = K)
  for (k in 1:K) {
    y_k <- numeric(M)
    for (l in 0:L) {

      beta_k <- beta[((k - 1)*P + 1):(k*P), l + 1]
      Z_l <- diag(Z[,l + 1])
      y_k <- y_k + Z_l %*% B %*% beta_k
    }
    yhat[,k] <- y_k

  }
  # Turn into RA version:
  yhat_ra <- exp(yhat)/rowSums(exp(yhat))
  return(yhat_ra)
}




# Return clusters ---------------------------------------------------------
#' Get cluster membership
#'
#' @param v final v (pairwise differences) from algorithm 1
#' @param K Number of responses
#' @param P Number of B-spline coefficients (order + nknots)
#'
#' @returns Cluster list information, indeces and number of clusters
#' @export
get_clusters <- function(v, K, P) {

  # Convert v into a matrix of P x length
  v_mat <- matrix(v, nrow = P)
  differences <- apply(v_mat, 2, FUN = function(x) {
    norm(as.matrix(x), "f")
  })
  connected_ix <- which(differences == 0)
  index <- t(utils::combn(K, 2))
  i <- index[connected_ix, 1]
  j <- index[connected_ix, 2]
  A <- matrix(0, nrow = K, ncol = K)
  A[(j - 1) * K + i] <- 1


  # Make graph from adjacency matrix
  graph <- igraph::graph_from_adjacency_matrix(
    A,
    mode = 'upper')
  #clustering membership
  clusters <- igraph::components(graph)
  return(clusters)
}
