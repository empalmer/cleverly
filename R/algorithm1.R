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
                       C,
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
  Kappa <- t(combn(K, 2))
  A <- get_A(Kappa = Kappa,
             K = K,
             P = P)

  # Initialize beta vector to 0s.
  beta <- initialize_beta(K = K, L = L, P = P)
  # loop initialization
  loop_list_beta <- list()
  loop_list_diff <- list()
  admm_beta_list <- list()
  admm_diffs <- list()
  phis_list <- list()
  diff <- 100
  s <- 1

  while ((diff > tol) & (s <= max_outer_iter)) {

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


  y_hat <- estimate_y(beta, B, Z, K)
  clusters <- get_clusters(beta)

  return(list(beta = beta,
              y_hat = y_hat,
              v = v,
              admm_diffs = admm_diffs,
              admm_beta_list = admm_beta_list,
              loop_list_beta = loop_list_beta,
              loop_list_diff = loop_list_diff,
              phis_list = phis_list))
}





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


get_clusters <- function(beta){
  L <- ncol(beta)/ncol(beta)
  clusters <- list()
  for (l in 1:L) {
    beta_l <- beta[,l]
    clusters[[l]] <- which(beta_l != 0)
  }
  return(clusters)
}


create_adjacency <- function(v, K) {
  differences <- apply(v, 2, FUN = function(x) {
    norm(as.matrix(x), "f")
  })
  connected_ix <- which(differences == 0)
  index <- t(combn(K, 2))
  i <- index[connected_ix, 1]
  j <- index[connected_ix, 2]
  A <- matrix(0, nrow = K, ncol = K)
  A[(j - 1) * K + i] <- 1
  return(A)
}
