#' Algorithm 1: Outer Loop for Parameter Estimation and Clustering
#'
#' Executes the main iterative algorithm for estimating parameters and identifying clusters based on a penalized objective function with ADMM and MCP.
#'
#' @param Y An \eqn{M \times K} matrix of response variables (e.g., counts). Each column corresponds to a different response, and each row to a subject-time observation.
#' @param Z An \eqn{M \times (L + 1)} matrix of covariates. The first column must be a vector of 1s (intercept). If there are no external variables, \code{Z} should be a single column of 1s.
#' @param time A numeric vector of length \eqn{M} representing the time points associated with each row of \code{Y}.
#' @param mi_vec A vector of length \eqn{n} specifying the number of timepoints for each subject.
#' @param lp Clustering index. An integer between 0 and \eqn{L}, indicating which covariate to use for clustering. Use 0 for baseline clustering.
#' @param gammas A numeric vector of length \eqn{L + 1} for penalizing components of the difference matrix \code{D}.
#' @param psi Clustering penalty hyperparameter. Larger values encourage more pairwise similarities (i.e., more shrinkage toward common clusters).
#' @param tau MCP (minimax concave penalty) hyperparameter.
#' @param theta ADMM (alternating direction method of multipliers) hyperparameter.
#' @param d Integer specifying the order of the difference matrix.
#' @param nknots Number of knots used in the B-spline basis.
#' @param order Order of the B-spline basis functions.
#' @param epsilon_b Convergence tolerance for Algorithm 1.
#' @param epsilon_r ADMM convergence tolerance for the primal residual.
#' @param epsilon_d ADMM convergence tolerance for the dual residual.
#' @param max_outer_iter Maximum number of iterations for the outer loop (Algorithm 1).
#' @param max_admm_iter Maximum number of iterations for the ADMM loop (Algorithm 3).
#' @param C Constant for determining when the Hessian should be updated.
#' @param max_2_iter Maximum number of iterations for the inner loop (Algorithm 2) per outer iteration.
#' @param epsilon_2 Convergence tolerance for Algorithm 2.
#' @param cor_str Correlation structure for the working covariance model. One of \code{"IND"}, \code{"CON"}, or \code{"AR1"}.
#' @param run_min Minimum number of runs to ensure stability (used for internal control).
#'
#' @return A list containing:
#' \describe{
#'   \item{beta}{Estimated regression coefficients.}
#'   \item{clusters}{Final cluster assignments.}
#'   \item{y_hat}{Estimated fitted values.}
#'   \item{v}{Latent variable estimates from ADMM.}
#'   \item{admm_diffs}{History of ADMM convergence diagnostics.}
#'   \item{admm_beta_list}{List of beta estimates from each ADMM iteration.}
#'   \item{loop_list_beta}{List of beta estimates from each outer loop iteration.}
#'   \item{loop_list_diff}{List of beta differences between iterations (for diagnostics).}
#'   \item{phis_list}{List of phi estimates across iterations.}
#' }
#'
#' @export

algorithm1 <- function(Y,
                       Z,
                       time,
                       mi_vec,
                       lp,
                       cor_str,
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
                       epsilon_2,
                       run_min,
                       max_outer_iter,
                       max_admm_iter,
                       max_2_iter) {

  # Get algorithm constants.
  P <- nknots + order
  L <- ncol(Z) - 1
  K <- ncol(Y)
  M <- nrow(Y)



  # Get indices for the non-cluster responses
  lp_minus <- setdiff(0:L , lp)

  # Get indeces for where each i starts (since different mi for each i)
  # Pre-calculate to save time
  i_index <- c(0, cumsum(mi_vec))

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
  #pre-calculate for computation speed.
  AtA <- crossprod(A)
  Y0 <- rowSums(Y)



  # Calculate helpers for correlation matrices:
  if (cor_str == "AR1" | cor_str == "AR1-d") {
    j1_j2_list <- lapply(mi_vec, function(mi) {
      times <- 1:mi
      diag_values <- rep(times, each = K)
      outer(diag_values, diag_values, "-")
    })
  } else {
    j1_j2_list <- NULL
  }


  # Define a list with blocks for each i what the correlation matrix structure is
  # We set rho = 1 just to give the 0/1 matrix of zero/nonzero elements of the
  # correlation matrix depending on the chosen strucutre.
  if (cor_str == "IND") {
    cor_blocks <- NULL
  } else {
    cor_blocks <- purrr::map(mi_vec, ~get_corR(cor_str = cor_str,
                                               mi = .x,
                                               K = K,
                                               rho = 1))

  }

  # Initialize phi to be 1
  #Dirichlet Multinomial over dispersion parameter.
  phi <- 1

  # Initialize using Algorithm 2 for ALL responses
  print(paste0("Initializing beta values for psi = ", round(psi, 2)))
  zeros_beta <- matrix(0, nrow = K * P, ncol = L + 1)
  beta_init <- algorithm2(Y = Y,
                          Y0 = Y0,
                          Z = Z,
                          mi_vec = mi_vec,
                          i_index = i_index,
                          lp = -1,
                          B = B,
                          beta = zeros_beta,
                          D = D,
                          gammas = gammas,
                          phi = phi,
                          C = C,
                          max_2_iter = max_2_iter,
                          epsilon_2 = epsilon_2,
                          time = time,
                          s = "Initial fit",
                          L = L,
                          P = P,
                          K = K,
                          M = M,
                          cor_str = cor_str,
                          cor_blocks = cor_blocks,
                          j1_j2_list = j1_j2_list)
  beta <- beta_init$beta

  # Initialize lambda to all 0
  lambda <- numeric(nrow(Kappa)*P)
  # loop initialization
  error <- NULL
  alg1_beta <- list()
  alg1_diff <- list()
  alg_2_beta_diff <- list()
  admm_beta_list <- list()
  admm_diffs <- list()
  phis_list <- list()
  cluster_list <- list()
  r_list <- list()
  d_list <- list()
  rs <- list()
  ts <- list()
  diff <- Inf
  s <- 1
  rs[[1]] <- beta_init$r


  for (s in 1:max_outer_iter) {
    print(paste0("Outer loop iteration: ", s))
    beta_old <- beta
    # Go straight to ADMM code if there is no external variables
    if (L > 0) {
      # Solution for l_p minus, with l_p fixed
      alg2 <- tryCatch({
        algorithm2(Y = Y,
                   Y0 = Y0,
                   Z = Z,
                   mi_vec = mi_vec,
                   i_index = i_index,
                   lp = lp,
                   B = B,
                   beta = beta,
                   D = D,
                   gammas = gammas,
                   phi = phi,
                   C = C,
                   max_2_iter = max_2_iter,
                   epsilon_2 = epsilon_2,
                   time = time,
                   s = s,
                   L = L,
                   P = P,
                   K = K,
                   M = M,
                   cor_str = cor_str,
                   cor_blocks = cor_blocks,
                   j1_j2_list = j1_j2_list)
      }, error = function(e) {
        print(paste0("ERROR alg2: ", e$message))
        return(e$message)
      })
      rs[[s + 1]] <- alg2$r
      alg_2_beta_diff[[s]] <- alg2$beta_diff
      # If there is an error, exit out of the loop
      if (is.character(alg2)) {
        error <- alg2
        break
      }
      beta_lp_minus <- alg2$beta
    } else {
      beta_lp_minus <- beta
    }


    # Solution for l p (ADMM) with beta l_p minus fixed
    alg3 <- tryCatch({
      algorithm3(Y = Y,
                 Y0 = Y0,
                 Z = Z,
                 lp = lp,
                 B = B,
                 beta = beta_lp_minus,
                 lambda = lambda,
                 A = A,
                 AtA = AtA,
                 P = P,
                 C = C,
                 D = D,
                 Kappa = Kappa,
                 mi_vec = mi_vec,
                 i_index = i_index,
                 gammas = gammas,
                 tau = tau,
                 theta = theta,
                 psi = psi,
                 phi = phi,
                 epsilon_r = epsilon_r,
                 epsilon_d = epsilon_d,
                 max_admm_iter = max_admm_iter,
                 s = s,
                 K = K,
                 L = L,
                 M = M,
                 cor_str = cor_str,
                 cor_blocks = cor_blocks,
                 j1_j2_list = j1_j2_list)
    }, error = function(e) {
      print(paste0("ERROR alg3: ", e$message))
      return(e$message)
    })
    if (is.character(alg3)) {
      error <- alg3
      break
    }
    if (!is.list(alg3)) {
      stop("Error in algorithm 3")
    }

    beta <- alg3$beta
    lambda <- alg3$lambda
    v <- alg3$v

    alpha <- get_alpha_list(beta = beta,
                            Z = Z,
                            B = B,
                            K = K,
                            i_index = i_index,
                            mi_vec = mi_vec,
                            L = L,
                            P = P)
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

    # Difference in betas between this loop and the last
    diff <- sum(abs(beta - beta_old)) # matrix difference
    # alg1_diff[[s]] <- diff

    # Exit outer loop early if convergence is reached
    if (diff < epsilon_b) {
      break
    }

    print(alg3$cluster_list[[length(alg3$cluster_list)]]$membership)
    cluster_list[[s]] <- alg3$cluster_list[[length(alg3$cluster_list)]]
    # Exit if constant cluster results for the past 3 iterations
    if (s >= run_min) {
      current <-  cluster_list[[s]]$membership
      past1 <- cluster_list[[s - 1]]$membership
      past2 <- cluster_list[[s - 2]]$membership

      if (identical(current, past1) && identical(current, past2)) {
        # If the clusters are the same as the last two iterations, break
        print("Clusters not changing, exiting")
        break
      }
    }
  }


  # After loop:
  # Calculate v for the last time (since v is updated before beta in admm)
  v <- update_v(beta = beta,
                lp = lp,
                lambda = lambda,
                Kappa = Kappa,
                P = P,
                gammas = gammas,
                tau = tau,
                theta = theta,
                psi = psi)$v
  # Calculte estimated ys and clusters
  clusters <- get_clusters(v = v,
                           K = K,
                           P = P)

  # This is in relative abundance
  y_hat <- estimate_y(beta = beta,
                      B = B,
                      Z = Z,
                      K = K,
                      Y = Y,
                      time = time)
  y_hat_baseline <- estimate_y(beta = beta,
                               B = B,
                               Z = Z,
                               K = K,
                               Y = Y,
                               time = time,
                               baseline = T)



  y_hat_init <- estimate_y(beta = beta_init$beta,
                      B = B,
                      Z = Z,
                      K = K,
                      Y = Y,
                      time = time)
  BIC <- BIC_cluster(y_ra_df = y_hat,
                     K = K,
                     L = L,
                     n_clusters = clusters$no,
                     mi_vec = mi_vec,
                     nknots = nknots,
                     order = order)


# New fxns ----------------------------------------------------------------


  beta_group <- beta_cluster_group_update(y = Y,
                                          Z = Z,
                                          beta = beta,
                                          lp = lp,
                                          lp_minus = lp_minus,
                                          B = B,
                                          clusters = clusters,
                                          K = K, P = P, M = M)

  y_hat_lp_group <- estimate_y(beta = beta_group,
                               B = B,
                               Z = Z,
                               K = K,
                               Y = Y,
                               time = time,
                               baseline = T)

  BIC_ra_group <- BIC_cluster(y_ra_df = y_hat_lp_group,
                              K = K,
                              L = L,
                              n_clusters = clusters$no,
                              mi_vec = mi_vec,
                              nknots = nknots,
                              order = order)

  y_hat_counts_group <- estimate_y_counts(beta = beta_group,
                                          B = B,
                                          Z = Z,
                                          K = K,
                                          Y = Y,
                                          time = time)

  BIC_group <- BIC_cluster_group(y_hat_counts = y_hat_counts_group,
                                 y_counts = Y,
                                 beta_group = beta_group,
                                 clusters = clusters,
                                 M = M,
                                 K = K,
                                 nknots = nknots,
                                 order = order)


  alpha <- get_alpha_list(beta = beta,
                          Z = Z,
                          B = B,
                          K = K,
                          i_index = i_index,
                          mi_vec = mi_vec,
                          L = L,
                          P = P)
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

  rho <- get_rho(pearson_residuals = pearson_residuals,
                 phi = phi,
                 K = K,
                 mi_vec = mi_vec,
                 M = M,
                 cor_str = cor_str,
                 cor_blocks = cor_blocks,
                 j1_j2_list = j1_j2_list)


  return(list(clusters = clusters,
              y_hat = y_hat,
              y_hat_init = y_hat_init,
              y_hat_lp_group = y_hat_lp_group,
              y_hat_baseline = y_hat_baseline,
              y_hat_counts_group = y_hat_counts_group,
              beta = beta,
              v = v,
              rho = rho,
              phi = phi,
              BIC = BIC,
              BIC_group = BIC_group,
              BIC_ra_group = BIC_ra_group,
              s = s,
              error = error))
}




# Initializing algorithm 1 ------------------------------------------------


#' Get Breaks for B-Spline Basis
#'
#' Computes the breakpoints (knots) needed to construct a B-spline basis given a vector of timepoints, the spline order, and the desired number of internal knots.
#'
#' @param t A numeric vector of timepoints over which the B-spline basis will be constructed.
#' @param k The order of the B-spline
#' @param m The number of internal knots to include.
#'
#' @return A numeric vector of knot positions, including boundary and internal knots, suitable for use in B-spline basis construction.
#'
#' @export

get_knots <- function(t, k, m) {
  # external knots are on boundary
  # return boundary with internal knots only
  breaks <- c(min(t), seq(from = min(t), to = max(t), length.out = m + 2)[-c(1, m + 2)], max(t))
  return(breaks)
}



#' Construct \eqn{d}th-Order Difference Operator Matrix
#'
#' Returns the difference operator matrix \eqn{D} used to penalize differences in coefficient values across responses.
#' **Note:** Currently implemented only for \eqn{d = 2}.
#'
#' @param K The number of response variables.
#' @param d The order of the difference operator (currently only \code{d = 2} is supported).
#' @param order The order of the B-spline basis functions.
#' @param nknots The number of internal knots used in the B-spline basis.
#'
#' @return A matrix \code{D} representing the second-order difference operator applied across response curves.
#'
#' @export
get_D <- function(K, d, order, nknots) {
  # P <- order + nknots
  C <- matrix(0, nrow = nknots + order - d, ncol = nknots + order)
  # dth order weighted difference operator
  for (j in 1:(nknots + order - 2)) {
    d_j <- c(rep(0, j - 1), 1, -2, 1, rep(0, (nknots + order) - 3 - (j - 1)))
    e_j <- c(rep(0, j - 1), 1, rep(0, (nknots + order) - 3 - (j - 1)))
    #C <- C + e_j %*% t(d_j)
    C <- C + tcrossprod(e_j, d_j)
  }
  #D <- t(C) %*% C
  D <- crossprod(C)
  diagD <- kronecker(diag(K), D)

  return(diagD)
}



#' Construct A Matrix for Pairwise Differences
#'
#' Constructs the \code{A} matrix used to track and apply pairwise differences across response curves. This matrix is used in the clustering penalty formulation.
#'
#' @param Kappa A matrix or index structure that encodes all unique response pairs to compute pairwise differences.
#' @param K The number of response variables.
#' @param P The number of B-spline coefficients (typically \code{order + nknots}).
#'
#' @return Matrix \code{A} that maps pairwise differences between responses in the model formulation.
#'
#' @export
get_A <- function(Kappa, K, P) {
  I <- diag(P)
  A <- matrix(ncol = K * P, nrow = P * nrow(Kappa))
  for (kappa in 1:nrow(Kappa)) {
    k1 <- Kappa[kappa, 1]
    k2 <- Kappa[kappa, 2]

    e_k1 <- rep(0, K)
    e_k2 <- rep(0, K)
    e_k1[k1] <- 1
    e_k2[k2] <- 1

    A_kappa <- kronecker(t(e_k1 - e_k2), I)

    A[(P * (kappa - 1) + 1):(P * kappa), ] <- A_kappa
  }

  return(A)
}


# Return clusters ---------------------------------------------------------
#' Get Cluster Membership from Pairwise Differences
#'
#' Derives cluster assignments based on the final pairwise difference estimates (\code{v}) from Algorithm 1.
#'
#' @param v A matrix of final pairwise differences from Algorithm 1, used to infer which subjects belong to the same cluster.
#' @param K The number of response variables.
#' @param P The number of B-spline coefficients (equal to \code{order + nknots}).
#'
#' @return A list containing:
#' \describe{
#'   \item{clusters}{A vector of cluster membership assignments.}
#'   \item{indices}{A list of indices corresponding to each cluster.}
#'   \item{n_clusters}{The total number of unique clusters identified.}
#' }
#'
#' @export

get_clusters <- function(v, K, P) {

  # Convert v into a matrix of P x length
  v_mat <- matrix(v, nrow = P)
  differences <- apply(v_mat, 2, function(x) {
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



# Calculate y_hat-------------------------------------------------
#' Estimate \eqn{\hat{y}} Given \eqn{\beta}
#'
#' Computes fitted values (\eqn{\hat{y}}) using estimated coefficients \eqn{\beta}, a B-spline basis matrix, and covariate information.
#'
#' @param beta A matrix of estimated coefficients (\eqn{\hat{\beta}}) with dimensions \eqn{(P \cdot K) \times L}, where \eqn{P} is the number of basis functions, \eqn{K} is the number of responses, and \eqn{L} is the number of covariates (excluding intercept).
#' @param B A B-spline basis matrix of dimensions \eqn{N \times P}, where \eqn{N} is the total number of observations and \eqn{P} is the number of spline basis functions.
#' @param Z A matrix of covariates of dimension \eqn{M \times (L + 1)}. The first column must be a vector of 1s (intercept). If no external variables are used, \code{Z} should contain only a single column of 1s.
#' @param K The number of response variables (i.e., the number of columns in \code{Y}).
#' @param Y A matrix of response values (e.g., counts), with \eqn{M} rows and \eqn{K} columns. Each row corresponds to a subject-time combination.
#' @param time A numeric vector of time values (length \eqn{M}) or a column reference if \code{Y} is a data frame.
#' @param baseline A T/F if we want to estimate only the baseline y hats, or the y hats overall
#'
#' @return A numeric vector of fitted values (\eqn{\hat{y}}) of length \eqn{M \cdot K}, representing the estimated responses.
#'
#' @export

estimate_y <- function(beta, B, Z, K, Y, time, baseline = F){

  Z_true <- Z
  if (baseline) {
    Z <- matrix(1, nrow = nrow(Z_true), ncol = 1)
  }

  P <- ncol(B)
  L <- ncol(Z) - 1
  M <- nrow(Z)
  yhat <- matrix(nrow = M, ncol = K)

  for (k in 1:K) {
    y_k <- numeric(M)
    for (l in 0:L) {

      beta_k <- beta[((k - 1)*P + 1):(k*P), l + 1, drop = F]
      Z_l <- diag(Z[,l + 1])
      #y_k <- y_k + Z_l %*% B %*% beta_k
      y_k <- y_k + fast_mat_mult3(Z_l, B, beta_k)
    }
    yhat[,k] <- y_k

  }
  # Turn into RA version:
  yhat_ra <- exp(yhat)/rowSums(exp(yhat))


  y_ra <- Y/rowSums(Y)



  # Check if Z is used, otherwise dont add
  if (identical(Z, matrix(1, nrow = nrow(Z_true), ncol = 1))) {
    Ys <- data.frame(time = time,
                     yhat_ra,
                     y_ra)
    colnames(Ys) <- c("time",
                      paste0("yhat_", 1:K),
                      paste0("y_", 1:K))

  } else {
    Ys <- data.frame(time = time,
                     Z = Z_true[, -1],
                     yhat_ra,
                     y_ra)
    colnames(Ys) <- c("time",
                      "Z",
                      paste0("yhat_", 1:K),
                      paste0("y_", 1:K))
  }

  Ys <- Ys %>%
    tidyr::pivot_longer(
      cols = tidyr::matches("yhat|y"),  # Selects both yhat and y columns
      names_to = c(".value", "response"),
      names_pattern = "(yhat|y)_(\\d+)"  # Splits into two parts: yhat/y and the number
    ) %>%
    dplyr::mutate(response = factor(response, levels = 1:K))

  return(Ys)
}

#' Estimate \eqn{\hat{y}} Given \eqn{\beta}
#'
#' Computes fitted values (\eqn{\hat{y}}) using estimated coefficients \eqn{\beta}, a B-spline basis matrix, and covariate information.
#'
#' @param beta A matrix of estimated coefficients (\eqn{\hat{\beta}}) with dimensions \eqn{(P \cdot K) \times L}, where \eqn{P} is the number of basis functions, \eqn{K} is the number of responses, and \eqn{L} is the number of covariates (excluding intercept).
#' @param B A B-spline basis matrix of dimensions \eqn{N \times P}, where \eqn{N} is the total number of observations and \eqn{P} is the number of spline basis functions.
#' @param Z A matrix of covariates of dimension \eqn{M \times (L + 1)}. The first column must be a vector of 1s (intercept). If no external variables are used, \code{Z} should contain only a single column of 1s.
#' @param K The number of response variables (i.e., the number of columns in \code{Y}).
#' @param Y A matrix of response values (e.g., counts), with \eqn{M} rows and \eqn{K} columns. Each row corresponds to a subject-time combination.
#' @param time A numeric vector of time values (length \eqn{M}) or a column reference if \code{Y} is a data frame.
#'
#' @return A numeric vector of fitted values (\eqn{\hat{y}}) of length \eqn{M \cdot K}, representing the estimated responses.
#'
#' @export

estimate_y_counts <- function(beta, B, Z, K, Y, time){
  P <- ncol(B)
  L <- ncol(Z) - 1
  M <- nrow(Z)
  yhat <- matrix(nrow = M, ncol = K)
  for (k in 1:K) {
    y_k <- numeric(M)
    for (l in 0:L) {

      beta_k <- beta[((k - 1)*P + 1):(k*P), l + 1, drop = F]
      Z_l <- diag(Z[,l + 1])
      #y_k <- y_k + Z_l %*% B %*% beta_k
      y_k <- y_k + fast_mat_mult3(Z_l, B, beta_k)
    }
    yhat[,k] <- y_k

  }

  # In case any negatives:
  yhat[yhat < 0] <- 0

  return(yhat)
}


# Re-fit group beta estimates ---------------------------------------------


beta_cluster_group <- function(y, Z, beta, lp,  lp_minus, B, clusters, K, P, M) {

  L <- ncol(Z) - 1
  # Initialize matrix to hold coefficient estimates
  num_features <- ncol(B)
  num_samples <- ncol(y)
  beta_group <- matrix(nrow = num_features, ncol = num_samples)

  log_y <- log(y + 1)

  y_minus_Z_B_beta <- matrix(nrow = M, ncol = K)
  for (k in 1:K) {
    Z_B_beta <- numeric(M)
    for (l in lp_minus) {
      beta_k <- beta[((k - 1)*P + 1):(k*P), l + 1, drop = F]
      Z_l <- diag(Z[,l + 1])
      #y_k <- y_k + Z_l %*% B %*% beta_k
      Z_B_beta <- Z_B_beta + fast_mat_mult3(Z_l, B, beta_k)
    }
    y_minus_Z_B_beta[,k] <- log_y[,k] - exp(Z_B_beta)
    #y_minus_Z_B_beta[,k] <- y[,k]
  }

  y_minus_Z_B_beta[y_minus_Z_B_beta < 0] <- 0

  # Log-transform the response
  log_y_ZB <-  y_minus_Z_B_beta


  # Loop over each unique cluster
  unique_clusters <- unique(clusters$membership)

  # transform beta for easier indexing.
  # Reshape beta_in: (P*K) x L --> (P, K, L)
  beta_array <- array(beta, dim = c(P, K, L + 1))
  # Permute to (P, L, K)
  beta_array_perm <- aperm(beta_array, c(1, 3, 2))
  # Reshape to matrix: (P*L) x K
  beta_group <- matrix(beta_array_perm, nrow = P * (L + 1), ncol = K)


  for (i in seq_along(clusters$csize)) {
    current_cluster <- unique_clusters[i]

    # Get indices of samples that belong to the current cluster
    class_indices <- which(clusters$membership == current_cluster)

    # Prepare design matrix and response vector for current cluster
    B_class <- do.call(rbind, replicate(length(class_indices), B, simplify = FALSE))
    log_y_class <- as.vector(log_y_ZB[, class_indices])

    # Estimate coefficients using least squares (for log(y + 1))
    beta_hat <- MASS::ginv(t(B_class) %*% B_class) %*% t(B_class) %*% log_y_class

    # Update beta with the estimated group beta values
    beta_group[(lp * P + 1):((lp + 1) * P), class_indices] <- beta_hat
  }


  # Step 1: Reshape to array of dimensions (P, L, K)
  beta_array_back <- array(beta_group, dim = c(P, (L + 1), K))

  # Step 2: Permute dimensions back to (P, K, L)
  beta_array_orig <- aperm(beta_array_back, c(1, 3, 2))

  # Step 3: Flatten to matrix of shape (P*K) x L
  beta_in <- matrix(beta_array_orig, nrow = P * K, ncol = (L + 1))




  return(beta_in)
}



beta_cluster_group_update <- function(y, Z, beta, lp,  lp_minus, B, clusters, K, P, M) {

  L <- ncol(Z) - 1

  # Loop over each unique cluster
  unique_clusters <- unique(clusters$membership)

  beta_group <- matrix(nrow = P * (L + 1), ncol = K)

  ZB_list <- list()
  for (l in 0:L) {
    Z_l <- diag(Z[,l + 1])
    ZB_list[[l + 1]] <- Z_l %*% B
  }


  for (c in seq_along(clusters$csize)) {
    current_cluster <- unique_clusters[c]
    class_indices <- which(clusters$membership == current_cluster)
    cluster_size <- length(class_indices)

    df_list <- list()
    beta_names <- c()
    for (l in 0:L) {
      I_k <- diag(cluster_size)
      one_k <- matrix(1, nrow = cluster_size)

      if (l == lp) {
        df_list[[l + 1]] <- kronecker(one_k, ZB_list[[l + 1]])
      } else {
        df_list[[l + 1]] <- kronecker(I_k, ZB_list[[l + 1]])
      }
    }
    df <- do.call(cbind, df_list)

    log_y_class <- as.vector(log(y[, class_indices] + 1))



    mod_g <- lm(log_y_class ~ 0 + df)
    coefs <- coef(mod_g)
    beta_c_mat <- matrix(coefs, nrow = P)



    col_index <- 1

    for (l in 0:L) {
      if (l == lp) {
        # Shared coefficients across all in cluster
        beta_group[(lp * P + 1):((lp + 1) * P), class_indices] <- beta_c_mat[, col_index]
        col_index <- col_index + 1
      } else {
        for (k in seq_along(class_indices)) {
          col_k <- col_index
          idx_k <- class_indices[k]
          beta_group[(l * P + 1):((l + 1) * P), idx_k] <- beta_c_mat[, col_k]
          col_index <- col_index + 1
        }
      }
    }

  }


  # Convert back to original structure
  beta_array_back <- array(beta_group, dim = c(P, (L + 1), K))
  beta_array_orig <- aperm(beta_array_back, c(1, 3, 2))
  beta_in <- matrix(beta_array_orig, nrow = P * K, ncol = (L + 1))

  return(beta_in)
}


