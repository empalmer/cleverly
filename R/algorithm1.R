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
#' @param max_2_iter Maximum number of iterations for algorithm 2 to run each loop
#' @param epsilon_2 Tolerance for convergence of algorithm 2
#' @param cor_str Type of correlation structure (IND, CON, AR1)
#' @param run_min
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
                       run_min,
                       nknots,
                       order,
                       epsilon_b,
                       epsilon_r,
                       epsilon_d,
                       max_outer_iter,
                       max_admm_iter,
                       max_2_iter,
                       epsilon_2,
                       cor_str) {

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
  j1_j2_list <- lapply(mi_vec, function(mi) {
    times <- 1:mi
    diag_values <- rep(times, each = K)
    outer(diag_values, diag_values, "-")
  })
  # off_bdiag_list <- lapply(j1_j2_list, function(block) {
  #   # return T if off the block diagonal
  #   !(block == 0 | lower.tri(block))
  #   browser()
  # })
  # Define a list with blocks for each i what the correlation matrix structure is
  # We set rho = 1 just to give the 0/1 matrix of zero/nonzero elements of the
  # correlation matrix depending on the chosen strucutre.
  cor_blocks <- purrr::map(mi_vec, ~get_corR(cor_str = cor_str,
                                             mi = .x,
                                             K = K,
                                             rho = 1))

  # Initialize phi to be 1
  #Dirichlet Multinomial over dispersion parameter.
  phi <- 1

  # Initialize using Algorithm 2 for ALL responses
  print("Initializing beta values")
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
                   M = M)



    # Difference in betas between this loop and the last
    diff <- sum(abs(beta - beta_old)) # matrix difference
    # alg1_diff[[s]] <- diff

    # Exit outer loop early if convergence is reached
    if (diff < epsilon_b) {
      break
    }

    cluster_list[[s]] <- alg3$cluster_list
    # Exit if constant cluster results for the past 3 iterations
    if (s >= run_min) {
      current <-  alg3$cluster_list[[length(alg3$cluster_list)]]$membership
      past1 <- cluster_list[[s - 1]][[length(cluster_list[[s - 1]])]]$membership
      past2 <- cluster_list[[s - 2]][[length(cluster_list[[s - 2]])]]$membership

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
  y_hat_init <- estimate_y(beta = beta_init$beta,
                      B = B,
                      Z = Z,
                      K = K,
                      Y = Y,
                      time = time)
  BIC <- BIC_cluster(y_ra_df = y_hat,
                     K = K,
                     n_clusters = clusters$no,
                     mi_vec = mi_vec,
                     nknots = nknots,
                     order = order)

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
              rho = rho,
              phi = phi,
              BIC = BIC,
              error = error))
}





# Helpers for Algorithm 1 -------------------------------------------------
#' Estimate y-hat given beta
#'
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param B B spline basis matrix of dimension (N x P)
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param K Number of responses
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param time either a vector of length(Y) or a column reference if Y is a data frame. Must be numeric
#'
#' @returns Vector of RA values for y
#' @export
estimate_y <- function(beta, B, Z, K, Y, time){
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

  Ys <- data.frame(time = time,
                   Z = Z[,L+1],
                   yhat_ra,
                   y_ra)

  colnames(Ys) <- c("time",
                    "Z",
                    paste0("yhat_", 1:K),
                    paste0("y_", 1:K))


  # Ys %>%
  #   tidyr::pivot_longer(cols = -c(Z, time),
  #                       names_to = c("type", "response"),
  #                       names_sep = "_",
  #                       values_to = "value")


  Ys <- Ys %>%
    tidyr::pivot_longer(
      cols = tidyr::matches("yhat|y"),  # Selects both yhat and y columns
      names_to = c(".value", "response"),
      names_pattern = "(yhat|y)_(\\d+)"  # Splits into two parts: yhat/y and the number
    ) %>%
    dplyr::mutate(response = factor(response, levels = 1:K))


  return(Ys)
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



# Initializing algorithm 1 ------------------------------------------------


#' Initialize beta (for count data)
#'
#' For each response independently, use the least squares estimate using the B-spline basis matrix and log Y
#'
#' @param K Number of responses
#' @param L Number of external variables
#' @param P Number of B-spline coefficients (order + nknots)
#'
#' @returns 0 matrix of dimension KP x (L + 1)
#' @export
initialize_beta_count <- function(K, L, P, B, Y, Z) {
  beta_init <- matrix(0, nrow = K * P, ncol = L + 1)

  # B_0 <- matrix(nrow = ncol(B), ncol = ncol(Y))
  # logy <- log(Y + 1)
  #
  # # Create design matrix that is c(B, diag(Z_l)B) etc
  # X <- do.call(cbind, lapply(1:(L + 1), function(i) diag(Z[, i]) %*% B))
  #
  #
  # # For each response, calculate linear model coefficients
  # for (k in 1:K) {
  #   coefs <- lm(logy[,k] ~ 0 + X)$coefficients
  #   beta_init[((k - 1)*P + 1):(k*P),] <- matrix(coefs, ncol = L + 1)
  # }
  #
  # # Try normalizing?
  #
  # beta_init <- beta_init/sd(beta_init)

  return(beta_init)
}




#' Format Z
#'
#' Add a column of 1s to Z if it doesn't already exist
#'
#' @param Z A matrix or data frame with columns of external variables for each subject/time
#'
#' @returns A matrix with a column of 1s representing L = 0, and values for the other external variables
#' @export
format_Z <- function(Z) {
  if (is.data.frame(Z) | is.matrix(Z)) {
    M <- nrow(Z)
    if (!identical(Z[, 1], rep(1, M))) {
      Z <- cbind(1, Z)
    }
    Z <- as.matrix(Z)
  } else if (is.vector(Z)) {
    if (!all(Z == 1)) {
      Z <- cbind(1, Z)
    }
  }
  return(Z)
}



#' Title
#'
#' @param lp either a numeric index of which external variable to cluster on, or the name of the column of Z that contains the clustering variable. Specify numeric 0 to cluster via baseline.
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#'
#' @returns lp index
#' @export
format_lp <- function(lp, Z) {
  return(lp)
}


#' Get breaks for B-spline basis
#'
#' @param t timepoints
#' @param k order of B-spline
#' @param m number of knots
#'
#' @returns vector of knot values
#' @export
get_knots <- function(t, k, m) {
  # external knots are on boundary
  # return boundary with internal knots only
  breaks <- c(min(t), seq(from = min(t), to = max(t), length.out = m + 2)[-c(1, m + 2)], max(t))
  return(breaks)
}



#' Return dth order differnece operator matrix
#'
#' Note that this only works for d = 2
#'
#' @param K Number of responses
#' @param d Order of the difference operator
#' @param order Order of the B-spline basis
#' @param nknots Number of knots for the B-spline basis
#'
#' @returns D matrix
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



#' Get A matrix
#'
#' This is the matrix that keeps track of the pairwise differences,
#'
#' @param Kappa keeps track of each possible response pair
#' @param K Number of responses
#' @param P Number of B-spline coefficients (order + nknots)
#'
#' @returns A matrix
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
