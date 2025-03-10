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

  # Initialize phi to be 1
  #Dirichlet Multinomial over dispersion parameter.
  phi <- 1

  # Initialize beta vector to 0s.
  # beta <- initialize_beta_count(K = K, L = L,
  #                               P = P, B = B,
  #                               Y = Y, Z = Z)
  # Initialize using Algorithm 2 for ALL responses
  zeros_beta <- matrix(0, nrow = K * P, ncol = L + 1)
  beta <- algorithm2(Y = Y,
                     Z = Z,
                     mi_vec = mi_vec,
                     i_index = i_index,
                     lp = -1,
                     B = B,
                     beta = zeros_beta,
                     D = D,
                     gammas = gammas,
                     phi = phi,
                     C = C)
  beta <- zeros_beta
  # Initialize lambda to all 0
  lambda <- numeric(nrow(Kappa)*P)
  # loop initialization
  error <- NULL
  loop_list_beta <- list()
  loop_list_diff <- list()
  admm_beta_list <- list()
  admm_diffs <- list()
  phis_list <- list()
  cluster_list <- list()
  r_list <- list()
  d_list <- list()
  diff <- Inf
  s <- 1


  # hyperparameters <- get_hyperparameters(
  #   Y = Y,
  #   Z, mi_vec, lp, B, D, K, P)

  # pb_outer <- utils::txtProgressBar(min = 0,
  #                             max = max_outer_iter,
  #                             style = 3)



  # pb_outer <- progress::progress_bar$new(
  #   format = "Outer Loop [:bar] :percent (:current/:total)",
  #   total = max_outer_iter, clear = FALSE, width = 50
  # )
  #
  # pb_admm <- progress::progress_bar$new(
  #   format = "  Inner Loop [:bar] :percent (:current/:total)",
  #   total = max_admm_iter, clear = TRUE, width = 50
  # )

  #pb_outer$tick(0)
  while ((diff > epsilon_b) & (s <= max_outer_iter)) {
    cat(paste0("\n","Backfitting algorithm iteration: ", s, "\n"))
    # Go straight to ADMM code if there is no external variables
    #cat(sprintf("\nOuter Loop (s = %d):\n", s))
    #setTxtProgressBar(pb_outer, s)
    if (L > 0) {
      # Solution for l_p minus, with l_p fixed

      beta_lp_minus <- tryCatch({
        algorithm2(Y = Y,
                   Z = Z,
                   mi_vec = mi_vec,
                   i_index = i_index,
                   lp = lp,
                   B = B,
                   beta = beta,
                   D = D,
                   gammas = gammas,
                   phi = phi,
                   C = C)
      }, error = function(e) {
        print(paste0("ERROR!!!!: ", e$message))
        return(e$message)
      })
    } else {
      beta_lp_minus <- beta
    }

    if (is.character(beta_lp_minus)) {
      error <- beta_lp_minus
      break
    }

    # Solution for l p (ADMM) with beta l_p minus fixed
    alg3 <- tryCatch({
      algorithm3(Y = Y,
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
                 max_admm_iter = max_admm_iter)
    }, error = function(e) {
      print(paste0("ERROR!!!!: ", e$message))
      return(e$message)
    })
    if (is.character(alg3)) {
      error <- alg3
      break
    }
    # phi <- estimate_phi()
    beta_lp <- alg3$beta
    admm_beta_list[[s]] <- alg3$beta_admm_track
    lambda <- alg3$lambda
    v <- alg3$v

    # Difference in betas between this loop and the last
    diff <- sum(abs(beta_lp - beta)) # matrix difference
    beta <- beta_lp

    loop_list_beta[[s]] <- beta
    loop_list_diff[[s]] <- diff
    admm_diffs[[s]] <- alg3$diff_admm
    phis_list[[s]] <- alg3$phi_track
    r_list[[s]] <- alg3$r_list
    d_list[[s]] <- alg3$d_list
    cluster_list[[s]] <- alg3$cluster_list


    # Increment loop
    s <- s + 1

    # pb_outer$tick()
    # pb_admm <- progress::progress_bar$new(  # Reset inner progress bar
    #   format = "  Inner Loop [:bar] :percent (:current/:total)",
    #   total = max_admm_iter, clear = TRUE, width = 50
    # )
  }
  if (!is.list(alg3)) {
    stop("Error in algorithm 3")
  }


  # After loop:
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


  BIC <- BIC_cluster(y_ra_df= y_hat,
                     K = K,
                     n_clusters = clusters$no,
                     mi_vec = mi_vec,
                     nknots = nknots,
                     order = order)

  return(list(beta = beta,
              clusters = clusters,
              y_hat = y_hat,
              v = v,
              B = B,
              admm_diffs = admm_diffs,
              admm_beta_list = admm_beta_list,
              loop_list_beta = loop_list_beta,
              loop_list_diff = loop_list_diff,
              cluster_list = cluster_list,
              phis_list = phis_list,
              r_list = r_list,
              d_list = d_list,
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
#' @param Y
#' @param time
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

      beta_k <- beta[((k - 1)*P + 1):(k*P), l + 1]
      Z_l <- diag(Z[,l + 1])
      y_k <- y_k + Z_l %*% B %*% beta_k
    }
    yhat[,k] <- y_k

  }
  # Turn into RA version:
  yhat_ra <- exp(yhat)/rowSums(exp(yhat))

  y_ra <- Y/rowSums(Y)


  Ys <- data.frame(time = time,
                         Z = Z[,2],
                         yhat_ra,
                         y_ra)

  colnames(Ys) <- c("time",
                    "Z",
                    paste0("yhat_", 1:K),
                    paste0("y_", 1:K))


  Ys %>%
    tidyr::pivot_longer(cols = -c(Z, time),
                        names_to = c("type", "response"),
                        names_sep = "_",
                        values_to = "value")


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



# Hyperparameter choosing -------------------------------------------------


