# Various -----------------------------------------------------------------
#' Get \eqn{Y_{ij0}}
#'
#' @param i subject
#' @param j time index
#' @param Y response
#' @param mi number of time points for subject i
#'
#' @returns Scalar of the total sum constraint for a given i, j
#' @export
#'
get_Y_ij0 <- function(i, j, Y, mi) {
  Y_ij0 <- sum(Y[((i - 1) * mi + j), ])
  return(Y_ij0)
}




#' Get \eqn{Y_i}
#'
#' @param i
#' @param Y
#' @param mi
#'
#' @returns Matrix of dimension mi x K
#' @export
#'
#' @examples
get_Y_i_mat <- function(i, mi, Y){
  Y_i <- Y[((i - 1) * mi + 1):(i * mi), ]
  rownames(Y_i) <- paste0("j=", 1:mi)
  colnames(Y_i) <- paste0("k=", 1:ncol(Y))
  return(Y_i)
}



#' Get Yi as a vector
#'
#' different k responses are listed after each other. so it goes
#' y(i,j=1,k=1), ... yi(j = mi, k = 1), yi(j = 1, k = 2)...
#' NOPE! Other way around...
#'
#' @param i
#' @param Y
#' @param mi
#' @param Y_mat
#'
#' @returns Vector of length K * mi
#' @export
#'
#' @examples
get_Y_i_vec <- function(i, mi, Y, Y_mat){
  if (missing(Y_mat)) {
    Y_i_mat <- get_Y_i_mat(i, mi, Y)
  }

  # This goes row by row.
  Y_i_vec <- as.vector(t(Y_i_mat))
  js <- rep(1:mi, each = ncol(Y))
  ks <- rep(1:ncol(Y),  mi)
  names(Y_i_vec) <- paste0("i=", i, ",j=", js, ",k=", ks)
  return(Y_i_vec)
}


#' Get \eqn{B_ij}
#'
#' @param i subject
#' @param j time index
#' @param B B-spline basis matrix
#' @param mi number of time points for subject i
#'
#' @returns A vector of length P
#' @export
#'
get_B_ij <- function(i, j, B, mi) {
  B_ij <- B[((i - 1) * mi + j), ]
  names(B_ij) <- paste0("i=", i, ",j=", j, ",l=", 1:ncol(B))
  return(B_ij)
}





#' Get \eqn{Z_{ijl}}
#'
#' @param i subject
#' @param j time index
#' @param l external variable index
#' @param Z matrix of external variables
#'
#' @returns Scalar for the lth external variable for subject i at time j
#' @export
#'
#' @examples
get_Z_ijl <- function(i, j, l, Z) {
  Z_ijl <- Z[i + (j - 1), (l + 1)]
  return(Z_ijl)
}



#' Get \eqn{\beta_{kl}}
#'
#' @param k
#' @param l
#' @param beta
#' @param P
#'
#' @returns
#' @export
#'
#' @examples
get_beta_kl <- function(k, l, beta, P) {
  beta_kl <- beta[((k - 1) * P + 1):(k * P), (l + 1)]
  return(beta_kl)
}

# Alphas ------------------------------------------------------------------

#' Calculate alpha for a single \eqn{i = 1, j = 1, k = 1} from \eqn{\beta}
#'
#' @param i subject
#' @param j time index
#' @param k
#' @param beta matrix of dimension (P*K) x L
#' @param Z
#' @param B
#' @param mi Number of timepoints for subject i
#'
#' @returns Numeric
#' @export
#'
#' @examples
get_alpha_ijk <- function(i, j, k, beta, Z, B, mi) {
  P <- ncol(B)
  L <- ncol(Z) - 1
  lsum <- numeric(L)
  for (l in 0:L) {
    Z_ijl <- get_Z_ijl(i, j, l, Z)
    beta_lk <- get_beta_kl(k, l, beta, P)
    B_ij <- get_B_ij(i, j, B, mi)

    lsum[l + 1] <- Z_ijl * t(B_ij) %*% beta_lk
  }
  alpha_ijk <- exp(sum(lsum))
  names(alpha_ijk) <- paste0("i=", i, ",j=", j, ",k=", k)
  return(alpha_ijk)
}


#' Calculate alpha for all k for a single i and j.
#'
#' @param i subject
#' @param j time index
#' @param beta beta matrix of dimension (P*K) x L
#' @param Z matrix of external variables of dimension (N x L)
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of parameters
#' @param mi Number of timepoints for subject i
#'
#' @returns Vector of length K
#' @export
#'
#' @examples
get_alpha_ij <- function(i, j, beta, Z, B, K, mi) {
  alphas <- c()
  for (k in 1:K) {
    alphas <- c(alphas, get_alpha_ijk(i, j, k, beta, Z, B, mi))
  }
  return(alphas)
}




# Mus ------------------------------------------------------------------

#' Get the ijth values for mu
#'
#' @param Y_ij0
#' @param i
#' @param j
#' @param beta
#' @param Z
#' @param B
#' @param K
#' @param mi
#' @param alpha_ij
#'
#' @returns Vector of length K
#' @export
#'
#' @examples
get_mu_ij <- function(Y_ij0, alpha_ij, i, j, beta, Z, B, K, mi) {
  if (missing(alpha_ij)) {
    # calculate with i, j, beta, Z, B
    alpha_ij <- get_alpha_ij(i = i,
                             j = j,
                             beta = beta,
                             Z = Z, B = B, K = K, mi = mi)
    mu_ij <- Y_ij0 / sum(alpha_ij) * alpha_ij
    names(mu_ij) <- paste0("i=", i, ",j=", j, ",k=", 1:K)
  } else {
    mu_ij <- Y_ij0 / sum(alpha_ij) * alpha_ij
    names(mu_ij) <- names(alpha_ij)
  }
  return(mu_ij)
}


#' Get mu_i mean vector for ith sample
#'
#' @param i
#' @param mi
#' @param Y
#' @param beta
#' @param Z
#' @param B
#' @param K
#'
#' @returns Vector of length K*mi with each K first
#' @export
#'
#' @examples
get_mu_i <- function(i, mi, Y, beta, Z, B, K) {
  mu_i <- c()
  K <- ncol(Y)
  for (j in 1:mi) {
    Y_ij0 <- get_Y_ij0(i, j, Y, mi)
    mu_ij <- get_mu_ij(Y_ij0 = Y_ij0, i = i, j = j,
                       beta = beta, Z = Z, B = B, K = K, mi = mi)
    mu_i <- c(mu_i, mu_ij)
  }
  return(mu_i)
}



# Variance -----------------------------------------------------

#' Get term U_ij
#'
#' @param Y_ij0
#' @param alpha_ij
#' @param i
#' @param j
#' @param beta
#' @param Z
#' @param B
#' @param K
#' @param mi
#'
#' @returns
#' @export
#'
#' @examples
get_U_ij <- function(Y_ij0, alpha_ij, i, j, beta, Z, B, K, mi) {
  if (missing(alpha_ij)) {
    # calculate with i, j, beta, Z, B
    alpha_ij <- get_alpha_ij(i, j, beta, Z, B, K, mi)
  }
  alpha_ij0 <- sum(alpha_ij)

  U_ij <- diag(alpha_ij / alpha_ij0) - alpha_ij %*% t(alpha_ij) / alpha_ij0^2
  return(U_ij)
}


#' Get DM variance for ith subject jth timepoint
#'
#' @param Y_ij0
#' @param phi
#' @param alpha_ij
#' @param i
#' @param j
#' @param beta
#' @param Z
#' @param B
#' @param K
#' @param mi
#'
#' @returns
#' @export
#'
#' @examples
get_V_ijj <- function(Y_ij0, phi, alpha_ij, i, j, beta, Z, B, K, mi) {
  if (missing(alpha_ij)) {
    # calculate with i, j, beta, Z, B
    alpha_ij <- get_alpha_ij(i, j, beta, Z, B, K, mi)
  }
  U_ij <- get_U_ij(Y_ij0, alpha_ij)

  V_ijj <- phi * Y_ij0 * U_ij
  return(V_ijj)
}


