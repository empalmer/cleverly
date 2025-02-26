# Various -----------------------------------------------------------------

get_indeces <- function(i, mi_vec) {
  start <- c(0, cumsum(mi_vec))[i] + 1
  end <- cumsum(mi_vec)[i]
  return(id = start:end)
}





#' Get \eqn{Y_{ij0}}
#'
#' @param i subject
#' @param j time index
#' @param Y response matrix of dimension N x K
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns Scalar of the total sum constraint for a given i, j
#' @export
#'
get_Y_ij0 <- function(i, j, Y, mi_vec) {
  i_start <- c(0, cumsum(mi_vec))[i] + 1
  Y_ij0 <- sum(Y[i_start + j - 1, ])
  # Y_ij0 <- sum(Y[((i - 1) * mi + j), ])
  return(Y_ij0)
}


#' Get data frame of mi values for each i
#'
#' @param Y
#' @param subject_ids
#' @param time_ids
#'
#' @returns
#' @export
#'
#' @examples
get_mi_vec <- function(Y, subject_ids, time_ids) {
  Y_wrapper <- get_Y_wrapper(Y, subject_ids, time_ids)

  mi_vec <- data.frame(subject_id = Y_wrapper$subject_id_values) %>%
    dplyr::group_by(subject_id) %>%
    dplyr::summarise(mi = dplyr::n())

  return(mi_vec)
}


#' Get Yi for multiple user input options
#'
#' Users can either specify Y as a matrix or data frame.
#'
#' The subject ids and time ids can either be references to columns in Y or external vectors that link the ids to the rows of Y
#'
#' @param Y either a data frame or matrix. Each response should be a separate column. Each row should be a separate subject/time combination. There should be M total rows.
#' @param subject_ids either a vector of length(Y) or a column reference if Y is a data frame
#' @param time_ids either a vector of length(Y) or a column reference if Y is a data frame
#'
#' @returns Matrix of dimension mi x K
#' @export
#'
#' @examples
get_Y_wrapper <- function(Y, subject_ids, time_ids) {
  id_quo <- rlang::enquo(subject_ids)
  time_quo <- rlang::enquo(time_ids)

  if (is.data.frame(Y)) {
    # For ID:
    if (rlang::quo_name(id_quo) %in% colnames(Y)) {
      subject_id_values <- dplyr::pull(Y, !!id_quo) # Extract the column
      Y <- dplyr::select(Y, -!!id_quo) # Remove the column
    } else if (length(subject_ids) == nrow(Y)) {
      subject_id_values <- subject_ids # Assume it's an external vector
    } else {
      stop("Invalid ID: must be a column name in Y or a vector of length nrow(Y)")
    }

    # For time:
    if (rlang::quo_name(time_quo) %in% colnames(Y)) {
      time_id_values <- dplyr::pull(Y, !!time_quo) # Extract the column
      Y <- dplyr::select(Y, -!!time_quo) # Remove the column
    } else if (length(time_ids) == nrow(Y)) {
      time_id_values <- time_ids # Assume it's an external vector
    } else {
      stop("Invalid ID: must be a column name in Y or a vector of length nrow(Y)")
    }

    # Convert Y into a matrix
    Y <- as.matrix(Y)
    dimnames(Y) <- NULL # To make the tests pass
  }
  # If Y is a matrix, expect ID and Time to be separate vectors
  else if (is.matrix(Y)) {
    if (length(subject_ids) != nrow(Y) || length(time_ids) != nrow(Y)) {
      stop("For matrices, ID and Time must be external vectors of length nrow(Y)")
    }
    subject_id_values <- subject_ids
    time_id_values <- time_ids
  } else {
    stop("Y must be either a data frame or a matrix")
  }


  list <- list(
    Y = Y,
    subject_id_values = subject_id_values,
    time_id_values = time_id_values
  )

  return(list)
}


#' Get \eqn{Y_i}
#'
#' @param i subject index
#' @param Y response matrix of dimension N x K
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns Matrix of dimension mi x K
#' @export
#'
#' @examples
get_Y_i_mat <- function(i, mi_vec, Y) {
  indeces <- get_indeces(i, mi_vec)
  Y_i <- Y[indeces, ]
  # Y_i <- Y[((i - 1) * mi + 1):(i * mi), ]
  rownames(Y_i) <- paste0("j=", 1:mi_vec[i])
  colnames(Y_i) <- paste0("i=", i, "k=", 1:ncol(Y))
  return(Y_i)
}



#' Get Yi as a vector
#'
#' different k responses are listed after each other. so it goes
#' y(i,j=1,k=1), ... yi(j = mi, k = 1), yi(j = 1, k = 2)...
#' NOPE! Other way around...
#'
#' @param i subject index
#' @param Y response matrix of dimension N x K
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param Y_mat (optional) pre-calculated Yi-matrix from Y
#'
#' @returns Vector of length K * mi
#' @export
#'
#' @examples
get_Y_i_vec <- function(i, mi_vec, Y, Y_mat) {
  if (missing(Y_mat)) {
    Y_i_mat <- get_Y_i_mat(
      i = i,
      mi_vec = mi_vec,
      Y = Y
    )
  }

  # This goes row by row.
  Y_i_vec <- as.vector(t(Y_i_mat))
  js <- rep(1:mi_vec[i], each = ncol(Y))
  ks <- rep(1:ncol(Y), mi_vec[i])
  names(Y_i_vec) <- paste0("i=", i, ",j=", js, ",k=", ks)
  return(Y_i_vec)
}







#' Get \eqn{Z_{ijl}}
#'
#' @param i subject index
#' @param j time index
#' @param l external variable index
#' @param Z matrix of external variables of dimension N x L
#'
#' @returns Scalar for the l-th external variable for subject i at time j
#' @export
#'
#' @examples
get_Z_ijl <- function(i, j, l, Z, mi_vec) {
  i_start <- c(0, cumsum(mi_vec))[i] + 1
  Z_ijl <- Z[i_start + j - 1, (l + 1)]
  return(Z_ijl)
}



#' Get \eqn{\beta_{kl}}
#'
#' @param k response index
#' @param l external variable index
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param P number of B-spline coefficients
#'
#' @returns Vector of length P
#' @export
#'
#' @examples
get_beta_kl <- function(k, l, beta, P) {
  beta_kl <- beta[((k - 1) * P + 1):(k * P), (l + 1)]
  return(beta_kl)
}




# B spline basis ----------------------------------------------------------
#' Get \eqn{B_ij}
#'
#' @param i subject index
#' @param j time index
#' @param B B-spline basis matrix of dimension N x P
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns A vector of length P
#' @export
#'
get_B_ij <- function(i, j, B, mi_vec) {
  i_start <- c(0, cumsum(mi_vec))[i] + 1
  i_start + j - 1
  B_ij <- B[i_start + j - 1, ]
  names(B_ij) <- paste0("i=", i, ",j=", j, ",l=", 1:ncol(B))
  return(B_ij)
}


#' Title
#'
#' @param time
#' @param order
#' @param nknots
#'
#' @returns
#' @export
#'
#' @examples
get_B <- function(time, order, nknots) {
  B <- fda::bsplineS(time,
    get_knots(time, k = order, m = nknots),
    norder = order
  )
  return(B)
}

# Alphas ------------------------------------------------------------------

#' Calculate alpha for a single \eqn{i = 1, j = 1, k = 1} from \eqn{\beta}
#'
#' @param i subject index
#' @param j time index
#' @param k response index
#' @param beta matrix of dimension (P*K) x L
#' @param Z matrix of external variables of dimension N x L
#' @param B B spline basis matrix of dimension (N x P)
#' @param mi Number of timepoints for subject i
#'
#' @returns Numeric (1x1)
#' @export
#'
#' @examples
get_alpha_ijk <- function(i, j, k, beta, Z, B, mi_vec) {
  P <- ncol(B)
  L <- ncol(Z) - 1
  lsum <- numeric(L)
  B_ij <- get_B_ij(i, j, B, mi_vec)

  for (l in 0:L) {
    Z_ijl <- get_Z_ijl(i, j, l, Z, mi_vec)
    beta_lk <- get_beta_kl(k, l, beta, P)
    lsum[l + 1] <- Z_ijl * t(B_ij) %*% beta_lk
  }
  alpha_ijk <- exp(sum(lsum))
  names(alpha_ijk) <- paste0("i=", i, ",j=", j, ",k=", k)

  if (sum(alpha_ijk < 0) != 0) {
    stop("Negative alpha")
  }
  return(alpha_ijk)
}


#' Calculate alpha for all k for a single i and j.
#'
#' @param i subject index
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
get_alpha_ij <- function(i, j, beta, Z, B, K, mi_vec) {
  alphas <- c()
  for (k in 1:K) {
    alphas <- c(alphas, get_alpha_ijk(i, j, k, beta, Z, B, mi_vec))
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
get_mu_ij <- function(Y_ij0, alpha_ij, i, j, beta, Z, B, K, mi_vec) {
  if (missing(alpha_ij)) {
    # calculate with i, j, beta, Z, B
    alpha_ij <- get_alpha_ij(
      i = i,
      j = j,
      beta = beta,
      Z = Z,
      B = B,
      K = K,
      mi_vec = mi_vec
    )
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
get_mu_i <- function(i, mi_vec, Y, beta, Z, B, K) {
  mu_i <- c()
  K <- ncol(Y)
  mi <- mi_vec[i]
  for (j in 1:mi) {
    Y_ij0 <- get_Y_ij0(i, j, Y, mi_vec)
    mu_ij <- get_mu_ij(
      Y_ij0 = Y_ij0,
      i = i,
      j = j,
      beta = beta,
      Z = Z,
      B = B,
      K = K,
      mi_vec = mi_vec
    )
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
get_U_ij <- function(alpha_ij, i, j, beta, Z, B, K, mi_vec) {
  if (missing(alpha_ij)) {
    # calculate with i, j, beta, Z, B
    alpha_ij <- get_alpha_ij(i, j, beta, Z, B, K, mi_vec)
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
#' @returns matrix of dimension K x K for a single j
#' @export
#'
#' @examples
get_V_ijj <- function(Y_ij0, phi, alpha_ij, i, j, beta, Z, B, K, mi) {
  if (missing(alpha_ij)) {
    # calculate with i, j, beta, Z, B
    alpha_ij <- get_alpha_ij(i, j, beta, Z, B, K, mi)
  }
  U_ij <- get_U_ij(alpha_ij)

  V_ijj <- phi * Y_ij0 * U_ij
  return(V_ijj)
}



#' Get Vi diagonal matrix of all js.
#'
#' @param i
#' @param Y
#' @param phi
#' @param beta
#' @param Z
#' @param B
#' @param K
#' @param mi
#'
#' @returns Diagonal matrix of Kmi x Kmi Each block is for a single j of block size K x K
#' @export
#'
#' @examples
get_V_i <- function(i, Y, Y_ij0, phi, beta, Z, B, K, mi_vec) {
  V_ij_list <- list()
  if (missing(Y_ij0)) {
    Y_ij0 <- get_Y_ij0(i, 1, Y, mi_vec)
  }
  mi <- mi_vec[i]
  for (j in 1:mi) {
    Y_ij0 <- get_Y_ij0(i, j, Y, mi_vec)
    V_ijj <- get_V_ijj(
      Y_ij0 = Y_ij0,
      phi = phi,
      alpha_ij = get_alpha_ij(
        i = i,
        j = j,
        beta = beta,
        Z = Z,
        B = B,
        K = K,
        mi_vec = mi_vec
      )
    )
    V_ij_list[[j]] <- V_ijj
  }
  V_i_bdiag <- Matrix::bdiag(V_ij_list)
  V_i <- as.matrix(V_i_bdiag)
  return(V_i)
}


# Initializing algorithm 1 ------------------------------------------------


#' Initialize beta
#'
#' @param K Number of responses
#' @param L Number of external variables
#' @param P Number of b-spline coefficients
#'
#' @returns 0 matrix of dimension KP x (L + 1)
#' @export
#'
#' @examples
initialize_beta <- function(K, L, P) {
  beta_init <- matrix(0, nrow = K * P, ncol = L + 1)
}




#' Format Z
#'
#' @param Z A matrix or data frame with columns of external variables for each subject/time
#'
#' @returns A matrix with a column of 1s representing L = 0, and values for the other external variables
#' @export
#'
#' @examples
format_Z <- function(Z) {
  M <- nrow(Z)
  if (!identical(Z[, 1], rep(1, M))) {
    Z <- cbind(1, Z)
  }
  return(Z)
}



#' Title
#'
#' @param lp
#' @param Z
#'
#' @returns
#' @export
#'
#' @examples
format_lp <- function(lp, Z) {
  return(lp)
}


#' Get breaks for B-spline basis
#'
#' @param t timepoints
#' @param k order of B-spline
#' @param m number of knots
#'
#' @returns
#' @export
#'
#' @examples
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
#' @param K
#' @param d
#' @param order
#' @param nknots
#'
#' @returns
#' @export
#'
#' @examples
get_D <- function(K, d, order, nknots) {
  P <- order + nknots
  C <- matrix(0, nrow = nknots + order - d, ncol = nknots + order)
  # dth order weighted difference operator
  for (j in 1:(nknots + order - 2)) {
    d_j <- c(rep(0, j - 1), 1, -2, 1, rep(0, (nknots + order) - 3 - (j - 1)))
    e_j <- c(rep(0, j - 1), 1, rep(0, (nknots + order) - 3 - (j - 1)))
    C <- C + e_j %*% t(d_j)
  }
  D <- t(C) %*% C
  diagD <- kronecker(diag(K), D)

  return(diagD)
}



#' Get A matrix
#'
#' This is the matrix that keeps track of the pairwise differences,
#'
#' @param Kappa
#' @param K
#' @param P
#'
#' @returns
#' @export
#'
#' @examples
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
