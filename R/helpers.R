# Various -----------------------------------------------------------------

#' Title
#'
#' @param i subject index
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns indeces for each ith subject repeated mi times
#' @export
get_indeces <- function(i, mi_vec) {
  start <- c(0, cumsum(mi_vec))[i] + 1
  end <- cumsum(mi_vec)[i]
  return(id = start:end)
}



#' Get \eqn{Y_{ij0}}
#'
#' @param i subject index
#' @param j time index
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
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
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param subject_ids either a vector of length(Y) or a column reference if Y is a data frame
#' @param time_ids either a vector of length(Y) or a column reference if Y is a data frame. Must be numeric
#'
#' @returns data frame of mi and subject id
#' @export
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
get_Y_wrapper <- function(Y,
                          subject_ids,
                          time_ids) {

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
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns Matrix of dimension mi x K
#' @export
get_Y_i_mat <- function(i, mi_vec, Y) {
  indeces <- get_indeces(i, mi_vec)
  Y_i <- Y[indeces, ]
  # Y_i <- Y[((i - 1) * mi + 1):(i * mi), ]
  # removing rownames since that takes a long time!
  #rownames(Y_i) <- paste0("j=", 1:mi_vec[i])
  #colnames(Y_i) <- paste0("i=", i, "k=", 1:ncol(Y))
  return(Y_i)
}



#' Get Yi as a vector
#'
#' different k responses are listed after each other. so it goes
#' y(i,j=1,k=1), ... yi(j = mi, k = 1), yi(j = 1, k = 2)...
#' NOPE! Other way around...
#'
#' @param i subject index
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param Y_mat (optional) pre-calculated Yi-matrix from Y
#'
#' @returns Vector of length K * mi
#' @export
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
  # removing rownames since that takes a long time!
  # names(Y_i_vec) <- paste0("i=", i, ",j=", js, ",k=", ks)
  return(Y_i_vec)
}


#' Get \eqn{Z_{ijl}}
#'
#' @param i subject index
#' @param j time index
#' @param l external variable index
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns Scalar for the l-th external variable for subject i at time j
#' @export
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
#' @param P Number of B-spline coefficients (order + nknots)
#'
#' @returns Vector of length P
#' @export
get_beta_kl <- function(k, l, beta, P) {
  beta_kl <- beta[((k - 1) * P + 1):(k * P), (l + 1)]
  return(beta_kl)
}




# B spline basis ----------------------------------------------------------
#' Get \eqn{B_ij}
#'
#' @param i subject index
#' @param j time index
#' @param B B spline basis matrix of dimension (N x P)
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns A vector of length P
#' @export
get_B_ij <- function(i, j, B, mi_vec) {
  i_start <- c(0, cumsum(mi_vec))[i] + 1
  i_start + j - 1
  B_ij <- B[i_start + j - 1, ]
  # removing rownames since that takes a long time!
  # names(B_ij) <- paste0("i=", i, ",j=", j, ",l=", 1:ncol(B))
  return(B_ij)
}


#' Title
#'
#' @param time vector of time values for each subject/time
#' @param order Order of the B-spline basis
#' @param nknots Number of knots for the B-spline basis
#'
#' @returns B spline basic matrix of dimension (M x P)
#' @export
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
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param B B spline basis matrix of dimension (N x P)
#'
#' @returns Numeric (1x1)
#' @export
get_alpha_ijk <- function(i, j, k, beta, Z, B, mi_vec) {
  P <- ncol(B)
  L <- ncol(Z) - 1
  lsum <- numeric(L)
  B_ij <- get_B_ij(i, j, B, mi_vec)

  for (l in 0:L) {
    Z_ijl <- get_Z_ijl(i, j, l, Z, mi_vec)
    beta_lk <- get_beta_kl(k, l, beta, P)
    #lsum[l + 1] <- Z_ijl * t(B_ij) %*% beta_lk
    lsum[l + 1] <- Z_ijl * crossprod(B_ij, beta_lk)
  }
  alpha_ijk <- exp(sum(lsum))
  if (any(is.infinite(alpha_ijk))) {
    stop("Infinite alpha")
  }
  # removing rownames since that takes a long time!
  # names(alpha_ijk) <- paste0("i=", i, ",j=", j, ",k=", k)

  if (sum(alpha_ijk < 0) != 0) {
    stop("Negative alpha")
  }
  return(alpha_ijk)
}


#' Calculate alpha for all k for a single i and j.
#'
#' @param i subject index
#' @param j time index
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param K Number of responses
#'
#' @returns Vector of length K
#' @export
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
#' @param Y_ij0 Total sum of counts across all K
#' @param i subject index
#' @param j time index
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of responses
#' @param alpha_ij Vector of DM parameters for i, j
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns Vector of length K
#' @export
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
    # removing rownames since that takes a long time!
    # names(mu_ij) <- paste0("i=", i, ",j=", j, ",k=", 1:K)
  } else {
    mu_ij <- Y_ij0 / sum(alpha_ij) * alpha_ij
    names(mu_ij) <- names(alpha_ij)
  }
  return(mu_ij)
}


#' Get mu_i mean vector for ith sample
#'
#' @param i subject index
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of responses
#' @returns Vector of length K*mi with each K first
#' @export
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
#' @param alpha_ij Vector of DM parameters for i, j
#' @param i subject index
#' @param j time index
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param K Number of responses
#'
#' @returns U ij part of partials/variance matrix
#' @export
get_U_ij <- function(alpha_ij, i, j, beta, Z, B, K, mi_vec) {
  if (missing(alpha_ij)) {
    # calculate with i, j, beta, Z, B
    alpha_ij <- get_alpha_ij(i, j, beta, Z, B, K, mi_vec)
  }
  alpha_ij0 <- sum(alpha_ij)

  #U_ij <- diag(alpha_ij / alpha_ij0) - alpha_ij %*% t(alpha_ij) / alpha_ij0^2
  U_ij <- diag(alpha_ij / alpha_ij0) - tcrossprod(alpha_ij) / alpha_ij0^2
  return(U_ij)
}


#' Get DM variance for ith subject jth timepoint
#'
#' @param Y_ij0 Total sum of counts across all K
#' @param phi Current value of overdispersion parameter
#' @param alpha_ij Vector of DM parameters for i, j
#' @param i subject index
#' @param j time index
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of responses
#' @param mi number of timepoints for ith subject
#'
#' @returns matrix of dimension K x K for a single j
#' @export
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
#' @param i subject index
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param phi Current value of overdispersion parameter
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param Y_ij0 Total sum of counts across all K
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param K Number of responses
#'
#' @returns Diagonal matrix of Kmi x Kmi Each block is for a single j of block size K x K
#' @export
get_V_i <- function(i, Y, Y_ij0, phi, beta, Z, B, K, mi_vec) {
  V_ij_list <- list()
  if (missing(Y_ij0)) {
    Y_ij0 <- get_Y_ij0(i, 1, Y, mi_vec)
  }
  mi <- mi_vec[i]
  for (j in 1:mi) {
    Y_ij0 <- get_Y_ij0(i, j, Y, mi_vec)
    alpha_ij <- get_alpha_ij(
      i = i,
      j = j,
      beta = beta,
      Z = Z,
      B = B,
      K = K,
      mi_vec = mi_vec
    )
    V_ijj <- get_V_ijj(
      Y_ij0 = Y_ij0,
      phi = phi,
      alpha_ij = alpha_ij
    )
    V_ij_list[[j]] <- V_ijj
  }
  V_i_bdiag <- Matrix::bdiag(V_ij_list)
  V_i <- as.matrix(V_i_bdiag)
  if (any(is.nan(V_i))) {
    browser()
  }
  return(V_i)
}

