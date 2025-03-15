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


#' Get \eqn{Z_{ijl}}
#'
#' @param i subject index
#' @param j time index
#' @param l external variable index
#' @param i_index
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#'
#' @returns Scalar for the l-th external variable for subject i at time j
#' @export
get_Z_ijl <- function(i, j, l, Z, i_index) {
  if (l == 0) {
    return(1)
  } else {
    i_start <- i_index[i] + 1
    Z_ijl <- Z[i_start + j - 1, (l + 1)]
  }
  # i_start <- i_index[i] + 1
  # Z_ijl <- Z[i_start + j - 1, (l + 1)]


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



# Y0 ----------------------------------------------------------------------



#' Get \eqn{Y_{ij0}}
#'
#' @param i subject index
#' @param Y0
#' @param i_index
#' @param j time index
#'
#' @returns Scalar of the total sum constraint for a given i, j
#' @export
#'
get_Y_ij0 <- function(i, j, Y0, i_index) {
  i_start <- i_index[i] + 1
  Y_ij0 <- Y0[i_start + j - 1]
  return(Y_ij0)
}


# Y - response -----------------------------------------------------------


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



# Mus ------------------------------------------------------------------

#' Get the ijth values for mu
#'
#' @param Y_ij0 Total sum of counts across all K
#' @param alpha_ij Vector of DM parameters for i, j
#'
#' @returns Vector of length K
#' @export
get_mu_ij <- function(Y_ij0, alpha_ij) {
    mu_ij <- Y_ij0 / sum(alpha_ij) * alpha_ij
    return(mu_ij)
}


#' Get mu_i mean vector for ith sample
#'
#' @param i subject index
#' @param alpha
#' @param i_index
#' @param Y0
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param K
#'
#' @returns Vector of length K*mi with each K first
#' @export
get_mu_i <- function(i, alpha, mi_vec, i_index, Y0, K) {
  mi <- mi_vec[i]
  mu_i <- numeric(mi*K)
  alpha_i <- alpha[[i]]
  for (j in 1:mi) {
    Y_ij0 <- get_Y_ij0(i = i,
                       j = j,
                       Y0 = Y0,
                       i_index = i_index)
    mu_ij <- get_mu_ij(Y_ij0 = Y_ij0,
                       alpha_ij = alpha_i[[j]])

    mu_i[(K*(j - 1) + 1):(K*j)] <-  mu_ij
  }
  return(mu_i)
}


#' Get the Yi minus mui term of the gradient
#'
#' @param i subject index
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param Y0
#' @param alpha
#' @param i_index
#' @param K
#'
#' @returns Vector of length Kmi
#' @export
#'
get_Yi_minus_mui <- function(i,
                             Y,
                             Y0,
                             alpha,
                             mi_vec,
                             i_index,
                             K){

  Y_i <- get_Y_i_vec(i = i,
                     mi_vec = mi_vec,
                     Y = Y)
  mu_i <- get_mu_i(i = i,
                   mi_vec = mi_vec,
                   alpha = alpha,
                   i_index = i_index,
                   Y0 = Y0,
                   K = K)
  term_i <- Y_i - mu_i
  return(term_i)
}



# B spline basis ----------------------------------------------------------
#' Get \eqn{B_ij}
#'
#' @param i subject index
#' @param j time index
#' @param i_index
#' @param B B spline basis matrix of dimension (N x P)
#'
#' @returns A vector of length P
#' @export
get_B_ij <- function(i, j, B, i_index) {
  i_start <- i_index[i] + 1
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
#' @param Z_ij vector of length l.
#' @param B_ij
#' @param i_index
#' @param P
#' @param L
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#'
#' @returns Numeric (1x1)
#' @export
get_alpha_ijk <- function(i, j, k, beta, Z_ij, B_ij, i_index, P, L) {
  # Extract all beta_lk at once for efficiency
  beta_k <- beta[((k - 1) * P + 1):(k * P), ]

  # Compute the weighted sum efficiently
  #lsum <- Z_ij * (t(B_ij) %*% beta_k)
  lsum <- Z_ij * crossprod(B_ij,beta_k)

  # Compute alpha_ijk
  alpha_ijk <- exp(sum(lsum))


  # alpha_ijk_cpp <- compute_alpha_ijk(beta,
  #                                    k = k,
  #                                    P = P,
  #                                    Z_ij = Z_ij,
  #                                    B_ij = B_ij)


  # if (any(is.infinite(alpha_ijk))) {
  #   stop("Infinite alpha")
  # }
  # if (sum(alpha_ijk < 0) != 0) {
  #   stop("Negative alpha")
  # }
  return(alpha_ijk)
}



#' Calculate alpha for all k for a single i and j.
#'
#' @param i subject index
#' @param j time index
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of responses
#' @param i_index
#' @param L
#' @param P
#'
#' @returns Vector of length K
#' @export
get_alpha_ij <- function(i, j, beta, Z, B, K, i_index, L, P) {
  B_ij <- get_B_ij(i = i,
                   j = j,
                   B = B,
                   i_index = i_index)

  # Z_ij is a vector of length L
  Z_ij <- Z[i_index[i] + j, 0:L + 1]


  alphas <- numeric(K)
  for (k in 1:K) {
    alphas[k] <- get_alpha_ijk(i = i,
                               j = j,
                               k = k,
                               beta = beta,
                               Z_ij = Z_ij,
                               B_ij = B_ij,
                               i_index = i_index,
                               L = L,
                               P = P)
  }
  return(alphas)
}


#' Title
#'
#' @param i
#' @param beta
#' @param Z
#' @param B
#' @param K
#' @param i_index
#' @param mi_vec
#' @param L
#' @param P
#'
#' @returns
#' @export
#'
#' @examples
get_alpha_i <- function(i, beta, Z, B, K, i_index, mi_vec, L, P){
  # L <- ncol(Z) - 1
  # P <- length(B_ij)
  mi <- mi_vec[i]
  alpha_i <- purrr::map(1:mi, ~get_alpha_ij(i = i,
                                           j = .x,
                                           beta = beta,
                                           Z = Z,
                                           B = B,
                                           K = K,
                                           i_index = i_index,
                                           L = L,
                                           P = P))

}

#' Title
#'
#' @param beta
#' @param Z
#' @param B
#' @param K
#' @param i_index
#' @param mi_vec
#' @param L
#' @param P
#'
#' @returns
#' @export
#'
#' @examples
get_alpha_list <- function(beta, Z, B, K, i_index, mi_vec, L, P){
  n <- length(mi_vec)
  alpha_list <- purrr::map(1:n, ~get_alpha_i(i = .x,
                                             beta = beta,
                                             Z = Z,
                                             B = B,
                                             K = K,
                                             i_index = i_index,
                                             mi_vec = mi_vec,
                                             L = L,
                                             P = P))
  return(alpha_list)

}



# Variance -----------------------------------------------------

#' Get term U_ij
#'
#' @param alpha_ij Vector of DM parameters for i, j
#'
#' @returns U ij part of partials/variance matrix
#' @export
get_U_ij <- function(alpha_ij) {
  alpha_ij0 <- sum(alpha_ij)
  U_ij <- diag(alpha_ij / alpha_ij0) - tcrossprod(alpha_ij) / alpha_ij0^2
  return(U_ij)
}

