# Variance -----------------------------------------------------

#' Get DM variance for ith subject jth timepoint
#'
#' @param Y_ij0 Total sum of counts across all K
#' @param phi Current value of overdispersion parameter
#' @param alpha_ij Vector of DM parameters for i, j
#'
#' @returns matrix of dimension K x K for a single j
#' @export
get_V_ijj <- function(Y_ij0,
                      phi,
                      alpha_ij) {

  U_ij <- get_U_ij(alpha_ij = alpha_ij)

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
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param K Number of responses
#' @param i_index starting index of the ith subject in the data
#' @param Y0 Vector of total count for each sample
#' @param L Number of external variables
#' @param P
#'
#' @returns Diagonal matrix of Kmi x Kmi Each block is for a single j of block size K x K
#' @export
get_V_i <- function(i,
                    Y,
                    Y0,
                    phi,
                    beta,
                    Z,
                    B,
                    K,
                    mi_vec,
                    i_index,
                    L,
                    P) {

  V_ij_list <- list()
  mi <- mi_vec[i]
  for (j in 1:mi) {
    Y_ij0 <- get_Y_ij0(i = i,
                       j = j,
                       Y0 = Y0,
                       i_index)
    alpha_ij <- get_alpha_ij(i = i,
                             j = j,
                             beta_ks = beta,
                             Z = Z,
                             B = B,
                             K = K,
                             i_index = i_index,
                             L = L,
                             P = P)
    V_ijj <- get_V_ijj(Y_ij0 = Y_ij0,
                       phi = phi,
                       alpha_ij = alpha_ij)
    V_ij_list[[j]] <- V_ijj
  }
  V_i_bdiag <- Matrix::bdiag(V_ij_list)
  V_i <- as.matrix(V_i_bdiag)
  # if (any(is.nan(V_i))) {
  #   stop()
  # }
  return(V_i)
}


# Inverse Variance-----------------------------------------------------------------


#' Get Vi inverse
#'
#' @param i subject index
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param phi Current value of overdispersion parameter
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param alpha
#' @param i_index starting index of the ith subject in the data
#' @param K Number of responses
#' @param Y0
#' @param cor_str
#' @param rho_cor
#'
#' @returns Matrix of dimension Kmi x Kmi
#' @export
#'
get_Vi_inv <- function(i,
                       Y,
                       Y0,
                       alpha,
                       mi_vec,
                       i_index,
                       phi,
                       beta,
                       Z,
                       B,
                       K,
                       cor_str,
                       rho_cor){



  if (cor_str == "IND") {
    # V_i_inv_list <- list()
    mi <- mi_vec[i]
    V_i_inv_mat <- matrix(0, nrow = mi*K, ncol = mi*K)
    for (j in 1:mi) {
      Y_ij0 <- get_Y_ij0(i = i,
                         j = j,
                         Y0 = Y0,
                         i_index)

      alpha_ij <- alpha[[i]][[j]]

      V_ijj <- get_V_ijj(Y_ij0 = Y_ij0,
                         phi = phi,
                         alpha_ij = alpha_ij)
      # If independent, can invert each block individually
      V_i_inv_mat[((j - 1)*K + 1):(j*K), ((j - 1)*K + 1):(j*K)] <- MASS::ginv(V_ijj)
    }
    V_i_inv <- V_i_inv_mat
  }
  else if (cor_str == "CON") {
    V_i_list <- list()
    R_i_list <- list()
    mi <- mi_vec[i]
    V_i <- matrix(nrow = 0, ncol = 0)
    R_i <- matrix(nrow = 0, ncol = 0)
    for (j in 1:mi) {
      Y_ij0 <- get_Y_ij0(i = i,
                         j = j,
                         Y0 = Y0,
                         i_index)
      alpha_ij <- alpha[[i]][[j]]


      V_ijj <- get_V_ijj(Y_ij0 = Y_ij0,
                         phi = phi,
                         alpha_ij = alpha_ij)
      A_ijj <- diag(1/sqrt(diag(V_ijj)))

      # Call C++ version of multiplying three matrices
      #R_ijj <- A_ijj %*% V_ijj %*% A_ijj
      R_ijj <- fast_mat_mult3(A_ijj, V_ijj, A_ijj)

      V_i <- magic::adiag(V_i, V_ijj)
      R_i <- magic::adiag(R_i, R_ijj)

    }

    A_inv <- diag(1/sqrt(diag(V_i)))

    # For CON, this is a matrix of 1s everywhere except on the block diagonal which has 0s
    # Then scaled by the calculated rho_cor
    corR <- rho_cor * kronecker(matrix(1, nrow = mi, ncol = mi) - diag(mi),
                                matrix(1, nrow = K, ncol = K))

    #V_i_inv <- (1/phi) * A_inv %*% MASS::ginv(R_i + corR) %*% A_inv
    V_i_inv <- (1/phi) * fast_mat_mult3(A_inv, MASS::ginv(R_i + corR), A_inv)
  }
  else if (cor_str == "AR1") {
    V_i_list <- list()
    R_i_list <- list()
    mi <- mi_vec[i]
    V_i <- matrix(nrow = 0, ncol = 0)
    R_i <- matrix(nrow = 0, ncol = 0)
    for (j in 1:mi) {
      Y_ij0 <- get_Y_ij0(i = i,
                         j = j,
                         Y0 = Y0,
                         i_index)
      alpha_ij <- alpha[[i]][[j]]

      V_ijj <- get_V_ijj(Y_ij0 = Y_ij0,
                         phi = phi,
                         alpha_ij = alpha_ij)
      A_ijj <- diag(1/sqrt(diag(V_ijj)))

      #R_ijj <- A_ijj %*% V_ijj %*% A_ijj
      R_ijj <- fast_mat_mult3(A_ijj, V_ijj, A_ijj)

      V_i <- magic::adiag(V_i, V_ijj)
      R_i <- magic::adiag(R_i, R_ijj)

    }
    A_inv <- diag(1/sqrt(diag(V_i)))

    # For AR1, this is rho^{|j1 - j2|}
    j1s <- matrix(rep(1:mi, each = mi), nrow = mi)
    j2s <- matrix(rep(1:mi, each = mi), ncol = mi, byrow = TRUE)
    corR_ar1 <- kronecker(rho_cor^abs(j1s - j2s), matrix(1, nrow = K, ncol = K))

    # Use C++ version of multiplying three matrices
    #V_i_inv <- (1/phi) * A_inv %*% MASS::ginv(R_i + corR) %*% A_inv
    V_i_inv <- (1/phi) * fast_mat_mult3(A_inv, MASS::ginv(R_i + corR_ar1), A_inv)

  }
  else{
    stop("Invalid corstr")
  }

  return(V_i_inv)
}

#' Get overall V inverse
#'
#' for all i
#'
#' @param Y
#' @param mi_vec
#' @param i_index starting index of the ith subject in the data
#' @param phi Current value of overdispersion parameter
#' @param beta
#' @param Z
#' @param B
#' @param K Number of responses
#' @param alpha
#' @param Y0 Vector of total count for each sample
#'
#' @returns
#' @export
#'
#' @examples
get_V_inv <- function(Y,
                      Y0,
                      alpha,
                      mi_vec,
                      i_index,
                      phi,
                      beta,
                      Z,
                      B,
                      K,
                      cor_str,
                      rho_cor){

  V_inv <- list()
  for (i in 1:length(mi_vec)) {
    V_inv[[i]] <- get_Vi_inv(i = i,
                             Y = Y,
                             Y0 = Y0,
                             alpha = alpha,
                             mi_vec = mi_vec,
                             i_index = i_index,
                             phi = phi,
                             beta = beta,
                             Z = Z,
                             B = B,
                             K = K,
                             cor_str = cor_str,
                             rho_cor = rho_cor)
  }
  return(V_inv)
}



# Correlation structures --------------------------------------------------


#' Title
#'
#' @param corstr
#' @param K Number of responses
#'
#' @returns
#' @export
#'
#' @examples
get_R_ijjp <- function(corstr, rho, K){
  if (corstr == "IND") {
    str <- diag(K) * corstr
    R_ijjp <- str * rho
  }
  else if (corstr == "AR1") {
    R_ijjp <- outer(1:K, 1:K, function(i, j) rho^abs(i - j))
  }
  else if (corstr == "CON") {
    str <- matrix(1, nrow = K, ncol = K) - diag(K)
    R_ijjp <- str * rho
  }
  else{
    stop("Invalid corstr")
  }
  return(R_ijjp)
}



#' Get rho
#'
#' @returns
#' @export
#'
#' @examples

get_rho <- function(Y,
                    Y0,
                    beta,
                    alpha,
                    Z,
                    B,
                    K,
                    mi_vec,
                    i_index,
                    M,
                    cor_str){

  if (cor_str == "IND") {
    return(0)
  }


  # Initialize list for each subject
  regression_data_list <- vector("list", length(mi_vec))

  for (i in seq_along(mi_vec)) {
    mi <- mi_vec[i]
    capture_number <- 1:mi_vec[i]

    # Compute the difference matrix efficiently
    diag_values <- rep(capture_number, each = K)
    # This lets us find the blocks of the matrix, only matters
    # if they are zero, not the actual value.
    block_diagonal_elements <- outer(diag_values, diag_values, "-")

    # Compute Pearson residual matrix
    pearson_residual_i <- get_pearson_residual_i(Y,
                                                 Y0,
                                                 i,
                                                 beta,
                                                 alpha,
                                                 Z,
                                                 B,
                                                 K,
                                                 mi_vec,
                                                 i_index)
    matrix_pearson_residual <- tcrossprod(pearson_residual_i)

    # We only want the upper non-block triangle of the matrix
    # Extract non-NA values (which are the upper triangular non-block)
    rijk_rijk <- matrix_pearson_residual[!(block_diagonal_elements == 0 |
                                            lower.tri(matrix_pearson_residual))]

    # Used for AR1 cor str
    abs_j1_j2 <- abs(block_diagonal_elements[!(block_diagonal_elements == 0 |
                                            lower.tri(matrix_pearson_residual))])

    # save for i
    regression_data_list[[i]] <- data.frame(rijk_rijk,
                                            abs_j1_j2)

  }
  # Combine all dataframes at the end for efficiency
  regressiondata <- dplyr::bind_rows(regression_data_list)


  if (cor_str == "CON") {
    # This is how we calculate it for the CON structure.
    rho_cor <- mean(regressiondata$rijk_rijk)
  }
  if (cor_str == "AR1") {
    objective_f <- purrr::map_dbl(seq(-1,1,0.1), ~sum((regressiondata$rijk_rijk -
                                                         .x^(regressiondata$abs_j1_j2))^2) )
    rho_cor <- seq(-1,1,0.1)[which.min(objective_f)]
  }

  # logic for edge cases, mean might not give a valid structure.
  if (rho_cor == 1) {
    rho_cor <- -1
  }
  if (rho_cor > 1) {
    rho_cor <- 1
  }
  return(rho_cor)

}



