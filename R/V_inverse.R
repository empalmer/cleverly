# Variance -----------------------------------------------------

#' Get term U_ij
#'
#' @param alpha_ij Vector of DM parameters for i, j
#'
#' @returns U ij part of partials/variance matrix
#' @export
get_U_ij <- function(alpha_ij) {
  alpha_ij0 <- sum(alpha_ij)
  if (any(alpha_ij == alpha_ij0)) {
    alpha_ij0 <- alpha_ij0 + 1e-10
    warning("alpha_ij = alphaij0, adding 1e-10")
    print("alpha_ij = alphaij0, adding 1e-10")
  }
  U_ij <- diag(alpha_ij / alpha_ij0) - tcrossprod(alpha_ij) / alpha_ij0^2
  return(U_ij)
}


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

  alpha_ij0 <- sum(alpha_ij)
  DM_overdispersion_ij <- (Y_ij0 + alpha_ij0)/(1 + alpha_ij0)
  #V_ijj <- phi * Y_ij0 * U_ij
  V_ijj <- Y_ij0 * U_ij * DM_overdispersion_ij



  return(V_ijj)
}



# Inverse Variance-----------------------------------------------------------------

#' Get overall V inverse
#'
#' for all i
#'
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param i_index starting index of the ith subject in the data
#' @param phi Current value of overdispersion parameter
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of responses
#' @param alpha list of alpha that can be subsetted by i and j
#' @param Y0 Vector of total count for each sample
#' @param cor_str correlation structures used ie. IND, CON, etc
#' @param rho_cor value of rho (correlation parameter)
#'
#' @returns V inverse list
#' @export
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


#' Get Vi inverse
#'
#' @param i subject index
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param phi Current value of overdispersion parameter
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param alpha list of alpha that can be subsetted by i and j
#' @param i_index starting index of the ith subject in the data
#' @param K Number of responses
#' @param Y0 Vector of total count for each sample
#' @param cor_str Type of correlation structure (IND, CON, AR1)
#' @param rho_cor Current value of rho
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


      # cond_num <- kappa(V_ijj, exact = F)
      # print(cond_num)

      # If independent, can invert each block individually
      # And each block is just the variance components, so we can
      #V_i_inv_mat[((j - 1)*K + 1):(j*K), ((j - 1)*K + 1):(j*K)] <- MASS::ginv(V_ijj)
      V_i_inv_mat[((j - 1)*K + 1):(j*K), ((j - 1)*K + 1):(j*K)] <- (1/phi) * MASS::ginv(V_ijj)

    }
    V_i_inv <- V_i_inv_mat
  }
  else if (cor_str == "CON" | cor_str == "AR1" | cor_str == "CON-d" | cor_str == "AR1-d") {
    # If not independent, we have to invert the entire matrix, and there are
    # cross correlations
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
      # R ijj is the independent part of Rijj, in the next code block we will add
      # corR to R, which will account for the correlation structure.
      R_ijj <- fast_mat_mult3(A_ijj, V_ijj, A_ijj)

      V_i <- magic::adiag(V_i, V_ijj)
      R_i <- magic::adiag(R_i, R_ijj)

    }
    corR <- get_corR(cor_str = cor_str,
                     mi = mi,
                     K = K,
                     rho = rho_cor)

    # Calculate the inverse after calculating the all terms:
    A_inv <- diag(1/sqrt(diag(V_i)))



    # R_inv <- MASS::ginv(Ri)
    R_sum <- R_i + corR
    R_inv <- MASS::ginv(R_sum)



    # R_sum_pd <- as.matrix(Matrix::nearPD(R_sum)$mat)
    # R_inv <- MASS::ginv(R_sum_pd)

    # cond_num <- kappa(R_sum, exact = F)
    # print(cond_num)
    #
    # epsilon <- 1e-2
    # R_sum_e <- R_i + corR + diag(epsilon, nrow(R_i))
    # cond_num <- kappa(R_sum, exact = F)
    # print(paste0("with epsilon: ", cond_num))
    #
    #R_inv <- MASS::ginv(R_sum)

    # V_i_inv <- (1/phi) * A_inv %*% MASS::ginv(R_i + corR) %*% A_inv
    V_i_inv <- (1/phi) * fast_mat_mult3(A_inv, R_inv , A_inv)
  }

  else{
    stop("Invalid cor_str")
  }

  return(V_i_inv)
}




# Correlation structures --------------------------------------------------

get_corR <- function(cor_str, mi, K, rho) {
  if (cor_str == "IND") {
    return(list())
  }
  # Calculate the longitudinal blocks of the correlation structure
  else if (cor_str == "CON") {
    # For CON, this is a matrix of 1s everywhere except on the block diagonal which has 0s
    # Then scaled by the calculated rho_cor
    # One block correlation matrix
    cor_str_matrix <- matrix(1, nrow = mi, ncol = mi) - diag(mi)
    # Make the full R matrix, which is blocks of the cor_str matrix (the same for all blocks)
    all_blocks <- matrix(1, nrow = K, ncol = K)
    block_cor_matrix <- kronecker(cor_str_matrix, all_blocks)
    corR <- rho * block_cor_matrix
  } else if (cor_str == "AR1") {
    # For AR1, this is rho^{|j1 - j2|}
    j1s <- matrix(rep(1:mi, each = mi), nrow = mi)
    j2s <- matrix(rep(1:mi, each = mi), ncol = mi, byrow = TRUE)
    cor_str_matrix <- rho^abs(j1s - j2s) - diag(mi)
    all_blocks <- matrix(1, nrow = K, ncol = K)
    corR <- kronecker(cor_str_matrix, all_blocks)
  } else if (cor_str == "CON-d") {
    # One block correlation matrix
    cor_str_matrix <- matrix(1, nrow = mi, ncol = mi) - diag(mi)
    diag_blocks <- diag(K)
    block_cor_matrix <- kronecker(cor_str_matrix, diag_blocks)
    corR <- rho * block_cor_matrix
  } else if (cor_str == "AR1-d") {
    j1s <- matrix(rep(1:mi, each = mi), nrow = mi)
    j2s <- matrix(rep(1:mi, each = mi), ncol = mi, byrow = TRUE)

    cor_str_matrix <- rho^abs(j1s - j2s) - diag(mi)
    diag_blocks <- diag(K)
    corR <- kronecker(cor_str_matrix, diag_blocks)
  }
  return(corR)
}

#' Get rho
#'
#' @param K Number of responses
#' @param mi_vec vector of the number of time points for each sample. Of length n
#' @param M Number of subjects times time points for each subject
#' @param cor_str Type of correlation structure (IND, CON, AR1)
#' @param pearson_residuals list of pearson residuals
#' @param phi value of phi
#' @param cor_blocks Pre-calculated to determine which pearson resids to use to estimate rho
#' @param j1_j2_list Used for Ar1 calculations j1-j2
#'
#' @returns Numeric estimate of rho
#' @export
get_rho <- function(pearson_residuals,
                    phi,
                    K,
                    mi_vec,
                    M,
                    cor_str,
                    cor_blocks,
                    j1_j2_list){

  if (cor_str == "IND") {
    return(0)
  }


  # Initialize list for each subject
  regression_data_list <- vector("list", length(mi_vec))

  # Calculate for each subject i
  for (i in seq_along(mi_vec)) {

    # Extract ith pearson residual vector
    pearson_residual_i <- pearson_residuals[[i]]
    # Create a matrix of rijk * rijk' residuals
    matrix_pearson_residual <- tcrossprod(pearson_residual_i)

    # Calculate corR just to select the right elements of the matrix to include in calculations
    # The non-zero elements of corR are what we want to select from the pearson resid matrix.
    # Here we set rho = 1 just so we can get the non-zero elements, really it can be any nonzero number
    # cor_blocks <- get_corR(cor_str = cor_str,
    #                  mi = mi_vec[i],
    #                  K = K,
    #                  rho = 1)

    # We only want the upper non-block triangle of the matrix
    # Extract non-NA values (which are the upper triangular non-block)
    #rijk_rijk <- matrix_pearson_residual[is_off_bdiag]
    rijk_rijk <- matrix_pearson_residual[cor_blocks[[i]] != 0]
    # Normalize by phi
    normalized_rijk_rijk <- 1/phi * rijk_rijk


    if (cor_str == "AR1" | cor_str == "AR1-d") {
      # For AR1, we need to get the absolute value of the j1 - j2
      # This is the same for both AR1 and AR1-d
      # Used for AR1 and AR1-d cor structure
      j1_j2 <- j1_j2_list[[i]]
      abs_j1_j2 <- abs(j1_j2[cor_blocks[[i]] != 0])

      # save for i
      regression_data_list[[i]] <- data.frame(normalized_rijk_rijk,
                                              abs_j1_j2)
    } else {
      # save for i
      regression_data_list[[i]] <- data.frame(normalized_rijk_rijk)
    }

  }
  # Combine all i data frames at the end
  regressiondata <- dplyr::bind_rows(regression_data_list)


  if (cor_str == "CON") {
    # This is how we calculate it for the CON structure.
    # The chosen rijk_rijk are different depending on if it is con or cond
    rho_cor <- mean(regressiondata$normalized_rijk_rijk)

  }
  if (cor_str == "AR1" | cor_str == "AR1-d") {
    # Find the minimum
    # Using seq() will speed up the process
    # The chosen rijk_rijk are different depending on if it is ar1 or ar1d
    objective_f <- purrr::map_dbl(seq(-1,1,0.1), ~sum((regressiondata$normalized_rijk_rijk -
                                                         .x^(regressiondata$abs_j1_j2))^2) )
    rho_cor <- seq(-1,1,0.1)[which.min(objective_f)]

  }
  if (cor_str == "CON-d") {
    # objective_f <- purrr::map_dbl(seq(-1,1,0.1),
    #                               ~sum((regressiondata$normalized_rijk_rijk - .x)^2) )
    # rho_cor <- seq(-1,1,0.1)[which.min(objective_f)]
    rho_cor <- mean(regressiondata$normalized_rijk_rijk)

  }




  # logic for edge cases, mean might not give a valid structure.
  if (rho_cor < (-1)) {
    rho_cor <- -1
  }
  if (rho_cor > 1) {
    rho_cor <- 1
  }

  return(rho_cor)

}



