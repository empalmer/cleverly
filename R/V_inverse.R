# Terms -----------------------------------------------------

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
#' @param i_index
#' @param Y0
#' @param L
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
                             beta = beta,
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
  if (any(is.nan(V_i))) {
    stop()
  }
  return(V_i)
}


# Inverse -----------------------------------------------------------------


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
#' @param i_index
#' @param K Number of responses
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
                       K){

  V_i_inv_list <- list()
  mi <- mi_vec[i]
  for (j in 1:mi) {
    Y_ij0 <- get_Y_ij0(i = i,
                       j = j,
                       Y0 = Y0,
                       i_index)

    alpha_ij <- alpha[[i]][[j]]

    V_ijj <- get_V_ijj(Y_ij0 = Y_ij0,
                       phi = phi,
                       alpha_ij = alpha_ij)

    V_i_inv_list[[j]] <- MASS::ginv(V_ijj)
    #V_ij_inv_list[[j]] <- V_ijj
  }

  # Invert each diagonal first!
  #V_i_inv_list <- purrr::map(V_ij_list, MASS::ginv)
  V_i_bdiag_inv <- Matrix::bdiag(V_i_inv_list)
  V_i_inv <- as.matrix(V_i_bdiag_inv)

  if (any(is.nan(V_i_inv))) {
    stop("V_i_inv has NaNs")
  }

  # Invert entire matrix.
  # Needed for non-independent case
  # start_time <- proc.time()
  # V_i_bdiag <- Matrix::bdiag(V_ij_list)
  # V_i <- as.matrix(V_i_bdiag)
  # V_i_inv_big <- MASS::ginv(V_i)
  # print(proc.time() - start_time)

  return(V_i_inv)
}

#' Get overal V inverse
#'
#' @param Y
#' @param mi_vec
#' @param i_index
#' @param phi
#' @param beta
#' @param Z
#' @param B
#' @param K
#' @param alpha
#' @param Y0
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
                      K){
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
                             K = K)
  }
  return(V_inv)
}


