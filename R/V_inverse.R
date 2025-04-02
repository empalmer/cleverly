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
  # if (any(is.nan(V_i))) {
  #   stop()
  # }
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
                       K,
                       corstr,
                       rho_cor){



  if (corstr == "IND") {
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
    #V_i_inv <- Matrix::bdiag(V_i_inv_list)
    V_i_bdiag_inv <- Matrix::bdiag(V_i_inv_list)
    V_i_inv <- as.matrix(V_i_bdiag_inv)

    # if (any(is.nan(V_i_inv))) {
    #   stop("V_i_inv has NaNs")
    # }

    # Invert entire matrix.
    # Needed for non-independent case
    # start_time <- proc.time()
    # V_i_bdiag <- Matrix::bdiag(V_ij_list)
    # V_i <- as.matrix(V_i_bdiag)
    # V_i_inv_big <- MASS::ginv(V_i)
    # print(proc.time() - start_time)
  }
  else if (corstr == "CON") {
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
      R_ijj <- diag(1/sqrt(diag(V_ijj))) %*% V_ijj %*% diag(1/sqrt(diag(V_ijj)))


      V_i <- magic::adiag(V_i, V_ijj)
      R_i <- magic::adiag(R_i, R_ijj)

    }

    A_inv <- diag(1/sqrt(diag(V_i)))

    # For CON, this is a matrix of 1s everywhere except on the block diagonal which has 0s
    # Then scaled by the calculated rho_cor
    corR <- rho_cor * kronecker(matrix(1, nrow = mi, ncol = mi) - diag(mi),
                                matrix(1, nrow = K, ncol = K))

    V_i_inv <- (1/phi) * A_inv %*% MASS::ginv(R_i + corR) %*% A_inv

  }
  else{
    stop("Invalid corstr")
  }


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
                      K,
                      corstr,
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
                             corstr = corstr,
                             rho_cor = rho_cor)
  }
  return(V_inv)
}



# Correlation structures --------------------------------------------------


#' Title
#'
#' @param corstr
#' @param K
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



#' Title
#'
#' @returns
#' @export
#'
#' @examples
get_rho_con <- function(Y,
                        Y0,
                        beta,
                        alpha,
                        Z,
                        B,
                        K,
                        mi_vec,
                        i_index,
                        M){


  regressiondata <- NULL
  #for (subject in unique(individual_time$individual)){
  for (i in 1:length(mi_vec)) {

    mi <- mi_vec[i]
    capture_number <- 1:mi_vec[i]


    j1minusj2 <- matrix(1, nrow = mi*K, ncol = mi*K)


    ####diagonal matrix#####
    diagnal2 <- kronecker(matrix(capture_number,nrow = 1),
                          matrix(rep(1,K),nrow = 1))
    diagnal3 <- matrix(rep(diagnal2,each=length(diagnal2)),
                       nrow=length(diagnal2)) -
      t(matrix(rep(diagnal2,each=length(diagnal2)),
               nrow=length(diagnal2)))


    ###response######
    pearson_residual_i <- get_pearson_residual_i(Y,
                                                 Y0,
                                                 i,
                                                 beta,
                                                 alpha = alpha,
                                                 Z,
                                                 B,
                                                 K,
                                                 mi_vec,
                                                 i_index)
    pearson_residual <- pearson_residual_i

    matrix_pearson_residual <- matrix(pearson_residual,ncol = 1) %*%
      matrix(pearson_residual,nrow = 1)

    pearson_resid_sq <- pearson_residual %*% t(pearson_residual)


    matrix_pearson_residual[round(diagnal3) == 0] = NA
    matrix_pearson_residual[lower.tri(matrix_pearson_residual)] = NA


    Rijkl <- as.vector( matrix_pearson_residual)
    newRijkl <- Rijkl[is.na(Rijkl) == FALSE]

    abs_j1_j2 <- round(abs(as.vector( j1minusj2))[is.na(Rijkl) == FALSE],0)
    fordata <- data.frame(newRijkl ,abs_j1_j2)
    regressiondata <- rbind(regressiondata, fordata)
  }

  # This is how we calculate it for the CON structure.
  rho_cor <- mean(regressiondata$newRijkl)
  if (rho_cor == 1) {
    rho_cor <- -1
  }
  if (rho_cor > 1) {
    rho_cor <- 1
  }
  return(rho_cor)

}



