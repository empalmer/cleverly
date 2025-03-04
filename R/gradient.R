# Terms -----------------------------------------------------


#' Get the Yi minus mui term of the gradient
#'
#' @param i subject index
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of responses
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns Vector of length Kmi
#' @export
#'
get_Yi_minus_mui <- function(i, Y, mi_vec, beta, Z, B, K){
  Y_i <- get_Y_i_vec(i = i, mi_vec = mi_vec, Y = Y)
  mu_i <- get_mu_i(i = i,
                   mi_vec = mi_vec,
                   Y = Y,
                   beta = beta,
                   Z = Z,
                   B = B,
                   K = K)

  term_i <- Y_i - mu_i

  return(term_i)
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
#' @param K Number of responses
#' @param V_i
#'
#' @returns Matrix of dimension Kmi x Kmi
#' @export
#'
get_Vi_inv <- function(V_i, i, Y, mi_vec, phi, beta, Z, B, K){
  if (missing(V_i)) {
    V_i <- get_V_i(i = i,
                   Y = Y,
                   phi = phi,
                   beta = beta,
                   Z = Z,
                   B = B,
                   K = K,
                   mi_vec = mi_vec)
  }
  #V_i_inv <- solve(V_i)
  if (any(is.nan(V_i))) {
    browser()
  }
  V_i_inv <- MASS::ginv(V_i)
  return(V_i_inv)
}



# Partials -----------------------------------------------------

#' Get partials matrix for i, j, l
#'
#' @param i subject index
#' @param j time index
#' @param l external variable index
#' @param K Number of responses
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns Matrix of dimension PK x K
#' @export
#'
get_partials_ijl <- function(i, j, l, mi_vec, K, Y, Z, B, beta){
  Y_ij0 <- get_Y_ij0(i = i,
                     j = j,
                     Y = Y,
                     mi_vec = mi_vec)
  U_ij <- get_U_ij(
    i = i,
    j = j,
    beta = beta,
    Z = Z,
    B = B,
    K = K,
    mi_vec = mi_vec
  )
  Z_ijl <- get_Z_ijl(
    i = i,
    j = j,
    l = l,
    Z = Z,
    mi_vec = mi_vec
  )
  B_ij <- get_B_ij(i = i,
                   j = j,
                   B = B,
                   mi_vec = mi_vec)

  partials_ijl <- Y_ij0 * Z_ijl * kronecker(U_ij, B_ij)
  return(partials_ijl)
}



#' Get full partials matrix for a given i, l
#'
#' @param i subject index
#' @param l external variable index
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#'
#' @returns Matrix of dimension KP x Kmi
#' @export
#'
get_partials_il <- function(i, l, Y, Z, B, beta, mi_vec){
  P <- ncol(B)
  K <- ncol(Y)
  mi <- mi_vec[i]
  partials_il <- matrix(nrow = K*P, ncol = K*mi)
  for (j in 1:mi) {
    partials_il[, ((j - 1)*K + 1):(j*K)] <-
      get_partials_ijl(i = i,
                       j = j,
                       l = l,
                       mi_vec = mi_vec,
                       K = K,
                       Y = Y,
                       Z = Z,
                       B = B,
                       beta = beta)
  }
  return(partials_il)
}



# Gradient ----------------------------------------------------------------


#' Title
#'
#' @param i subject index
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param phi Current value of overdispersion parameter
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param l external variable index
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param B B spline basis matrix of dimension (N x P)
#'
#' @returns Vector of length PK
#' @export
#'
get_gradient_il <- function(i, l, Y, mi_vec, phi, beta, Z, B){
  K <- ncol(Y)
  partials_il <- get_partials_il(i = i,
    l = l,
    Y = Y,
    Z = Z,
    B = B,
    beta = beta,
    mi_vec = mi_vec
  )
  V_i_inv <- get_Vi_inv(i = i,
    Y = Y,
    mi_vec = mi_vec,
    phi = phi,
    beta = beta,
    Z = Z,
    B = B,
    K = K
  )
  Yi_minus_mui <- get_Yi_minus_mui(i = i,
    Y = Y,
    mi_vec = mi_vec,
    beta = beta,
    Z = Z,
    B = B,
    K = K
  )
  gradient_il <-  partials_il %*% V_i_inv %*% Yi_minus_mui
}



#' Get the gradient for a given l
#'
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param l external variable index
#' @param phi Current value of overdispersion parameter
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#'
#' @returns vector of length PK x 1
#' @export
#'
get_gradient_l <- function(Y, mi_vec, l, phi, beta, Z, B){
  P <- ncol(B)
  K <- ncol(Y)
  gradient_sum <- numeric(P*K)
  for (i in 1:length(mi_vec)) {
    gradient_il <- get_gradient_il(
      i = i,
      l = l,
      Y = Y,
      mi_vec = mi_vec,
      phi = phi,
      beta = beta,
      Z = Z,
      B = B
    )
    gradient_sum <- gradient_sum + gradient_il
  }
  return(gradient_sum)
}


get_gradient <- function(){
  NULL
}


# Nuisance parameters ----------------------------------------------------

#' Title
#'
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param beta beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of responses
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns scalar phi
#' @export
get_phi <- function(Y, beta, Z, B, K, mi_vec){

  r <- get_pearson_residuals(Y = Y,
                             beta = beta,
                             Z = Z,
                             B = B,
                             K = K,
                             mi_vec = mi_vec)
  M <- sum(mi_vec)
  phi <- sum(r^2) / (K*M - 1)

  return(phi)
}


#' Get pearson residuals
#'
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param beta beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of responses
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns Pearson residuals vector for all i, j, k
#' @export
get_pearson_residuals <- function(Y, beta, Z, B, K, mi_vec){
  r <- c()
  for (i in 1:length(mi_vec)) {
    ri <- get_pearson_residual_i(Y, i, beta, Z, B, K, mi_vec)
    r <- c(r, ri)
  }
  return(r)
}

#' Get the pearson residual for a given i
#'
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param phi Current value of overdispersion parameter
#' @param i subject index
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of responses
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns pearson residual vector for each i
#' @export
get_pearson_residual_i <- function(Y,
                                   i,
                                   beta,
                                   Z,
                                   B,
                                   K,
                                   mi_vec){

  Yi_minus_mui <- get_Yi_minus_mui(i = i,
                                   Y = Y,
                                   mi_vec = mi_vec,
                                   beta = beta,
                                   Z = Z,
                                   B = B,
                                   K = K)
  ri <- numeric(mi_vec[i]*K)
  for (j in 1:mi_vec[i]) {
    Yij_minus_muij <- Yi_minus_mui[((j - 1)*K + 1):(K*j)]
    alpha_ij <- get_alpha_ij(beta = beta,
                             i = i,
                             j = j,
                             Z = Z,
                             B = B,
                             K = K,
                             mi_vec = mi_vec)
    Y_ij0 <- get_Y_ij0(i = i,
                       j = j,
                       Y = Y,
                       mi_vec = mi_vec)
    alpha_ij0 <- sum(alpha_ij)
    # Should be of length K.
    rij <- Yij_minus_muij /
      sqrt( Y_ij0 *
             (Y_ij0 + alpha_ij0) / (1 + alpha_ij0) *
             alpha_ij/alpha_ij0 * (1 - alpha_ij/alpha_ij0))
    ri[((j - 1)*K + 1):(K*j)] <- rij
  }

  return(ri)
}
