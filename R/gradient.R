# Terms -----------------------------------------------------


#' Get the Yi minus mui term of the gradient
#'
#' @param i
#' @param Y
#' @param mi
#' @param beta
#' @param Z
#' @param B
#' @param K
#'
#' @returns
#' @export
#'
#' @examples
get_Yi_minus_mui <- function(i, Y, mis, beta, Z, B, K){
  Y_i <- get_Y_i_vec(i = 1, mis = mis, Y = Y)
  mu_i <- get_mu_i(i = 1,
                   mis = mis,
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
#' @param i
#' @param Y
#' @param mi
#' @param phi
#' @param beta
#' @param Z
#' @param B
#' @param K
#'
#' @returns
#' @export
#'
#' @examples
get_Vi_inv <- function(V_i, i, Y, mis, phi, beta, Z, B, K){
  if (missing(V_i)) {
    V_i <- get_V_i(i = i,
                   Y = Y,
                   phi = phi,
                   beta = beta,
                   Z = Z,
                   B = B,
                   K = K,
                   mis = mis)
  }
  #V_i_inv <- solve(V_i)
  V_i_inv <- MASS::ginv(V_i)
  return(V_i_inv)
}



# Partials -----------------------------------------------------

#' Get partials matrix for i, j, l
#'
#' @param i
#' @param j
#' @param l
#' @param mi
#' @param K
#' @param Y
#' @param Z
#' @param B
#' @param beta
#'
#' @returns Matrix of dimension PK x K
#' @export
#'
#' @examples
get_partials_ijl <- function(i, j, l, mis, K, Y, Z, B, beta){
  Y_ij0 <- get_Y_ij0(i = i,
                     j = j,
                     Y = Y,
                     mis = mis)
  U_ij <- get_U_ij(
    i = i,
    j = j,
    beta = beta,
    Z = Z,
    B = B,
    K = K,
    mis = mis
  )
  Z_ijl <- get_Z_ijl(
    i = i,
    j = j,
    l = l,
    Z = Z,
    mis = mis
  )
  B_ij <- get_B_ij(i = i,
                   j = j,
                   B = B,
                   mis = mis)

  partials_ijl <- Y_ij0 * Z_ijl * kronecker(U_ij, B_ij)
  return(partials_ijl)
}



#' Get full partials matrix for a given i, l
#'
#' @param i
#' @param l
#' @param Y
#' @param Z
#' @param B
#' @param beta
#' @param mi
#' @param K
#' @param P
#'
#' @returns
#' @export
#'
#' @examples
get_partials_il <- function(i, l, Y, Z, B, beta, mis){
  P <- ncol(B)
  K <- ncol(Y)
  mi <- mis[i]
  partials_il <- matrix(nrow = K*P, ncol = K*mi)
  for (j in 1:mi) {
    partials_il[, ((j - 1)*K + 1):(j*K)] <-
      get_partials_ijl(i = i,
                       j = j,
                       l = l,
                       mis = mis,
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
#' @param i
#' @param Y
#' @param mi
#' @param phi
#' @param beta
#' @param Z
#' @param B
#' @param K
#'
#' @returns
#' @export
#'
#' @examples
get_gradient_il <- function(i, l, Y, mis, phi, beta, Z, B){
  K <- ncol(Y)
  partials_il <- get_partials_il(i = i,
    l = l,
    Y = Y,
    Z = Z,
    B = B,
    beta = beta,
    mis = mis
  )
  V_i_inv <- get_Vi_inv(i = i,
    Y = Y,
    mis = mis,
    phi = phi,
    beta = beta,
    Z = Z,
    B = B,
    K = K
  )
  Yi_minus_mui <- get_Yi_minus_mui(i = i,
    Y = Y,
    mis = mis,
    beta = beta,
    Z = Z,
    B = B,
    K = K
  )
  gradient_il <-  partials_il %*% V_i_inv %*% Yi_minus_mui
}



#' Title
#'
#' @param Y Matrix of dimension M x K
#' @param subject_ids numeric
#' @param time_ids numeric
#' @param l external variable index
#' @param phi variance parameter
#' @param beta
#' @param Z
#' @param B
#' @param K
#'
#' @returns vector of length PK x 1
#' @export
#'
#' @examples
get_gradient_l <- function(Y, subject_ids, time_ids, l, phi, beta, Z, B){
  P <- ncol(B)
  K <- ncol(Y)
  gradient_sum <- numeric(P*K)
  Y_values <- get_Y_wrapper(Y = Y, subject_ids = subject_ids, time_ids = time_ids)
  Y_use <- Y_values$Y
  subject_ids_use <- Y_values$subject_id_values


  mis <- get_mis(Y, subject_ids, time_ids)$mi

  for (i in 1:length(mis)) {
    gradient_il <- get_gradient_il(
      i = i,
      l = l,
      Y = Y_use,
      mis = mis,
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
