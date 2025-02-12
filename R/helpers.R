# Various -----------------------------------------------------------------
#' Get \eqn{Y_{ij0}}
#'
#' @param i subject
#' @param j time index
#' @param Y response
#' @param mi number of time points for subject i
#'
#' @returns Scalar of the total sum constraint for a given i, j
#'
get_Yij0 <- function(i, j, Y, mi) {
  Yij0 <- sum(Y[((i - 1) * mi + j), ])
  return(Yij0)
}



#' Get \eqn{B_ij}
#'
#' @param i subject
#' @param j time index
#' @param B B-spline basis matrix
#' @param mi number of time points for subject i
#'
#' @returns
#'
get_Bij <- function(i, j, B, mi) {
  Bij <- B[((i - 1) * mi + j), ]
  return(Bij)
}





#' Get \eqn{Z_{ijl}}
#'
#' @param i subject
#' @param j time index
#' @param l external variable index
#' @param Z matrix of external variables
#'
#' @returns
#'
#' @examples
get_Z_ijl <- function(i, j, l, Z) {
  Z_ijl <- Z[i + (j - 1), (l + 1)]
  return(Z_ijl)
}


# Alphas ------------------------------------------------------------------

#' Calculate alpha for a single \eqn{i = 1, j = 1, k = 1} from \eqn{\beta}
#'
#' @param i
#' @param j
#' @param k
#' @param beta
#' @param Z
#' @param B
#'
#' @returns
#'
#' @examples
get_alpha_ijk <- function(i, j, k, beta, Z, B){
  lsum <- numeric(L)
  for (l in 0:L) {
    Z_ijl <- Z[i + (j - 1), (l + 1)]
    beta_lk <- beta[((k - 1)*P + 1):(k*P), (l + 1)]
    B_ij <- B[((j - 1)*mi + i), ]

    lsum <- Z_ijl * t(B_ij) %*% beta_lk
  }
  return(exp(sum(lsum)))
}


#' Calculate alpha for all k for a single i and j.
#'
#' @param i
#' @param j
#' @param beta
#' @param Z
#' @param B
#'
#' @returns
#' @export
#'
#' @examples
get_alpha_ij <- function(i, j, beta, Z, B){
  alphas <- numeric(K)
  for (k in 1:K) {
    alphas[k] <- get_alpha_ijk(i, j, k, beta, Z, B)
  }
  return(alphas)
}




# Mus ------------------------------------------------------------------

get_mu_ij <- function(Yij0, alpha_ij, ... ){
  if (missing(alpha_ij)) {
    # calculate with i, j, beta, Z, B
    alpha_ij <- get_alpha_ij(...)
  }
  mu_ij <- Yij0/sum(alpha_ij) * alpha_ij
  return(mu_ij)
}

get_mu_ij(Yij0 = 100, alpha_ij = get_alpha_ij(i = 1, j = 2, beta, Z, B))
get_mu_ij(Yij0 = 100, i = 1, j = 2, beta = beta, Z = Z, B = B)


get_mu_i <- function(i, mi, Y, beta, Z, B){
  for (j in 1:mi) {
    Y_ij0 <- get_Yij0(i, j, Y)
    mu_ij <- get_mu_ij(Y_ij0, i, j, beta, Z, B)
  }
}


# Variance -----------------------------------------------------

get_U_ij <- function(Y_ij0, alpha_ij,...){
  if (missing(alpha_ij)) {
    # calculate with i, j, beta, Z, B
    alpha_ij <- get_alpha_ij(...)
  }
  alpha_ij0 <- sum(alpha_ij)

  U_ij <- diag(alpha_ij/alpha_ij0) - alpha_ij %*% t(alpha_ij)/alpha_ij0^2
  return(U_ij)
}

get_U_ij(Yij0 = 100, alpha_ij = get_alpha_ij(i = 1, j = 2, beta, Z, B))


get_V_ijj <- function(Y_ij0, phi, alpha_ij, ...){
  if (missing(alpha_ij)) {
    # calculate with i, j, beta, Z, B
    alpha_ij <- get_alpha_ij(...)
  }
  U_ij <- get_U_ij(Y_ij0, alpha_ij)

  V_ijj <- phi * Y_ij0 * U_ij
  return(V_ijj)
}

get_V_i <- function(i)


  # Partials -----------------------------------------------------

get_partials_ijl <- function(i, j, l, Y_ij0, Z, B, beta){
  # U_ij will calculate alphas based on i, j, Z, beta.
  U_ij <- get_U_ij(i = i, j = j, beta = beta, Z = Z, B = B)
  Z_ijl <- get_Z_ijl(i, j, l, Z)
  B_ij <- get_Bij(i, j, B)

  partials_ijl <- Y_ij0 * Z_ijl * kronecker(U_ij, B_ij)

  return(partials_ijl)
}

get_partials_ijl(i = 1, j = 2, l = 1, Y_ij0 = 100, Z, B, beta)


get_partials_il <- function(i, l, Y, Z, B, beta, mi = 3, K, P){
  partials_il <- matrix(nrow = K*P, ncol = K*mi)
  for (j in 1:mi) {
    partials_il[, ((j - 1)*K + 1):(j*K)] <-
      get_partials_ijl(i, j, l, Y_ij0 = get_Yij0(i, j, Y), Z, B, beta)
  }
  return(partials_il)
}

get_partials_il(i = 1, l = 1, Y, Z, B, beta, mi = 3, K, P) %>% View()




# Gradient -----------------------------------------------------


get_Yi_minus_mui <- function(i, Y){
  NULL
}


get_gradient_i <- function(){
  partials_i <- get_partials_i(){

  }
}

get_gradient <- function(){
  NULL
}





# Hessian -----------------------------------------------------------------

get_Hessian <- function(){
  NULL
}


