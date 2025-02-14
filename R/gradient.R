# Gradient -----------------------------------------------------


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
get_Yi_minus_mui <- function(i, Y, mi, beta, Z, B, K){
  Y_i <- get_Y_i_vec(i, Y, mi)
  mu_i <- get_mu_i(i, mi, Y, beta, Z, B, K)

  term_i <- Y_i - mu_i

  return(term_i)
}


get_gradient_i <- function(){
  NULL
}

get_gradient <- function(){
  NULL
}

# Partials -----------------------------------------------------

get_partials_ijl <- function(i, j, l, Y_ij0, Z, B, beta){
  # U_ij will calculate alphas based on i, j, Z, beta.
  U_ij <- get_U_ij(i = i, j = j, beta = beta, Z = Z, B = B)
  Z_ijl <- get_Z_ijl(i, j, l, Z)
  B_ij <- get_B_ij(i, j, B)

  partials_ijl <- Y_ij0 * Z_ijl * kronecker(U_ij, B_ij)

  return(partials_ijl)
}

#get_partials_ijl(i = 1, j = 2, l = 1, Y_ij0 = 100, Z, B, beta)


get_partials_il <- function(i, l, Y, Z, B, beta, mi = 3, K, P){
  partials_il <- matrix(nrow = K*P, ncol = K*mi)
  for (j in 1:mi) {
    partials_il[, ((j - 1)*K + 1):(j*K)] <-
      get_partials_ijl(i, j, l, Y_ij0 = get_Y_ij0(i, j, Y), Z, B, beta)
  }
  return(partials_il)
}

#get_partials_il(i = 1, l = 1, Y, Z, B, beta, mi = 3, K, P) %>% View()

