# Gradient ----------------------------------------------------------------


#' Get gradient for a given i, l
#'
#' @param i subject index
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param phi Current value of overdispersion parameter
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param l external variable index
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param B B spline basis matrix of dimension (N x P)
#' @param i_index starting index of the ith subject in the data
#' @param Y0
#' @param alpha
#' @param Vi_inv
#' @param partials_il
#' @param K
#'
#' @returns Vector of length PK
#' @export
#'
get_gradient_il <- function(i,
                            l,
                            Y,
                            Y0,
                            mi_vec,
                            i_index,
                            phi,
                            beta,
                            alpha,
                            Z,
                            B,
                            Vi_inv,
                            partials_il,
                            K){

  Yi_minus_mui <- get_Yi_minus_mui(i = i,
                                   Y = Y,
                                   Y0 = Y0,
                                   alpha = alpha,
                                   mi_vec = mi_vec,
                                   i_index = i_index,
                                   K = K)

  #gradient_il <-  partials_il %*% Vi_inv %*% Yi_minus_mui

  #term1 <- fast_mat_mult2(partials_il, Vi_inv)
  #gradient_il <- fast_mat_mult2(term1, matrix(Yi_minus_mui))


  gradient_il <- fast_mat_mult3(partials_il, Vi_inv, matrix(Yi_minus_mui))

  return(gradient_il)
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
#' @param i_index
#' @param V_inv
#' @param partials_l
#' @param Y0
#' @param alpha
#' @param K
#' @param P
#'
#' @returns vector of length PK x 1
#' @export
#'
get_gradient_l <- function(Y,
                           Y0,
                           mi_vec,
                           i_index,
                           l,
                           phi,
                           beta,
                           alpha,
                           Z,
                           B,
                           V_inv,
                           partials_l,
                           K,
                           P){

  gradient_sum <- numeric(P*K)
  for (i in 1:length(mi_vec)) {
    gradient_il <- get_gradient_il(i = i,
                                   l = l,
                                   Y = Y,
                                   Y0 = Y0,
                                   mi_vec = mi_vec,
                                   i_index = i_index,
                                   phi = phi,
                                   beta = beta,
                                   alpha = alpha,
                                   Z = Z,
                                   B = B,
                                   Vi_inv = V_inv[[i]],
                                   partials_il = partials_l[[i]],
                                   K = K)
    gradient_sum <- gradient_sum + gradient_il
  }
  return(gradient_sum)
}





