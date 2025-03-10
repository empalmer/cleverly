#' Algorithm 2: Updating the non clustering betas
#'
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#' @param lp clustering index (integer between 0 and L)
#' @param B B spline basis matrix of dimension (N x P)
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param D Matrix of dth order weighted difference operator
#' @param gammas Vector of dimension L + 1 for penalizing the D matrix
#' @param phi Current value of overdispersion parameter
#' @param C Constant for determining the hessian change.
#' @param i_index
#' @param V_inv
#'
#' @returns Vector of length PK x L
#' @export
algorithm2 <- function(Y,
                       Z,
                       mi_vec,
                       i_index,
                       lp,
                       B,
                       beta, D, gammas, phi, C){
  L <- ncol(Z) - 1
  K <- ncol(Y)
  lp_minus <- NULL

  # Loop only through non-clustering values
  l_loop <- setdiff(0:L , lp)

  for (l in l_loop) {
    gamma_l <- gammas[l + 1]
    beta_l <- beta[, l + 1]


    # Things to be calculated once per loop/beta update
    # Which are: alpha, V inverse, and partials
    alpha <- get_alpha_list(beta = beta,
                            Z = Z,
                            B = B,
                            K = K,
                            i_index = i_index,
                            mi_vec = mi_vec)
    # Get V inverse for all is (and ls... )
    V_inv <- get_V_inv(Y = Y,
                       alpha = alpha,
                       mi_vec = mi_vec,
                       i_index = i_index,
                       phi = phi,
                       beta = beta,
                       Z = Z,
                       B = B,
                       K = ncol(Y))
    partials_l <- get_partials_l_list(Y = Y,
                                      l = l,
                                      mi_vec = mi_vec,
                                      i_index = i_index,
                                      beta = beta,
                                      Z = Z,
                                      B = B,
                                      alpha = alpha)



    # function of beta_l_list not beta_l
    gradient_l <- get_gradient_l(Y = Y,
                                 mi_vec = mi_vec,
                                 i_index = i_index,
                                 l = l,
                                 phi = phi,
                                 beta = beta,
                                 Z = Z,
                                 B = B,
                                 V_inv = V_inv,
                                 partials_l = partials_l) # function of beta_l_list not beta_l
    Hessian_l <- get_Hessian_l(l = l,
                               Y = Y,
                               mi_vec = mi_vec,
                               i_index = i_index,
                               beta = beta,
                               Z = Z,
                               B = B,
                               phi = phi,
                               C = C,
                               V_inv = V_inv,
                               partials_l = partials_l,
                               alpha = alpha)
    # function of beta_l_list not beta_l

    first_term <- -Hessian_l + gamma_l*D
    second_term <- gradient_l - Hessian_l %*% beta_l

    beta_l_s <- MASS::ginv(first_term) %*% second_term

    # Update:
    beta[,l + 1] <- beta_l_s
  }

  return(beta)
}

