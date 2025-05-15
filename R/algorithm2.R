#' Algorithm 2: Updating Non-Clustering coefficients
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
#' @param Y0 Vector of total count for each sample
#' @param max_2_iter Maximum number of iterations for algorithm 2 to run each loop
#' @param epsilon_2 Tolerance for convergence of algorithm 2
#' @param time either a vector of length(Y) or a column reference if Y is a data frame. Must be numeric
#' @param s Current iteration of algorithm 1
#' @param L Number of external variables
#' @param P Number of B-spline coefficients (order + nknots)
#' @param K Number of responses
#' @param M Number of samples times timepoints for each sample
#' @param i_index starting index of the ith subject in the data
#' @param cor_blocks
#' @param j1_j2_list
#' @param cor_str
#'
#' @returns Vector of length PK x L
#' @export
algorithm2 <- function(Y,
                       Y0,
                       Z,
                       mi_vec,
                       i_index,
                       lp,
                       B,
                       beta,
                       D,
                       gammas,
                       phi,
                       C,
                       max_2_iter,
                       epsilon_2,
                       time,
                       s,
                       L,
                       P,
                       K,
                       M,
                       cor_str,
                       cor_blocks,
                       j1_j2_list){

  lp_minus <- NULL
  beta_diffs <- list()

  # Loop only through non-clustering values
  l_loop <- setdiff(0:L , lp)
  # Initialize progress bar
  # cat(paste0("Algorithm 2 iteration: ", s, "\n"))
  # pb_alg2 <- utils::txtProgressBar(min = 0,
  #                                  max = max_2_iter,
  #                                  style = 3)

  for (r in 1:max_2_iter) {
    beta_old <- beta
    for (l in l_loop) {
      gamma_l <- gammas[l + 1]

      beta_l <- beta[, l + 1, drop = F]


      # Things to be calculated once per loop/beta update
      # Which are: alpha, V inverse, and partials
      alpha <- get_alpha_list(beta = beta,
                              Z = Z,
                              B = B,
                              K = K,
                              i_index = i_index,
                              mi_vec = mi_vec,
                              L = L,
                              P = P)
      # Calculate pearson residuals (used to calculate phi and rho)
      # Creates a list of length i
      pearson_residuals <- get_pearson_residuals(Y = Y,
                                                 Y0 = Y0,
                                                 beta = beta,
                                                 alpha = alpha,
                                                 Z = Z,
                                                 B = B,
                                                 K = K,
                                                 mi_vec = mi_vec,
                                                 i_index = i_index,
                                                 M = M)
      # Update dispersion parameter
      phi <- get_phi(pearson_residuals = pearson_residuals,
                     K = K,
                     M = M,
                     L = L,
                     P = P)

      # Update correlation parameter
      rho_cor <- get_rho(pearson_residuals = pearson_residuals,
                         phi = phi,
                         K = K,
                         mi_vec = mi_vec,
                         M = M,
                         cor_str = cor_str,
                         cor_blocks = cor_blocks,
                         j1_j2_list = j1_j2_list)


      # Get V inverse for all is (and ls... )
      V_inv <- get_V_inv(Y = Y,
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
      partials_l <- get_partials_l_list(Y0 = Y0,
                                        l = l,
                                        mi_vec = mi_vec,
                                        i_index = i_index,
                                        beta = beta,
                                        Z = Z,
                                        B = B,
                                        alpha = alpha,
                                        K = K,
                                        P = P)

      # functions of beta_l_list not beta_l
      gradient_l <- get_gradient_l(Y = Y,
                                   Y0 = Y0,
                                   mi_vec = mi_vec,
                                   i_index = i_index,
                                   l = l,
                                   phi = phi,
                                   beta = beta,
                                   alpha = alpha,
                                   Z = Z,
                                   B = B,
                                   V_inv = V_inv,
                                   partials_l = partials_l,
                                   K = K,
                                   P = P)
      Hessian_l <- get_Hessian_l(l = l,
                                 Y0 = Y0,
                                 mi_vec = mi_vec,
                                 i_index = i_index,
                                 beta = beta,
                                 Z = Z,
                                 B = B,
                                 phi = phi,
                                 C = C,
                                 V_inv = V_inv,
                                 partials_l = partials_l,
                                 alpha = alpha,
                                 P = P,
                                 K = K)

      # # function of beta_l_list not beta_l
      # first_term <- -Hessian_l + gamma_l*D
      # first_term_inv <- MASS::ginv(first_term)
      # #second_term <- gradient_l - Hessian_l %*% beta_l
      # second_term <- gradient_l - fast_mat_mult2(Hessian_l,beta_l)
      # #beta_l_s <- Matrix::Matrix(first_term_inv, sparse = T) %*% second_term
      # #beta_l_s <- MASS::ginv(first_term) %*% second_term
      # beta_l_s <- fast_mat_mult2(MASS::ginv(first_term), second_term)
      # Use C++ to calculate the update for beta_l
      beta_l_s <- calculate_alg2(H = Hessian_l,
                                 gamma = gamma_l,
                                 D = D,
                                 Q = gradient_l,
                                 beta_l = beta_l)


      # Update:
      beta[,l + 1] <- as.numeric(beta_l_s)
    }

    beta_diff <- sum((beta - beta_old)^2)
    beta_diffs[[r]] <- beta_diff
    if (beta_diff < epsilon_2) {
      break
    }
    # utils::setTxtProgressBar(pb_alg2, r)
  }
  # utils::setTxtProgressBar(pb_alg2, max_2_iter)
  # close(pb_alg2)

  return(list(beta = beta,
              r = r,
              beta_diffs = beta_diffs))
}

