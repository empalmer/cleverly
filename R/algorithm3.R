#' Update v
#'
#' Update v step of ADMM algorithm. V is the pairwise differences vector of dimension |Kappa|P
#'
#' @param beta
#' @param lp
#' @param lambda
#' @param Kappa
#' @param gammas
#' @param tau
#' @param theta
#' @param psi
#'
#' @returns
#' @export
#'
#' @examples
update_v <- function(beta,
                     lp,
                     lambda,
                     Kappa,
                     P,
                     gammas,
                     tau,
                     theta,
                     psi){

  mcp <- tau * theta / (tau * theta - 1)
  sigma <- psi/theta
  v <- numeric(nrow(Kappa) * P)

  for (kappa in 1:nrow(Kappa)) {
    k1 <- Kappa[kappa, 1]
    k2 <- Kappa[kappa, 2]
    beta_k1 <- beta[(P * (k1 - 1) + 1):(P * k1), lp + 1]
    beta_k2 <- beta[(P * (k2 - 1) + 1):(P * k2), lp + 1]
    lambda_kappa <- lambda[(P * (kappa - 1) + 1):(P * kappa)]
    u <- beta_k1 - beta_k2 - lambda_kappa/theta

    norm_u <- sum(u^2) # or norm(u, type = "2")
    if (norm_u >= tau * psi | norm_u == 0) {
      v[((kappa - 1) * P + 1):(kappa * P)] <- u
    } else {
      v[((kappa - 1) * P + 1):(kappa * P)] <- mcp * max(0, ( 1 - sigma/norm_u)) * u
    }
  }
  return(v)
}


#' Update beta_lp in the ADMM algorithm
#'
#' @param Y
#' @param beta
#' @param lp
#' @param v
#' @param lambda
#' @param theta
#' @param gammas
#' @param D
#' @param A
#' @param Z
#' @param B
#' @param phi
#' @param C
#' @param mi_vec
#'
#' @returns
#' @export
#'
#' @examples
update_beta_admm <- function(Y,
                             beta,
                             lp,
                             v,
                             lambda,
                             theta,
                             gammas,
                             D,
                             A,
                             Z,
                             B,
                             phi,
                             C,
                             mi_vec){

  gamma <- gammas[lp + 1]
  H <- get_Hessian_l(l = lp,
                     Y = Y,
                     mi_vec = mi_vec,
                     beta = beta,
                     Z = Z,
                     B = B,
                     phi = phi,
                     C = C)
  Q <- get_gradient_l(Y = Y,
                      mi_vec = mi_vec,
                      l = lp,
                      phi = phi,
                      beta = beta,
                      Z = Z,
                      B = B)

  v_tilde <- v - lambda/theta
  beta_lp <- beta[,lp + 1]

  first_term <- -H + gamma * D + theta * t(A) %*% A
  second_term <- Q - H %*% beta_lp + theta * t(A) %*% v_tilde


  beta_lp_new <- MASS::ginv(first_term) %*% second_term

  # All other beta elements are fixed, only the lp column is updated.
  beta[,lp + 1] <- beta_lp_new
  return(beta)
}

#' Update lambda step of the admm algorithm
#'
#' @param beta
#' @param v
#' @param lambda
#' @param A
#' @param theta
#'
#' @returns
#' @export
#'
#' @examples
update_lambda <- function(beta_lp, v, lambda, A, theta){
  lambda_new <- lambda + theta * (A %*% beta_lp - v)
  return(lambda_new)
}


#' ADMM algorithm for updating the clustering beta
#'
#' @param Y
#' @param Z
#' @param lp
#' @param B
#' @param beta
#' @param A
#' @param kappa
#' @param gamma
#' @param tau
#' @param theta
#' @param psi
#' @param max_admm_iter
#'
#' @returns
#' @export
#'
#' @examples
algorithm3 <- function(Y,
                       Z,
                       lp,
                       B,
                       beta,
                       A,
                       P,
                       C,
                       D,
                       Kappa,
                       mi_vec,
                       gammas,
                       tau,
                       theta,
                       psi,
                       phi,
                       max_admm_iter = 100) {
  lambda <- numeric(nrow(Kappa)*P)

  t <- 1


  beta_admm_track <- list()
  v_admm_track <- list()
  lambda_admm_track <- list()

  while (t <= max_admm_iter) {
    v_new <- update_v(beta = beta,
                      lp = lp,
                      lambda = lambda,
                      Kappa = Kappa,
                      P = P,
                      gammas = gammas,
                      tau = tau,
                      theta = theta,
                      psi = psi)
    beta_new <- update_beta_admm(Y = Y,
                                 beta = beta,
                                 lp = lp,
                                 v = v_new,
                                 lambda = lambda,
                                 theta = theta,
                                 gammas = gammas,
                                 D = D,
                                 A = A,
                                 Z = Z,
                                 B = B,
                                 phi = phi,
                                 C = C,
                                 mi_vec = mi_vec)
    lambda_new <- update_lambda(beta_lp = beta_new[,lp + 1],
                                v = v_new,
                                lambda = lambda,
                                A = A,
                                theta = theta)


    beta_admm_track[[t]] <- beta_new
    v_admm_track[[t]] <- v_new
    lambda_admm_track[[t]] <- lambda_new
    # Check for convergence
    converged <- FALSE
    if (converged) {
      break
    }
    # Prepare for next iteration
    v <- v_new
    beta <- beta_new
    lambda <- lambda_new
    t <- t + 1
  }

  return(list(beta = beta,
              lambda = lambda,
              v = v,
              beta_admm_track = beta_admm_track,
              lambda_admm_track = lambda_admm_track,
              v_admm_track = v_admm_track))
}
