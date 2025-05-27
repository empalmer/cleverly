
#' CLR cluster
#'
#' @param y
#' @param Z
#' @param beta
#' @param lp
#' @param lp_minus
#' @param B
#' @param clusters
#' @param K
#' @param P
#' @param M
#'
#' @returns
#' @export
CLR_cluster <- function(y, Z, lp,  lp_minus, B, K, P, M) {

  L <- ncol(Z) - 1


  gamma <- matrix(nrow = P * (L + 1), ncol = K)

  ZB_list <- list()
  for (l in 0:L) {
    Z_l <- diag(Z[,l + 1])
    ZB_list[[l + 1]] <- Z_l %*% B
  }

  y_clr <- compositions::clr(y)



  for(k in 1:K)



  for (c in seq_along(clusters$csize)) {

    df_list <- list()
    beta_names <- c()
    for (l in 0:L) {

    df <- do.call(cbind, df_list)

    log_y_class <- as.vector(log(y[, class_indices] + 1))



    mod_g <- lm(log_y_class ~ 0 + df)
    coefs <- coef(mod_g)
    beta_c_mat <- matrix(coefs, nrow = P)



    col_index <- 1

    for (l in 0:L) {
      if (l == lp) {
        # Shared coefficients across all in cluster
        beta_group[(lp * P + 1):((lp + 1) * P), class_indices] <- beta_c_mat[, col_index]
        col_index <- col_index + 1
      } else {
        for (k in seq_along(class_indices)) {
          col_k <- col_index
          idx_k <- class_indices[k]
          beta_group[(l * P + 1):((l + 1) * P), idx_k] <- beta_c_mat[, col_k]
          col_index <- col_index + 1
        }
      }
    }

  }


  # Convert back to original structure
  beta_array_back <- array(beta_group, dim = c(P, (L + 1), K))
  beta_array_orig <- aperm(beta_array_back, c(1, 3, 2))
  beta_in <- matrix(beta_array_orig, nrow = P * K, ncol = (L + 1))

  return(beta_in)
}

}
