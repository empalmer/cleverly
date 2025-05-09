# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

calculate_alg2 <- function(H, gamma, D, Q, beta_l) {
    .Call(`_cleverly_calculate_alg2`, H, gamma, D, Q, beta_l)
}

fast_mat_mult2 <- function(A, B) {
    .Call(`_cleverly_fast_mat_mult2`, A, B)
}

calculate_beta_lp_new <- function(H, gamma, D, theta, AtA, Q, beta_lp, A, v_tilde) {
    .Call(`_cleverly_calculate_beta_lp_new`, H, gamma, D, theta, AtA, Q, beta_lp, A, v_tilde)
}

compute_alpha_ijk <- function(beta, k, P, Z_ij, B_ij) {
    .Call(`_cleverly_compute_alpha_ijk`, beta, k, P, Z_ij, B_ij)
}

get_U_ij_cpp <- function(alpha_ij) {
    .Call(`_cleverly_get_U_ij_cpp`, alpha_ij)
}

fast_mat_mult3 <- function(A, B, C) {
    .Call(`_cleverly_fast_mat_mult3`, A, B, C)
}

