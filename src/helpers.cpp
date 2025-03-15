#include <RcppArmadillo.h>  // Required for Armadillo support
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat fast_mat_mult2(const arma::mat& A, const arma::mat& B) {
  return A * B;  // Uses Armadillo's optimized matrix multiplication
}


// Function to calculate beta_lp_new
// [[Rcpp::export]]
arma::vec calculate_beta_lp_new(const arma::mat& H, double gamma,
                                const arma::mat& D, double theta,
                                const arma::mat& AtA, const arma::mat& Q,
                                const arma::mat& beta_lp, const arma::mat& A,
                                const arma::mat& v_tilde) {

  // Calculate first term
  arma::mat first_term = -H + gamma * D + theta * AtA;

  // Calculate second term
  arma::vec second_term = Q - H * beta_lp + theta * arma::trans(A) * v_tilde;

  // Compute beta_lp_new using the generalized inverse of first_term
  arma::mat first_term_inv = arma::pinv(first_term); // Generalized inverse (Moore-Penrose)
  arma::vec beta_lp_new = first_term_inv * second_term;

  return beta_lp_new;
}


// Function to compute alpha_ijk
// [[Rcpp::export]]
double compute_alpha_ijk(const arma::mat& beta, int k, int P,
                         const arma::vec& Z_ij, const arma::vec& B_ij) {

  // Extract beta_k (subset of beta matrix)
  arma::mat beta_k = beta.rows((k - 1) * P, k * P - 1);

  // Compute (B_ij' * beta_k), which results in a row vector
  arma::rowvec Bbeta = B_ij.t() * beta_k;

  // Element-wise multiplication Z_ij * Bbeta (broadcasting across rows)
  arma::vec lsum = Z_ij % Bbeta.t();  // Convert row vector to column vector

  // Compute alpha_ijk
  double alpha_ijk = std::exp(arma::sum(lsum));

  return alpha_ijk;


}


// Function to compute U_ij
// [[Rcpp::export]]
arma::mat get_U_ij_cpp(const arma::vec& alpha_ij) {

  // Step 1: Compute alpha_ij0 (sum of alpha_ij)
  double alpha_ij0 = arma::sum(alpha_ij);

  // Step 2: Compute the scaled diagonal matrix directly
  arma::mat diag_matrix = arma::diagmat(alpha_ij) / alpha_ij0;

  // Step 3: Compute the outer product (cross product), scale it by 1 / alpha_ij0^2
  arma::mat tcrossprod_matrix = alpha_ij * alpha_ij.t();
  tcrossprod_matrix /= (alpha_ij0 * alpha_ij0);

  // Step 4: Compute U_ij as diag(alpha_ij / alpha_ij0) - tcrossprod(alpha_ij) / alpha_ij0^2
  arma::mat U_ij = diag_matrix - tcrossprod_matrix;

  return U_ij;
}
