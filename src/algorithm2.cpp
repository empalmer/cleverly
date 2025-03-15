#include <RcppArmadillo.h>  // Required for Armadillo support
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Function to calculate beta_lp_new
// [[Rcpp::export]]
arma::vec calculate_alg2(const arma::mat& H,
                         double gamma,
                         const arma::mat& D,
                         const arma::mat& Q,
                         const arma::mat& beta_l) {

  // Calculate first term
  arma::mat first_term = -H + gamma * D;

  // Calculate second term
  arma::vec second_term = Q - H * beta_l;

  // Compute beta_lp_new using the generalized inverse of first_term
  arma::mat first_term_inv = arma::pinv(first_term); // Generalized inverse (Moore-Penrose)
  arma::vec beta_lp_new = first_term_inv * second_term;


  return beta_lp_new;
}
