#include <RcppArmadillo.h>  // Required for Armadillo support
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat fast_mat_mult(const arma::mat& A, const arma::mat& B) {
  return A * B;  // Uses Armadillo's optimized matrix multiplication
}
