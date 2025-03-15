#include <RcppArmadillo.h>  // Required for Armadillo support
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat fast_mat_mult3(const arma::mat& A, const arma::mat& B, const arma::mat& C) {
  return A * B * C;  // Uses Armadillo's optimized matrix multiplication
}
