#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::colvec trajectory(arma::mat G, Rcpp::NumericVector current_position, Function W, Function grad, double diff_cut_off, int max_gens, int max_stand_still) {
  int k;
  arma::colvec z(current_position.begin(), current_position.size(), false);
  Rcpp::NumericVector beta;
  for(k = 0; k < max_gens; ++k){
    beta = grad(W, NumericVector(z.begin(),z.end()));
    arma::colvec b(beta.begin(), beta.size(), false);
    z = G * b + z;
  }
  return z;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
diff_cut_off = 1e-1
max_gens = 100
max_stand_still = 10
G = RandomMatrix(2)
init = c(0, 0)
random_peaks = matrix(runif(n_traits*n_peaks, -8, 8), n_peaks, 2, byrow = T)
W_bar = W_bar_factory(random_peaks)
grad(W_bar, t(init))
trajectory(G, init, W_bar, grad, diff_cut_off, max_gens, max_stand_still)
*/
