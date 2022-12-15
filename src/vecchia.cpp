
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

/*
 * Code derived from GPvecchia package (Katzfuss et al.)
 */

#include "cov.h"
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::mat rev_matrix(arma::mat x) {
  return reverse(x, 1);
}

double d2_vector(arma::rowvec x1, arma::rowvec x2) { 
  int n = x1.size();
  double d2 = 0.0;
  for(int k = 0; k < n; k++) {
    d2 += (x1[k] - x2[k]) * (x1[k] - x2[k]);
  }
  return d2;
}

arma::mat d2_matrix(arma::mat x) {
  int outrows = x.n_rows;
  int outcols = x.n_rows;
  arma::mat d2(outrows, outcols);
  for (int i = 0 ; i < outrows ; i++) {
    for (int j = 0 ; j < outcols ; j++) {
      d2(i, j) = d2_vector(x.row(i), x.row(j));
    }
  }
  return d2;
}

// [[Rcpp::export]]
arma::mat U_entries (const int Ncores, const arma::mat& x, 
                     const arma::umat& revNNarray, const arma::mat& revCondOnLatent, 
                     const double tau2, const double theta, const double g, const double v){
  
  const int m = revNNarray.n_cols - 1;
  const int n = x.n_rows;
  arma::mat Lentries = zeros(n, m + 1);
  
#ifdef _OPENMP
  
#pragma omp parallel for num_threads(Ncores) shared(Lentries) schedule(static)
  for (int k = 0; k < n; k++) {
    arma::uvec inds = revNNarray.row(k).t();
    arma::vec revCon_row = revCondOnLatent.row(k).t();
    arma::uvec inds00 = inds.elem(find(inds)) - 1;
    uword n0 = inds00.n_elem;
    arma::mat dist = d2_matrix(x.rows(inds00));
    arma::mat covmat = Matern(dist, tau2, theta, g, v);
    arma::vec onevec = zeros(n0);
    onevec[n0 - 1] = 1;
    arma::vec M = solve(chol(covmat, "upper"), onevec);
    Lentries(k, span(0, n0 - 1)) = M.t();
  }
  
#else
  
  for (int k = 0; k < n; k++) {
    arma::uvec inds = revNNarray.row(k).t();
    arma::vec revCon_row = revCondOnLatent.row(k).t();
    arma::uvec inds00 = inds.elem(find(inds)) - 1;
    uword n0 = inds00.n_elem;
    arma::mat dist = d2_matrix(x.rows(inds00));
    arma::mat covmat = Matern(dist, tau2, theta, g, v);
    arma::vec onevec = zeros(n0);
    onevec[n0 - 1] = 1;
    arma::vec M = solve(chol(covmat, "upper"), onevec);
    Lentries(k, span(0, n0 - 1)) = M.t();
  }
  
#endif
  
  return Lentries;
}

// [[Rcpp::export]]
arma::mat U_entries_sep (const int Ncores, const arma::mat& x, 
                     const arma::umat& revNNarray, const arma::mat& revCondOnLatent, 
                     const double tau2, const arma::vec theta, const double g, const double v){
  
  const int m = revNNarray.n_cols - 1;
  const int n = x.n_rows;
  arma::mat Lentries = zeros(n, m + 1);
  
#ifdef _OPENMP
  
#pragma omp parallel for num_threads(Ncores) shared(Lentries) schedule(static)
  for (int k = 0; k < n; k++) {
    arma::uvec inds = revNNarray.row(k).t();
    arma::vec revCon_row = revCondOnLatent.row(k).t();
    arma::uvec inds00 = inds.elem(find(inds)) - 1;
    uword n0 = inds00.n_elem;
    arma::mat covmat = MaternSep(x.rows(inds00), x.rows(inds00), tau2, theta, g, v);
    arma::vec onevec = zeros(n0);
    onevec[n0 - 1] = 1;
    arma::vec M = solve(chol(covmat, "upper"), onevec);
    Lentries(k, span(0, n0 - 1)) = M.t();
  }
  
#else
  
  for (int k = 0; k < n; k++) {
    arma::uvec inds = revNNarray.row(k).t();
    arma::vec revCon_row = revCondOnLatent.row(k).t();
    arma::uvec inds00 = inds.elem(find(inds)) - 1;
    uword n0 = inds00.n_elem;
    arma::mat covmat = MaternSep(x.rows(inds00), x.rows(inds00), tau2, theta, g, v);
    arma::vec onevec = zeros(n0);
    onevec[n0 - 1] = 1;
    arma::vec M = solve(chol(covmat, "upper"), onevec);
    Lentries(k, span(0, n0 - 1)) = M.t();
  }
  
#endif
  
  return Lentries;
}


// [[Rcpp::export]]
void check_omp () {
  #ifdef _OPENMP
    // DO NOTHING
  #else 
    Rcout << "NOTE: OpenMP install suggested for best results; see ?fit_two_layer for details \n";
  #endif
}
  
