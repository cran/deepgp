
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
NumericVector forward_solve_raw(NumericMatrix U, NumericVector z,
                            NumericMatrix NNarray) {
  // Solves U * y = z for y
  // Uses raw form of U (in create_U use raw_form = TRUE)
  
  int n = U.nrow();
  NumericVector y(n);
  int mp1 = NNarray.ncol(); // m plus 1

  y(0) = z(0) / U(0, 0);
  
  for (int i = 1; i < n; i++) {
    int B = min(i + 1, mp1);
    y(i) = z(i);
    for (int j = 1; j < B; j++) {
      y(i) -= U(i, j) * y(NNarray(i, j) - 1);
    }
    y(i) = y(i) / U(i, 0);
  }
  return y;
}

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
arma::mat U_entries (const int Ncores, const arma::mat& x, const arma::umat& revNNarray,
                     const double tau2, const double theta, const double g, const double v){
  
  const int m = revNNarray.n_cols - 1;
  const int n = x.n_rows;
  arma::mat Lentries = zeros(n, m + 1);
  arma::mat covmat;
  
  #ifdef _OPENMP
    
  #pragma omp parallel for num_threads(Ncores) shared(Lentries) schedule(static)
    for (int k = 0; k < n; k++) {
      arma::uvec inds = revNNarray.row(k).t();
      arma::uvec inds00 = inds.elem(find(inds)) - 1;
      uword n0 = inds00.n_elem;
      arma::mat dist = d2_matrix(x.rows(inds00));
      arma::mat covmat(n0, n0);
      if (v == 999) {
        covmat = Exp2(dist, tau2, theta, g);
      } else {
        covmat = Matern(dist, tau2, theta, g, v);
      }
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(k, span(0, n0 - 1)) = M.t();
    }
    
  #else
    
    for (int k = 0; k < n; k++) {
      arma::uvec inds = revNNarray.row(k).t();
      arma::uvec inds00 = inds.elem(find(inds)) - 1;
      uword n0 = inds00.n_elem;
      arma::mat dist = d2_matrix(x.rows(inds00));
      arma::mat covmat(n0, n0);
      if (v == 999) {
        covmat = Exp2(dist, tau2, theta, g);
      } else {
        covmat = Matern(dist, tau2, theta, g, v);
      }
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(k, span(0, n0 - 1)) = M.t();
    }
    
  #endif
  
  return Lentries;
}

// [[Rcpp::export]]
arma::mat U_entries_sep (const int Ncores, const arma::mat& x, const arma::umat& revNNarray, 
                     const double tau2, const arma::vec theta, const double g, const double v){
  
  const int m = revNNarray.n_cols - 1;
  const int n = x.n_rows;
  arma::mat Lentries = zeros(n, m + 1);
  arma::mat covmat;
  
  #ifdef _OPENMP
    
  #pragma omp parallel for num_threads(Ncores) shared(Lentries) schedule(static)
    for (int k = 0; k < n; k++) {
      arma::uvec inds = revNNarray.row(k).t();
      arma::uvec inds00 = inds.elem(find(inds)) - 1;
      uword n0 = inds00.n_elem;
      arma::mat covmat(n0, n0);
      if (v == 999) {
        covmat = Exp2Sep(x.rows(inds00), x.rows(inds00), tau2, theta, g);
      } else {
        covmat = MaternSep(x.rows(inds00), x.rows(inds00), tau2, theta, g, v);
      }
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(k, span(0, n0 - 1)) = M.t();
    }
    
  #else
    
    for (int k = 0; k < n; k++) {
      arma::uvec inds = revNNarray.row(k).t();
      arma::uvec inds00 = inds.elem(find(inds)) - 1;
      uword n0 = inds00.n_elem;
      arma::mat covmat(n0, n0);
      if (v == 999) {
        covmat = Exp2Sep(x.rows(inds00), x.rows(inds00), tau2, theta, g);
      } else {
        covmat = MaternSep(x.rows(inds00), x.rows(inds00), tau2, theta, g, v);
      }
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
  
// [[Rcpp::export]]
arma::mat row_col_pointers(const arma::umat& NNarray) {
    
  const int m = NNarray.n_cols- 1;
  const int n = NNarray.n_rows;
  int start, col_count;
    
  int length = (n - m) * (m + 1);
  for (int i = 1; i <= m; i ++)
    length += i;
    
  arma::mat pointers = zeros(length, 2);
    
  start = 0;
  for (int i = 1; i <= n; i++) {
    if (i <= m) {
      col_count = i - 1;
      for (int j = start; j < start + i; j++) {
        pointers(j, 0) = i;
        pointers(j, 1) = NNarray(i - 1, col_count);
        col_count -= 1;
      }
      start += i;
    } else {
      col_count = m;
      for (int j = start; j < start + m + 1; j++) {
        pointers(j, 0) = i;
        pointers(j, 1) = NNarray(i - 1, col_count);
        col_count -= 1;
      }
      start += m + 1;
    }
  }
  return pointers;
}

