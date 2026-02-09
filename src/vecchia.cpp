
/*
 * This file contains functions for vecchia calculations
 * THESE FUNCTIONS HAVE BEEN STREAMLINED FOR VERSION 1.2.0
 * Internal:
 *    d2_vector
 *    d2_matrix
 * External:
 *    forward_solve_raw
 *    U_entries
 *    U_entries_grad
 *    U_entries_sep
 *    U_entries_sep_grad
 *    check_omp
 *    row_col_pointers
 */    

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

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
arma::vec forward_solve_raw(const arma::mat& U, const arma::vec z,
                            const arma::umat& NNarray) {
  // Solves t(U) * y = z for y
  // Uses original form of U from the U_entries function which has NOT been 
  // turned into a sparse matrix with row_col_pointers yet
  
  int n = U.n_cols;
  arma::vec y = zeros(n);
  int mp1 = NNarray.n_cols; // m plus 1
  arma::vec Urow_with_na = zeros(mp1);
  
  y(0) = z(0) / U(0, 0);
  
  for (int i = 1; i < n; i++) {
    Urow_with_na = reverse(U.col(i)); // grab the ith column, NA values will be leading zeros
    arma::vec Urow = Urow_with_na.elem(find(Urow_with_na)); // finds nonzero entries
    int B = min(i + 1, mp1);
    y(i) = z(i);
    for (int j = 1; j < B; j++) {
      y(i) -= Urow(j) * y(NNarray(i, j) - 1);
    }
    y(i) = y(i) / Urow(0);
  }
  return y;
}

double d2_vector(arma::rowvec x1, arma::rowvec x2) { 
  // Gets squared distances between two vectors
  int n = x1.size();
  double d2 = 0.0;
  for(int k = 0; k < n; k++) {
    d2 += (x1[k] - x2[k]) * (x1[k] - x2[k]);
  }
  return d2;
}

arma::mat d2_matrix(arma::mat x) {
  // Gets pairwise squared distances for all elements of a matrix
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
arma::mat U_entries(const int cores, const arma::mat& x, const arma::umat& NNarray,
                     const double tau2, const double theta, const double g, const double v) {
  
  // This function has been derived from that of the GPvecchia package (Katzfuss et al.)
  // and the GpGP package (Guinness et al.)
  // Note: x must be ordered to align with NNarray
  
  const int m = NNarray.n_cols - 1;
  const int n = x.n_rows;
  arma::mat Lentries = zeros(n, m + 1);

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(cores) shared(Lentries) schedule(static)
    for (int i = 0; i < n; i++) {
      // Grab indices of conditioning set for observation i
      arma::uvec inds_with_na = reverse(NNarray.row(i)).t(); // turns NA into zero
      arma::uvec inds = inds_with_na.elem(find(inds_with_na)); // finds nonzero entries
      inds = inds - 1; // convert R indexing which starts at 1 to cpp indexing which starts at 0
      uword n0 = inds.n_elem; // 1 + size of conditioning set for observation i
      
      // Get covariance matrix for point i and its conditioning set
      arma::mat dist = d2_matrix(x.rows(inds));
      arma::mat covmat(n0, n0);
      if (v == 999) {
        covmat = Exp2(dist, tau2, theta, g);
      } else {
        covmat = Matern(dist, tau2, theta, g, v);
      }
      
      // Solve for the entries of L
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(i, arma::span(0, n0 - 1)) = M.t();
    }
  #else
    for (int i = 0; i < n; i++) {
      // Grab indices of conditioning set for observation i
      arma::uvec inds_with_na = reverse(NNarray.row(i)).t(); // turns NA into zero
      arma::uvec inds = inds_with_na.elem(find(inds_with_na)); // finds nonzero entries
      inds = inds - 1; // convert R indexing which starts at 1 to cpp indexing which starts at 0
      uword n0 = inds.n_elem; // 1 + size of conditioning set for observation i
      
      // Get covariance matrix for point i and its conditioning set
      arma::mat dist = d2_matrix(x.rows(inds));
      arma::mat covmat(n0, n0);
      if (v == 999) {
        covmat = Exp2(dist, tau2, theta, g);
      } else {
        covmat = Matern(dist, tau2, theta, g, v);
      }
      
      // Solve for the entries of L
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(i, arma::span(0, n0 - 1)) = M.t();
    }
  #endif
  
  return Lentries.t(); // U = transpose of L
}

// [[Rcpp::export]]
arma::mat U_entries_grad(const int cores, const arma::mat& x, const arma::umat& NNarray,
                         const arma::vec grad, const double tau2, 
                         const double theta, const double g, const double v) {
  
  // ONLY implemented for "Exp2" with v = 999
  
  const int m = NNarray.n_cols - 1;
  const int n = x.n_rows;
  arma::mat Lentries = zeros(n, m + 1);
  
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(cores) shared(Lentries) schedule(static)
    for (int i = 0; i < n; i++) {
      // Grab indices of conditioning set for observation i
      arma::uvec inds_with_na = reverse(NNarray.row(i)).t(); // turns NA into zero
      arma::uvec inds = inds_with_na.elem(find(inds_with_na)); // finds nonzero entries
      inds = inds - 1; // convert R indexing which starts at 1 to cpp indexing which starts at 0
      uword n0 = inds.n_elem; // 1 + size of conditioning set for observation i
      
      // Get covariance matrix for point i and its conditioning set
      arma::mat covmat(n0, n0);
      covmat = Exp2Grad(x.rows(inds), x.rows(inds), grad(inds), grad(inds), tau2, theta, g);
      
      // Solve for the entries of L
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(i, arma::span(0, n0 - 1)) = M.t();
    }
  #else
    for (int i = 0; i < n; i++) {
      // Grab indices of conditioning set for observation i
      arma::uvec inds_with_na = reverse(NNarray.row(i)).t(); // turns NA into zero
      arma::uvec inds = inds_with_na.elem(find(inds_with_na)); // finds nonzero entries
      inds = inds - 1; // convert R indexing which starts at 1 to cpp indexing which starts at 0
      uword n0 = inds.n_elem; // 1 + size of conditioning set for observation i
      
      // Get covariance matrix for point i and its conditioning set
      arma::mat covmat(n0, n0);
      covmat = Exp2Grad(x.rows(inds), x.rows(inds), grad(inds), grad(inds), tau2, theta, g);
      
      // Solve for the entries of L
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(i, arma::span(0, n0 - 1)) = M.t();
    }
  #endif
  
  return Lentries.t(); // U = transpose of L
}

// [[Rcpp::export]]
arma::mat U_entries_sep(const int cores, const arma::mat& x, const arma::umat& NNarray,
                     const double tau2, const arma::vec theta, const double g, const double v) {
  
  // Same as U_entries function, but accepts a vector of theta instead of a scalar

  const int m = NNarray.n_cols - 1;
  const int n = x.n_rows;
  arma::mat Lentries = zeros(n, m + 1);

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(cores) shared(Lentries) schedule(static)
    for (int i = 0; i < n; i++) {
      // Grab indices of conditioning set for observation i
      arma::uvec inds_with_na = reverse(NNarray.row(i)).t(); // turns NA into zero
      arma::uvec inds = inds_with_na.elem(find(inds_with_na)); // finds nonzero entries
      inds = inds - 1; // convert R indexing which starts at 1 to cpp indexing which starts at 0
      uword n0 = inds.n_elem; // 1 + size of conditioning set for observation i
      
      // Get covariance matrix for point i and its conditioning set
      arma::mat covmat(n0, n0);
      if (v == 999) {
        covmat = Exp2Sep(x.rows(inds), x.rows(inds), tau2, theta, g);
      } else {
        covmat = MaternSep(x.rows(inds), x.rows(inds), tau2, theta, g, v);
      }
      
      // Solve for the entries of L
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(i, arma::span(0, n0 - 1)) = M.t();
    }
  #else
    for (int i = 0; i < n; i++) {
      // Grab indices of conditioning set for observation i
      arma::uvec inds_with_na = reverse(NNarray.row(i)).t(); // turns NA into zero
      arma::uvec inds = inds_with_na.elem(find(inds_with_na)); // finds nonzero entries
      inds = inds - 1; // convert R indexing which starts at 1 to cpp indexing which starts at 0
      uword n0 = inds.n_elem; // 1 + size of conditioning set for observation i
      
      // Get covariance matrix for point i and its conditioning set
      arma::mat covmat(n0, n0);
      if (v == 999) {
        covmat = Exp2Sep(x.rows(inds), x.rows(inds), tau2, theta, g);
      } else {
        covmat = MaternSep(x.rows(inds), x.rows(inds), tau2, theta, g, v);
      }
      
      // Solve for the entries of L
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(i, arma::span(0, n0 - 1)) = M.t();
    }
  #endif
  
  return Lentries.t(); // U = transpose of L
}

// [[Rcpp::export]]
arma::mat U_entries_sep_grad(const int cores, const arma::mat& x, const arma::umat& NNarray,
                             const arma::vec grad, const double tau2, 
                             const arma::vec theta, const double g, const double v) {
  
  // ONLY implemented for "Exp2" with v = 999
  // Same as U_entries_grad function, but accepts a vector of theta instead of a scalar
  
  const int m = NNarray.n_cols - 1;
  const int n = x.n_rows;
  arma::mat Lentries = zeros(n, m + 1);
  
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(cores) shared(Lentries) schedule(static)
    for (int i = 0; i < n; i++) {
      // Grab indices of conditioning set for observation i
      arma::uvec inds_with_na = reverse(NNarray.row(i)).t(); // turns NA into zero
      arma::uvec inds = inds_with_na.elem(find(inds_with_na)); // finds nonzero entries
      inds = inds - 1; // convert R indexing which starts at 1 to cpp indexing which starts at 0
      uword n0 = inds.n_elem; // 1 + size of conditioning set for observation i
      
      // Get covariance matrix for point i and its conditioning set
      arma::mat covmat(n0, n0);
      covmat = Exp2SepGrad(x.rows(inds), x.rows(inds), grad(inds), grad(inds), tau2, theta, g);
      
      // Solve for the entries of L
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(i, arma::span(0, n0 - 1)) = M.t();
    }
  #else
    for (int i = 0; i < n; i++) {
      // Grab indices of conditioning set for observation i
      arma::uvec inds_with_na = reverse(NNarray.row(i)).t(); // turns NA into zero
      arma::uvec inds = inds_with_na.elem(find(inds_with_na)); // finds nonzero entries
      inds = inds - 1; // convert R indexing which starts at 1 to cpp indexing which starts at 0
      uword n0 = inds.n_elem; // 1 + size of conditioning set for observation i
      
      // Get covariance matrix for point i and its conditioning set
      arma::mat covmat(n0, n0);
      covmat = Exp2SepGrad(x.rows(inds), x.rows(inds), grad(inds), grad(inds), tau2, theta, g);
      
      // Solve for the entries of L
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(i, arma::span(0, n0 - 1)) = M.t();
    }
  #endif
  
  return Lentries.t(); // U = transpose of L
}

// [[Rcpp::export]]
void check_omp() {
  #ifdef _OPENMP
    // DO NOTHING
  #else 
    Rcout << "NOTE: OpenMP install suggested for best results; see ?fit_two_layer for details \n";
  #endif
}
  
// [[Rcpp::export]]
arma::mat row_col_pointers(const arma::umat& NNarray) {
    
  // Returns row and column pointers for the entries of U to be placed in a sparse matrix
  
  const int m = NNarray.n_cols - 1;
  const int n = NNarray.n_rows;
  int start, col_count;
    
  int length = (n-m) * (m+1);
  for (int i = 1; i <= m; i ++)
    length += i;
    
  arma::mat pointers = zeros(length, 2);
    
  start = 0;
  for (int i = 1; i <= n; i++) {
    if (i <= m) {
      col_count = i - 1;
      for (int j = start; j < start + i; j++) {
        pointers(j, 0) = NNarray(i - 1, col_count);
        pointers(j, 1) = i;
        col_count -= 1;
      }
      start += i;
    } else {
      col_count = m;
      for (int j = start; j < start + m + 1; j++) {
        pointers(j, 0) = NNarray(i - 1, col_count);
        pointers(j, 1) = i;
        col_count -= 1;
      }
      start += m + 1;
    }
  }
  return pointers;
}

