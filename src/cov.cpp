
/*
 * This file contains functions for kernel calculations
 * External:
 *    diag_quad_mat
 *    Exp2
 *    Exp2Grad
 *    Exp2Sep
 *    Exp2SepGrad
 *    Matern
 *    MaternSep
 */


// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <math.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
arma::vec diag_quad_mat(const arma::mat A, const arma::mat B) {
  // Returns diag(A %*% B %*% t(A)) for s2 calculations
  int Arow = A.n_rows;
  int Brow = B.n_rows;
  arma::vec s2(Arow);
  for (int i = 0; i < Arow; i++) {
    s2(i) = 0.0;
    for (int j = 0; j < Brow; j++) {
      double temp_sum = 0.0;
      for (int n = 0; n < Brow; n++)
        temp_sum += A(i, n) * B(n, j);
      s2(i) += temp_sum * A(i, j);
    }
  }
  return s2;
}

// [[Rcpp::export]]
arma::mat Exp2(const arma::mat distmat, const double tau2, const double theta,
                  const double g) { 
  // distmat = matrix of SQUARED distances
  int n1 = distmat.n_rows;
  int n2 = distmat.n_cols;
  arma::mat covmat(n1, n2);
  double r;
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      r = distmat(i, j) / theta;
      covmat(i, j) = tau2 * exp(-r);
    }
  }
  if (n1 == n2) {
    for (int i = 0; i < n1; i++) 
      covmat(i, i) += tau2 * g;
  }
  return covmat;
}

// [[Rcpp::export]]
arma::mat Exp2Grad(const arma::mat x1, const arma::mat x2, const arma::vec grad1, 
                   const arma::vec grad2, const double tau2, 
                   const double theta, const double g) { 
  // grad1 = integers indicating which type of observations each x1 value is
  // grad2 = integers indicating which type of observations each x2 value is
  // 0 is normal observation, 1 is gradient with respect to dimension 1, etc.
  if (x1.n_cols != x2.n_cols) 
    stop("dimension of x1 and x2 do not match");
  if (grad1.n_elem != x1.n_rows)
    stop("length of grad1 does not match dimension of x1");
  if (grad2.n_elem != x2.n_rows)
    stop("length of grad2 does not match dimension of x2");
  int n1 = x1.n_rows;
  int n2 = x2.n_rows;
  int d = x1.n_cols;
  arma::mat covmat(n1, n2);
  double r;
  double K_00;
  int d_index;
  int f_index;
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      // First get Euclidean distance
      r = 0.0;
      for (int k = 0; k < d; k++) {
        r += (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k));
      }
      // Then get K_00 (all methods use it)
      K_00 = exp(-r/theta);
      
      if (grad1(i) == 0 and grad2(j) == 0) { // K_00
        // Get Euclidean distance
        covmat(i, j) = tau2 * K_00;
      } else if (grad1(i) == 0 and grad2(j) != 0) { // K_0d
        d_index = grad2(j)-1; // change indexing to start at zero
        covmat(i, j) = tau2 * K_00 * (2/theta) * (x1(i, d_index) - x2(j, d_index));
      } else if (grad1(i) != 0 and grad2(j) == 0) { // K_d0
        d_index = grad1(i)-1; // change indexing to start at zero
        covmat(i, j) = tau2 * K_00 * (-2/theta) * (x1(i, d_index) - x2(j, d_index));
      } else if (grad1(i) == grad2(j)) { // K_dd
        d_index = grad1(i)-1; // change indexing to start at zero
        covmat(i, j) = tau2 * K_00 * (2/theta) * (1 - (2/theta)*(x1(i, d_index) - x2(j, d_index))*(x1(i, d_index) - x2(j, d_index)));
      } else { // K_df
        d_index = grad1(i)-1; // change indexing to start at zero
        f_index = grad2(j)-1; // change indexing to start at zero
        covmat(i, j) = tau2 * K_00 * (-4/(theta * theta))*(x1(i, d_index) - x2(j, d_index))*(x1(i, f_index) - x2(j, f_index));
      }
    }
  }
  // Add nugget to all diagonal elements (acts as jitter for numerical stability)
  if (n1 == n2) {
    for (int i = 0; i < n1; i++) {
      covmat(i, i) += tau2*g;
    }
  }
  return covmat;
}

// [[Rcpp::export]]
arma::mat Exp2Sep(const arma::mat x1, const arma::mat x2, const double tau2, 
                     const arma::vec theta, const double g) { 
  if (x1.n_cols != x2.n_cols) 
    stop("dimension of x1 and x2 do not match");
  if (theta.n_elem != x1.n_cols)
    stop("length of theta does not match dimension of x");
  int n1 = x1.n_rows;
  int n2 = x2.n_rows;
  int d = x1.n_cols;
  arma::mat covmat(n1, n2);
  double r;
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      r = 0.0;
      for (int k = 0; k < d; k++)
        r += (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k);
      covmat(i, j) = tau2 * exp(-r);
    }
  }
  if (n1 == n2) {
    for (int i = 0; i < n1; i++) 
      covmat(i, i) += tau2 * g;
  }
  return covmat;
}

// [[Rcpp::export]]
arma::mat Exp2SepGrad(const arma::mat x1, const arma::mat x2, const arma::vec grad1,
                      const arma::vec grad2, const double tau2, 
                      const arma::vec theta, const double g) { 
  // grad1 = integers indicating which type of observations each x1 value is
  // grad2 = integers indicating which type of observations each x2 value is
  // 0 is normal observation, 1 is gradient with respect to dimension 1, etc.
  if (x1.n_cols != x2.n_cols) 
    stop("dimension of x1 and x2 do not match");
  if (grad1.n_elem != x1.n_rows)
    stop("length of grad1 does not match dimension of x1");
  if (grad2.n_elem != x2.n_rows)
    stop("length of grad2 does not match dimension of x2");
  if (theta.n_elem != x1.n_cols)
    stop("length of theta does not match dimension of x");
  int n1 = x1.n_rows;
  int n2 = x2.n_rows;
  int d = x1.n_cols;
  arma::mat covmat(n1, n2);
  double r;
  double K_00;
  int d_index;
  int f_index;
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      // First get scaled Euclidean distance
      r = 0.0;
      for (int k = 0; k < d; k++) {
        r += (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k);
      }
      // Then get K_00 (all methods use it)
      K_00 = exp(-r);
      if (grad1(i) == 0 and grad2(j) == 0) { // K_00
        covmat(i, j) = tau2 * K_00;
      } else if (grad1(i) == 0 and grad2(j) != 0) { // K_0d
        d_index = grad2(j)-1; // change indexing to start at zero
        covmat(i, j) = tau2 * K_00 * (2/theta(d_index)) * (x1(i, d_index) - x2(j, d_index));
      } else if (grad1(i) != 0 and grad2(j) == 0) { // K_d0
        d_index = grad1(i)-1; // change indexing to start at zero
        covmat(i, j) = tau2 * K_00 * (-2/theta(d_index)) * (x1(i, d_index) - x2(j, d_index));
      } else if (grad1(i) == grad2(j)) { // K_dd
        d_index = grad1(i)-1; // change indexing to start at zero
        covmat(i, j) = tau2 * K_00 * (2/theta(d_index)) * (1 - (2/theta(d_index))*(x1(i, d_index) - x2(j, d_index))*(x1(i, d_index) - x2(j, d_index)));
      } else { // K_df
        d_index = grad1(i)-1; // change indexing to start at zero
        f_index = grad2(j)-1; // change indexing to start at zero
        covmat(i, j) = tau2 * K_00 * (-4/(theta(d_index) * theta(f_index)))*(x1(i, d_index) - x2(j, d_index))*(x1(i, f_index) - x2(j, f_index));
      }
    }
  }
  // Add nugget to all diagonal elements (acts as jitter for numerical stability)
  if (n1 == n2) { 
    for (int i = 0; i < n1; i++) {
      covmat(i, i) += tau2*g;
    }
  }
  return covmat;
}

// [[Rcpp::export]]
arma::mat Matern(const arma::mat distmat, const double tau2, const double theta,
                    const double g, const double v) { 
  // distmat = matrix of SQUARED distances
  int n1 = distmat.n_rows;
  int n2 = distmat.n_cols;
  arma::mat covmat(n1, n2);
  double r;
  if (v == 0.5) { 
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = sqrt(distmat(i, j) / theta);
        covmat(i, j) = tau2 * exp(-r);
      }
    }
  } else if (v == 1.5) {
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = sqrt(3 * distmat(i, j) / theta);
        covmat(i, j) = tau2 * (1 + r) * exp(-r);
      }
    }
  } else if (v == 2.5) {
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = sqrt(5 * distmat(i, j) / theta);
        covmat(i, j) = tau2 * (1 + r + pow(r, 2) / 3) * exp(-r);
      }
    }
  } 
  if (n1 == n2) {
    for (int i = 0; i < n1; i++) 
      covmat(i, i) += tau2 * g;
  }
  return covmat;
}

// [[Rcpp::export]]
arma::mat MaternSep(const arma::mat x1, const arma::mat x2, const double tau2, 
                    const arma::vec theta, const double g, const double v) { 
  if (x1.n_cols != x2.n_cols) 
    stop("dimension of x1 and x2 do not match");
  if (theta.n_elem != x1.n_cols)
    stop("length of theta does not match dimension of x");
  int n1 = x1.n_rows;
  int n2 = x2.n_rows;
  int d = x1.n_cols;
  arma::mat covmat(n1, n2);
  double r;
  if (v == 0.5) { 
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = 0.0;
        for (int k = 0; k < d; k++)
          r += (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k);
        covmat(i, j) = tau2 * exp(-sqrt(r));
      }
    }
  } else if (v == 1.5) {
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = 0.0;
        for (int k = 0; k < d; k++)
          r += 3 * (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k);
        covmat(i, j) = tau2 * (1 + sqrt(r)) * exp(-sqrt(r));
      }
    }
  } else if (v == 2.5) {
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = 0.0;
        for (int k = 0; k < d; k++)
          r += 5 * (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k);
        covmat(i, j) = tau2 * (1 + sqrt(r) + r / 3) * exp(-sqrt(r));
      }
    }
  } 
  if (n1 == n2) {
    for (int i = 0; i < n1; i++) 
      covmat(i, i) += tau2 * g;
  }
  return covmat;
}
  
