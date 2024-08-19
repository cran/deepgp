#include <Rcpp.h>
using namespace Rcpp;

/*
 * This file contains helper functions for the "fo_approx" R function
 */

int less_than_index(const NumericVector xg, const double x) {
  // Returns the index of xg that is just ABOVE x location
  // A value of 1 indicates that x is below all values of xg
  // A value of xg.length + 1 indicates that x is above all values of xg
  // Indexing starts at 1 (there will not be any zeros)
  int ng = xg.length();
  int index = 0;
  for (int i = 1; i < (ng + 1); i++) {
    if (xg(i - 1) <= x) 
      index = i;
  }
  return index + 1; // to change indexing from 0 to 1
}

// [[Rcpp::export]]
NumericMatrix fo_approx_init(const NumericMatrix xg, const NumericMatrix x) {
  // xg: grid locations
  // x: observed locations
  // Returns matrix whose entries contain the index of xg which
  // is just ABOVE that observation in x
  // Indexing starts at 1
  // Note: xg must be pre-sorted from lowest to highest
  int d = x.ncol();
  int n = x.nrow();
  NumericMatrix lower_ind(n, d);
  for (int j = 0; j < d; j++) {
    for (int i = 0; i < n; i++) {
      lower_ind(i, j) = less_than_index(xg(_, j), x(i, j));
    }
  }
  return lower_ind;
}
