
#define BOOST_DISABLE_ASSERTS
#define ARMA_DONT_PRINT_ERRORS

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(BH)]]

#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <boost/math/special_functions/gamma.hpp>

using namespace Rcpp;
using namespace arma;
using namespace std;

double sqdist(rowvec l1, rowvec l2) { 
  double ssq = 0.0;
  for(arma::uword k = 0; k < l1.size(); ++k) {
    ssq += (l1[k] - l2[k])*(l1[k] - l2[k]);
  }
  return ssq;
}

arma::mat calc_sqdist(arma::mat x) {
  arma::uword outrows = x.n_rows;
  arma::uword outcols = x.n_rows;
  arma::mat out(outrows, outcols);
  for (arma::uword arow = 0 ; arow < outrows ; arow++) {
    for (arma::uword acol = 0 ; acol < outcols ; acol++) {
      out(arow, acol) = sqdist(x.row(arow), x.row(acol));
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat ExpFun(arma::mat distmat, arma::vec covparms) { 
  // distmat = matrix of SQUARED distances
  // covparms = c(tau2, theta, g, v = NULL)
  int d1 = distmat.n_rows;
  int d2 = distmat.n_cols;
  int j1;
  int j2;
  arma::mat covmat(d1, d2);
  double scaledist;
  for (j1 = 0; j1 < d1; j1++) {
    for (j2 = 0; j2 < d2; j2++) {
      if (distmat(j1, j2) == 0) {
        covmat(j1, j2) = covparms(0)*(1 + covparms(2));
      } else {
        scaledist = distmat(j1, j2)/covparms(1);
        covmat(j1, j2) = covparms(0)*exp(-scaledist);
      }
    }
  }
  return covmat;
}

// [[Rcpp::export]]
arma::mat MaternFun(arma::mat distmat, arma::vec covparms) { 
  // distmat = matrix of SQUARED distances
  // covparms = c(tau2, theta, g, v)
  int d1 = distmat.n_rows;
  int d2 = distmat.n_cols;
  int j1;
  int j2;
  arma::mat covmat(d1, d2);
  double scaledist;
  double alpha = sqrt(covparms(1));
  if (covparms(3) == 0.5) { 
    for (j1 = 0; j1 < d1; j1++) {
      for (j2 = 0; j2 < d2; j2++) {
        if (distmat(j1, j2) == 0) {
          covmat(j1, j2) = covparms(0)*(1 + covparms(2));
        } else {
          scaledist = sqrt(distmat(j1, j2))/alpha;
          covmat(j1, j2) = covparms(0)*exp(-scaledist);
        }
      }
    }
  } else if(covparms(3) == 1.5) {
    double normcon = (pow(alpha, -3)/sqrt(M_PI))*(boost::math::tgamma(2)/boost::math::tgamma(1.5));
    for (j1 = 0; j1 < d1; j1++) {
      for (j2 = 0; j2 < d2; j2++) {
        if (distmat(j1, j2) == 0) {
          covmat(j1,j2) = covparms(0)*(1 + covparms(2));
        } else {
          scaledist = sqrt(distmat(j1, j2))/alpha;
          covmat(j1, j2) = normcon*covparms(0)*0.5*M_PI*pow(alpha,3)*exp(-scaledist)*(1+scaledist);
        }
      }
    }
  } else if(covparms(3) == 2.5) {
    double normcon = (pow(alpha, -5)/sqrt(M_PI))*(boost::math::tgamma(3)/boost::math::tgamma(2.5));
    for (j1 = 0; j1 < d1; j1++) {
      for (j2 = 0; j2 < d2; j2++) {
        if (distmat(j1, j2) == 0) {
          covmat(j1,j2) = covparms(0)*(1 + covparms(2));
        } else {
          scaledist = sqrt(distmat(j1, j2))/alpha;
          covmat(j1, j2) = normcon*covparms(0)*0.125*M_PI*pow(alpha,5)*exp(-scaledist)*(3+3*scaledist+scaledist*scaledist);
        }
      }
    }
  } 
  return covmat;
}

