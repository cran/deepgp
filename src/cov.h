
#include <iostream>
#include <math.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

NumericVector diag_quad_mat(NumericMatrix A, NumericMatrix B);
  
arma::mat Exp2(arma::mat distmat, const double tau2, const double theta,
                  const double g);

arma::mat Exp2Grad(const arma::mat x1, const arma::mat x2, const arma::vec grad1, 
                   const arma::vec grad2, const double tau2, 
                   const double theta, const double g);

arma::mat Exp2Sep(arma::mat x1, arma::mat x2, const double tau2, 
                     arma::vec theta, const double g);

arma::mat Exp2SepGrad(const arma::mat x1, const arma::mat x2, const arma::vec grad1,
                      const arma::vec grad2, const double tau2, 
                      const arma::vec theta, const double g);

arma::mat Matern(arma::mat distmat, const double tau2, const double theta,
                    const double g, const double v);

arma::mat MaternSep(arma::mat x1, arma::mat x2, const double tau2, arma::vec theta,
                 const double g, const double v);

