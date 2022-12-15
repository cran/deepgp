
#include <iostream>
#include <math.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat Exp2(arma::mat distmat, const double tau2, const double theta,
                  const double g);

arma::mat Exp2Sep(arma::mat x1, arma::mat x2, const double tau2, 
                     arma::vec theta, const double g);

arma::mat Matern(arma::mat distmat, const double tau2, const double theta,
                    const double g, const double v);

arma::mat MaternSep(arma::mat x1, arma::mat x2, const double tau2, arma::vec theta,
                 const double g, const double v);

