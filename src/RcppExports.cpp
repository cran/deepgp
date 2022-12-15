#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rev_matrix
arma::mat rev_matrix(const arma::mat x);
RcppExport SEXP _deepgp_rev_matrix(SEXP xSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rev_matrix(x));
    return rcpp_result_gen;
    END_RCPP
}

// Exp2
arma::mat Exp2(const arma::mat distmat, const double tau2, const double theta,
                  const double g);
RcppExport SEXP _deepgp_Exp2(SEXP distmatSEXP, SEXP tau2SEXP, SEXP thetaSEXP,
                                SEXP gSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type distmat(distmatSEXP);
    Rcpp::traits::input_parameter< const double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type g(gSEXP);
    rcpp_result_gen = Rcpp::wrap(Exp2(distmat, tau2, theta, g));
    return rcpp_result_gen;
    END_RCPP
}

// Exp2Sep
arma::mat Exp2Sep(const arma::mat x1, const arma::mat x2, const double tau2, arma::vec theta,
               const double g);
RcppExport SEXP _deepgp_Exp2Sep(SEXP x1SEXP, SEXP x2SEXP, SEXP tau2SEXP, SEXP thetaSEXP,
                             SEXP gSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type g(gSEXP);
    rcpp_result_gen = Rcpp::wrap(Exp2Sep(x1, x2, tau2, theta, g));
    return rcpp_result_gen;
    END_RCPP
}

// Matern
arma::mat Matern(const arma::mat distmat, const double tau2, const double theta,
                    const double g, const double v);
RcppExport SEXP _deepgp_Matern(SEXP distmatSEXP, SEXP tau2SEXP, SEXP thetaSEXP,
                                  SEXP gSEXP, SEXP vSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type distmat(distmatSEXP);
    Rcpp::traits::input_parameter< const double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type g(gSEXP);
    Rcpp::traits::input_parameter< const double >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(Matern(distmat, tau2, theta, g, v));
    return rcpp_result_gen;
    END_RCPP
}

// MaternSep
arma::mat MaternSep(const arma::mat x1, const arma::mat x2, const double tau2, const arma::vec theta,
                 const double g, const double v);
RcppExport SEXP _deepgp_MaternSep(SEXP x1SEXP, SEXP x2SEXP, SEXP tau2SEXP, SEXP thetaSEXP,
                               SEXP gSEXP, SEXP vSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type g(gSEXP);
    Rcpp::traits::input_parameter< const double >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(MaternSep(x1, x2, tau2, theta, g, v));
    return rcpp_result_gen;
    END_RCPP
}

// U_entries
arma::mat U_entries(const int Ncores, const arma::mat& x, 
                    const arma::umat& revNNarray, const arma::mat& revCondOnLatent, 
                    const double tau2, const double theta, const double g, const double v);
RcppExport SEXP _deepgp_U_entries(SEXP NcoresSEXP, SEXP xSEXP, SEXP 
                                      revNNarraySEXP, SEXP revCondOnLatentSEXP, SEXP tau2SEXP,
                                      SEXP thetaSEXP, SEXP gSEXP, SEXP vSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type Ncores(NcoresSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type revNNarray(revNNarraySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type revCondOnLatent(revCondOnLatentSEXP);
    Rcpp::traits::input_parameter< const double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type g(gSEXP);
    Rcpp::traits::input_parameter< const double >::type v(vSEXP);    
    rcpp_result_gen = Rcpp::wrap(U_entries(Ncores, x, revNNarray, revCondOnLatent, tau2, 
                                           theta, g, v));
    return rcpp_result_gen;
    END_RCPP
}


// U_entries_sep
arma::mat U_entries_sep(const int Ncores, const arma::mat& x, 
                    const arma::umat& revNNarray, const arma::mat& revCondOnLatent, 
                    const double tau2, const arma::vec theta, const double g, const double v);
RcppExport SEXP _deepgp_U_entries_sep(SEXP NcoresSEXP, SEXP xSEXP, SEXP 
                                      revNNarraySEXP, SEXP revCondOnLatentSEXP, SEXP tau2SEXP,
                                      SEXP thetaSEXP, SEXP gSEXP, SEXP vSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type Ncores(NcoresSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type revNNarray(revNNarraySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type revCondOnLatent(revCondOnLatentSEXP);
    Rcpp::traits::input_parameter< const double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type g(gSEXP);
    Rcpp::traits::input_parameter< const double >::type v(vSEXP);    
    rcpp_result_gen = Rcpp::wrap(U_entries_sep(Ncores, x, revNNarray, revCondOnLatent, tau2, 
                                           theta, g, v));
    return rcpp_result_gen;
    END_RCPP
}

// check_omp
void check_omp();
RcppExport SEXP _deepgp_check_omp() {
    BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    check_omp();
    return R_NilValue;
    END_RCPP
}

extern "C" {
#include <stdlib.h>
#include <R_ext/Rdynload.h>
    extern void alcGP_R(void *, void *, void *, void *, void *, void *, void *,
                        void *, void *, void *, void *, void *, void *);
    extern void inv_det_R(void *, void *, void *, void *);
    extern void distance_R(void *, void *, void *, void *, void *, void *);
    extern void distance_symm_R(void *, void *, void *, void *);
    extern void Wij_R(void *, void *, void *, void *, void *, void *, void *,
                      void *, void *);
    
    static const R_CMethodDef CEntries[] = {
        {"alcGP_R",           (DL_FUNC) &alcGP_R,          13},
        {"inv_det_R",         (DL_FUNC) &inv_det_R,         4},
        {"distance_R",        (DL_FUNC) &distance_R,        6},
        {"distance_symm_R",   (DL_FUNC) &distance_symm_R,   4},
        {"Wij_R",             (DL_FUNC) &Wij_R,             9},
        {NULL, NULL, 0}
    };
}

static const R_CallMethodDef CallEntries[] = {
    {"_deepgp_rev_matrix",     (DL_FUNC) &_deepgp_rev_matrix,     1},
    {"_deepgp_Exp2",           (DL_FUNC) &_deepgp_Exp2,           4},
    {"_deepgp_Exp2Sep",        (DL_FUNC) &_deepgp_Exp2Sep,        5},
    {"_deepgp_Matern",         (DL_FUNC) &_deepgp_Matern,         5},
    {"_deepgp_MaternSep",      (DL_FUNC) &_deepgp_MaternSep,      6},
    {"_deepgp_U_entries",      (DL_FUNC) &_deepgp_U_entries,      8},
    {"_deepgp_U_entries_sep",  (DL_FUNC) &_deepgp_U_entries_sep,  8},
    {"_deepgp_check_omp",      (DL_FUNC) &_deepgp_check_omp,      0},
    {NULL, NULL, 0}
};

RcppExport void R_init_deepgp(DllInfo *dll) {
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
