#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ExpFun
arma::mat ExpFun(const arma::mat distmat, const arma::vec covparms);
RcppExport SEXP _deepgp_ExpFun(SEXP distmatSEXP, SEXP covparmsSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type distmat(distmatSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(ExpFun(distmat, covparms));
    return rcpp_result_gen;
    END_RCPP
}

// MaternFun
arma::mat MaternFun(const arma::mat distmat, const arma::vec covparms);
RcppExport SEXP _deepgp_MaternFun(SEXP distmatSEXP, SEXP covparmsSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type distmat(distmatSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(MaternFun(distmat, covparms));
    return rcpp_result_gen;
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
    {"_deepgp_ExpFun",     (DL_FUNC) &_deepgp_ExpFun,     2},
    {"_deepgp_MaternFun",  (DL_FUNC) &_deepgp_MaternFun,  2},
    {NULL, NULL, 0}
};

RcppExport void R_init_deepgp(DllInfo *dll) {
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
