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

// Exp2Fun
arma::mat Exp2Fun(const arma::mat distmat, const arma::vec covparms);
RcppExport SEXP _deepgp_Exp2Fun(SEXP distmatSEXP, SEXP covparmsSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type distmat(distmatSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(Exp2Fun(distmat, covparms));
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

// U_entries
arma::mat U_entries(const int Ncores, const arma::uword n, const arma::mat& locs, const arma::umat& revNNarray, const arma::mat& revCondOnLatent, const arma::vec covparms);
RcppExport SEXP _deepgp_U_entries(SEXP NcoresSEXP, SEXP nSEXP, SEXP locsSEXP, SEXP revNNarraySEXP, SEXP revCondOnLatentSEXP, SEXP covparmsSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type Ncores(NcoresSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type locs(locsSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type revNNarray(revNNarraySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type revCondOnLatent(revCondOnLatentSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(U_entries(Ncores, n, locs, revNNarray, revCondOnLatent, covparms));
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
    {"_deepgp_rev_matrix",  (DL_FUNC) &_deepgp_rev_matrix,  1},
    {"_deepgp_Exp2Fun",     (DL_FUNC) &_deepgp_Exp2Fun,     2},
    {"_deepgp_MaternFun",   (DL_FUNC) &_deepgp_MaternFun,   2},
    {"_deepgp_U_entries",   (DL_FUNC) &_deepgp_U_entries,   6},
    {"_deepgp_check_omp",   (DL_FUNC) &_deepgp_check_omp,   0},
    {NULL, NULL, 0}
};

RcppExport void R_init_deepgp(DllInfo *dll) {
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
