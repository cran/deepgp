#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// forward_solve
arma::vec forward_solve_raw(const arma::mat& U, const arma::vec z,
                            const arma::umat& NNarray);
RcppExport SEXP _deepgp_forward_solve_raw(SEXP USEXP, SEXP zSEXP,
                                      SEXP NNarraySEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type NNarray(NNarraySEXP);
    rcpp_result_gen = Rcpp::wrap(forward_solve_raw(U, z, NNarray));
    return rcpp_result_gen;
    END_RCPP
}

// diag_quad_mat
arma::vec diag_quad_mat(const arma::mat A, const arma::mat B);
RcppExport SEXP _deepgp_diag_quad_mat(SEXP ASEXP, SEXP BSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(diag_quad_mat(A, B));
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

// Exp2Grad
arma::mat Exp2Grad(const arma::mat x1, const arma::mat x2, const arma::vec grad1, 
                   const arma::vec grad2, const double tau2, 
                   const double theta, const double g);
RcppExport SEXP _deepgp_Exp2Grad(SEXP x1SEXP, SEXP x2SEXP, SEXP grad1SEXP,
                                 SEXP grad2SEXP, SEXP tau2SEXP, SEXP thetaSEXP,
                                 SEXP gSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< const arma::mat >::type x1(x1SEXP);
  Rcpp::traits::input_parameter< const arma::mat >::type x2(x2SEXP);
  Rcpp::traits::input_parameter< const arma::vec >::type grad1(grad1SEXP);
  Rcpp::traits::input_parameter< const arma::vec >::type grad2(grad2SEXP);
  Rcpp::traits::input_parameter< const double >::type tau2(tau2SEXP);
  Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
  Rcpp::traits::input_parameter< const double >::type g(gSEXP);
  rcpp_result_gen = Rcpp::wrap(Exp2Grad(x1, x2, grad1, grad2, tau2, theta, g));
  return rcpp_result_gen;
  END_RCPP
}

// Exp2Sep
arma::mat Exp2Sep(const arma::mat x1, const arma::mat x2, const double tau2, 
                  const arma::vec theta, const double g);
RcppExport SEXP _deepgp_Exp2Sep(SEXP x1SEXP, SEXP x2SEXP, SEXP tau2SEXP, 
                                SEXP thetaSEXP, SEXP gSEXP) {
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

// Exp2SepGrad
arma::mat Exp2SepGrad(const arma::mat x1, const arma::mat x2, const arma::vec grad1,
                      const arma::vec grad2, const double tau2, const arma::vec theta,
                      const double g);
RcppExport SEXP _deepgp_Exp2SepGrad(SEXP x1SEXP, SEXP x2SEXP, SEXP grad1SEXP,
                                    SEXP grad2SEXP, SEXP tau2SEXP, SEXP thetaSEXP,
                                    SEXP gSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< const arma::mat >::type x1(x1SEXP);
  Rcpp::traits::input_parameter< const arma::mat >::type x2(x2SEXP);
  Rcpp::traits::input_parameter< const arma::vec >::type grad1(grad1SEXP);
  Rcpp::traits::input_parameter< const arma::vec >::type grad2(grad2SEXP);
  Rcpp::traits::input_parameter< const double >::type tau2(tau2SEXP);
  Rcpp::traits::input_parameter< const arma::vec >::type theta(thetaSEXP);
  Rcpp::traits::input_parameter< const double >::type g(gSEXP);
  rcpp_result_gen = Rcpp::wrap(Exp2SepGrad(x1, x2, grad1, grad2, tau2, theta, g));
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
arma::mat MaternSep(const arma::mat x1, const arma::mat x2, const double tau2, 
                    const arma::vec theta, const double g, const double v);
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
arma::mat U_entries(const int cores, const arma::mat& x, const arma::umat& NNarray,
                    const double tau2, const double theta, const double g, const double v);
RcppExport SEXP _deepgp_U_entries(SEXP coresSEXP, SEXP xSEXP, SEXP NNarraySEXP, 
                                  SEXP tau2SEXP, SEXP thetaSEXP, SEXP gSEXP, SEXP vSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type cores(coresSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type NNarray(NNarraySEXP);
    Rcpp::traits::input_parameter< const double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type g(gSEXP);
    Rcpp::traits::input_parameter< const double >::type v(vSEXP);    
    rcpp_result_gen = Rcpp::wrap(U_entries(cores, x, NNarray, tau2, theta, g, v));
    return rcpp_result_gen;
    END_RCPP
}

// U_entries_grad
arma::mat U_entries_grad(const int cores, const arma::mat& x, const arma::umat& NNarray,
                         const arma::vec grad, const double tau2, 
                         const double theta, const double g, const double v);
RcppExport SEXP _deepgp_U_entries_grad(SEXP coresSEXP, SEXP xSEXP, SEXP NNarraySEXP, 
                                      SEXP gradSEXP, SEXP tau2SEXP, 
                                      SEXP thetaSEXP, SEXP gSEXP, SEXP vSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< const int >::type cores(coresSEXP);
  Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
  Rcpp::traits::input_parameter< const arma::umat& >::type NNarray(NNarraySEXP);
  Rcpp::traits::input_parameter< const arma::vec >::type grad(gradSEXP);
  Rcpp::traits::input_parameter< const double >::type tau2(tau2SEXP);
  Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
  Rcpp::traits::input_parameter< const double >::type g(gSEXP);
  Rcpp::traits::input_parameter< const double >::type v(vSEXP);    
  rcpp_result_gen = Rcpp::wrap(U_entries_grad(cores, x, NNarray, grad, tau2, theta, g, v));
  return rcpp_result_gen;
  END_RCPP
}

// U_entries_sep
arma::mat U_entries_sep(const int cores, const arma::mat& x, 
                        const arma::umat& NNarray, const double tau2, 
                        const arma::vec theta, const double g, const double v);
RcppExport SEXP _deepgp_U_entries_sep(SEXP coresSEXP, SEXP xSEXP, SEXP NNarraySEXP, 
                                      SEXP tau2SEXP, SEXP thetaSEXP, SEXP gSEXP, SEXP vSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type cores(coresSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type NNarray(NNarraySEXP);
    Rcpp::traits::input_parameter< const double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type g(gSEXP);
    Rcpp::traits::input_parameter< const double >::type v(vSEXP);    
    rcpp_result_gen = Rcpp::wrap(U_entries_sep(cores, x, NNarray, tau2, 
                                               theta, g, v));
    return rcpp_result_gen;
    END_RCPP
}

// U_entries_sep_grad
arma::mat U_entries_sep_grad(const int cores, const arma::mat& x, 
                             const arma::umat& NNarray, const arma::vec grad,
                             const double tau2, const arma::vec theta, 
                             const double g, const double v);
RcppExport SEXP _deepgp_U_entries_sep_grad(SEXP coresSEXP, SEXP xSEXP, SEXP NNarraySEXP, 
                                           SEXP gradSEXP, SEXP tau2SEXP, 
                                           SEXP thetaSEXP, SEXP gSEXP, SEXP vSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< const int >::type cores(coresSEXP);
  Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
  Rcpp::traits::input_parameter< const arma::umat& >::type NNarray(NNarraySEXP);
  Rcpp::traits::input_parameter< const arma::vec >::type grad(gradSEXP);
  Rcpp::traits::input_parameter< const double >::type tau2(tau2SEXP);
  Rcpp::traits::input_parameter< const arma::vec >::type theta(thetaSEXP);
  Rcpp::traits::input_parameter< const double >::type g(gSEXP);
  Rcpp::traits::input_parameter< const double >::type v(vSEXP);    
  rcpp_result_gen = Rcpp::wrap(U_entries_sep_grad(cores, x, NNarray, grad, tau2, 
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

// row_col_pointers
arma::mat row_col_pointers(const arma::umat& NNarray);
RcppExport SEXP _deepgp_row_col_pointers(SEXP NNarraySEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat& >::type NNarray(NNarraySEXP);
    rcpp_result_gen = Rcpp::wrap(row_col_pointers(NNarray));
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
    {"_deepgp_forward_solve_raw",  (DL_FUNC) &_deepgp_forward_solve_raw,  3},
    {"_deepgp_diag_quad_mat",      (DL_FUNC) &_deepgp_diag_quad_mat,      2},
    {"_deepgp_Exp2",               (DL_FUNC) &_deepgp_Exp2,               4},
    {"_deepgp_Exp2Grad",           (DL_FUNC) &_deepgp_Exp2Grad,           7},
    {"_deepgp_Exp2Sep",            (DL_FUNC) &_deepgp_Exp2Sep,            5},
    {"_deepgp_Exp2SepGrad",        (DL_FUNC) &_deepgp_Exp2SepGrad,        7},
    {"_deepgp_Matern",             (DL_FUNC) &_deepgp_Matern,             5},
    {"_deepgp_MaternSep",          (DL_FUNC) &_deepgp_MaternSep,          6},
    {"_deepgp_U_entries",          (DL_FUNC) &_deepgp_U_entries,          7},
    {"_deepgp_U_entries_grad",     (DL_FUNC) &_deepgp_U_entries_grad,     8},
    {"_deepgp_U_entries_sep",      (DL_FUNC) &_deepgp_U_entries_sep,      7},
    {"_deepgp_U_entries_sep_grad", (DL_FUNC) &_deepgp_U_entries_sep_grad, 8},
    {"_deepgp_check_omp",          (DL_FUNC) &_deepgp_check_omp,          0},
    {"_deepgp_row_col_pointers",   (DL_FUNC) &_deepgp_row_col_pointers,   1},
    {NULL, NULL, 0}
};

RcppExport void R_init_deepgp(DllInfo *dll) {
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
