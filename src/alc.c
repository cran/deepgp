#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>
#include "misc.h"

#ifndef CBLAS_ENUM_DEFINED_H
#define CBLAS_ENUM_DEFINED_H
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, AtlasConj=114};
#endif

#define dgemm dgemm_
extern void dgemm(char*, char*,  size_t*, size_t*, size_t*, double*,
                  double*, size_t*, double*, size_t*, double*, double*, size_t*);

#define dsymv dsymv_
extern void dsymv(char*, size_t*, double*, double*, size_t*,
                  double*, size_t*, double*, double*, size_t*);

#define ddot ddot_
extern double ddot(size_t*, double*, size_t*, double*, size_t*);

/*
 * linalg_dsymv:
 *
 * analog of dsymv in blas
 * assumed column major representation
 */

void linalg_dsymv(int n, double alpha, double **A, int lda, double *X,
                  int ldx, double beta, double *Y, int ldy)
{
  char uplo = 'U';
  size_t n64, lda64, ldy64, ldx64;
  n64 = n; lda64 = lda; ldx64 = ldx; ldy64 = ldy;
  dsymv(&uplo,&n64,&alpha,*A,&lda64,X,&ldx64,&beta,Y,&ldy64);
}

/*
 * linalg_dgemm:
 *
 * analog of dgemm in blas
 * assumed column major representation
 */

void linalg_dgemm(const enum CBLAS_TRANSPOSE TA, const enum CBLAS_TRANSPOSE TB,
                  int m, int n, int k, double alpha, double **A, int lda,
                  double **B, int ldb, double beta, double **C, int ldc)
{
  size_t m64, n64, k64, lda64, ldb64, ldc64;
  char ta, tb;
  m64 = m; n64 = n; k64 = k; lda64 = lda; ldb64 = ldb; ldc64 = ldc;
  if(TA == CblasTrans) ta = 'T'; else ta = 'N';
  if(TB == CblasTrans) tb = 'T'; else tb = 'N';
  dgemm(&ta,&tb,&m64,&n64,&k64,&alpha,*A,&lda64,*B,&ldb64,&beta,*C,&ldc64);
}

/*
 * linalg_ddot:
 *
 * analog of ddot in blas
 */

double linalg_ddot(int n, double *X, int ldx, double *Y, int ldy)
{
  double result;
  size_t n64,ldx64,ldy64;
  n64 = n; ldx64 = ldx; ldy64=ldy;
  result = ddot(&n64,X,&ldx64,Y,&ldy64);
  return result;
}

/*
 * covar:
 *
 * calculate the correlation (K) between X1 and X2 with
 * an isotropic power exponential correlation function
 * with range d
 */

void covar(const int col, double **X1, const int n1, double **X2,
           const int n2, double d, double **K)
{
  int i, j, k;
  
  /* calculate the covariance */
  for(i=0; i<n1; i++)
    for(j=0; j<n2; j++) {
      K[i][j] = 0.0;
      for(k=0; k<col; k++) K[i][j] += sq(X1[i][k] - X2[j][k]);
      K[i][j] = exp(0.0 - K[i][j]/d);
    }
}

/*
 * calc_g_mui_kxy:
 *
 * function for calculating the g vector, mui scalar, and
 * kxy vector for the IECI calculation; kx is length-n
 * utility space -- only implements isotropic case
 */

void calc_g_mui_kxy(const int col, double *x, double **X,
                    const int n, double **Ki, double **Xref,
                    const int m, double d, double g, double *gvec,
                    double *mui, double *kx, double *kxy)
{
  double mu_neg;
  int i;
  
  /* sanity check */
  if(m == 0) assert(!kxy && !Xref);
  
  covar(col, &x, 1, X, n, d, &kx);
  if(m > 0) covar(col, &x, 1, Xref, m, d, &kxy);
  
  linalg_dsymv(n,1.0,Ki,n,kx,1,0.0,gvec,1);
  
  *mui = 1.0 + g - linalg_ddot(n, kx, 1, gvec, 1);
  
  mu_neg = 0.0 - 1.0/(*mui);
  for(i=0; i<n; i++) gvec[i] *= mu_neg;
}

/*
 * calc_ktKikx:
 *
 * function for calculating the ktKikx vector used in the
 * IECI calculation -- writes over the KtKik input
 */

void calc_ktKikx(double *ktKik, const int m, double **k, const int n,
                 double *g, const double mui, double *kxy, double **Gmui,
                 double *ktGmui, double *ktKikx)
{
  int i;
  
  /* loop over all of the m candidates */
  for(i=0; i<m; i++) {
    
    ktKikx[i] = sq(linalg_ddot(n, k[i], 1, g, 1))*mui;
    
    ktKikx[i] += 2.0*linalg_ddot(n, k[i], 1, g, 1)*kxy[i];
    
    ktKikx[i] += sq(kxy[i])/mui;
  }
}

/*
 * calc_alc:
 *
 * function that iterates over the m Xref locations, and the
 * stats calculated by previous calc_* function in order to
 * calculate the reduction in variance
 */

double calc_alc(const int m, double *ktKik, double *s2p, const double phi,
                double *badj, const double tdf, double *w)
{
  int i;
  double alc;
  
  alc = 0.0;
  for(i=0; i<m; i++) {
    alc += phi*ktKik[i];
  }
  
  return (alc/m);
}

/*
 * alcGP:
 *
 * return s2' component of the ALC calculation of the
 * expected reduction in variance calculation at locations
 * Xcand averaging over reference locations Xref:
 * ds2 = s2 - s2', where the s2s are at Xref and the
 * s2' incorporates Xcand, and everything is averaged
 * over Xref.
 */

void alcGP(int col, double **X, double **Ki, double d, double g,
           unsigned int n, double phi,
           unsigned int ncand, double **Xcand, unsigned int nref,
           double **Xref,  int verb, double *alc)
{
  int i;
  double **k;
  double *kx, *kxy, *gvec, *ktKikx;
  double mui, df;
  double s2p[2] = {0, 0};
  
  
  /* degrees of freedom */
  df = (double) n;
  
  /* allocate g, kxy, and ktKikx vectors */
  gvec = new_vector(n);
  kxy = new_vector(nref);
  kx = new_vector(n);
  ktKikx = new_vector(nref);
  
  k = new_matrix(nref, n);
  covar(col, Xref, nref, X, n, d, k);    //// CHANGED ncand to col
  
  /* calculate the ALC for each candidate */
  for(i=0; i<ncand; i++) {
    
    /* calculate the g vector, mui, and kxy */
    calc_g_mui_kxy(col, Xcand[i], X, n, Ki, Xref, nref, d,
                   g, gvec, &mui, kx, kxy);
    
    /* use g, mu, and kxy to calculate ktKik.x */
    calc_ktKikx(NULL, nref, k, n, gvec, mui, kxy, NULL, NULL, ktKikx);
    
    /* calculate the ALC */
    alc[i] = calc_alc(nref, ktKikx, s2p, phi, NULL, df, NULL);
  }
  
  /* clean up */
  free(ktKikx);
  free(gvec);
  free(kx);
  free(kxy);
  delete_matrix(k);
}

/*
 * alcGP_R:
 *
 * R interface to C-side function that returns the
 * s2' component of the ALC calculation of the
 * expected reduction in variance calculation given
 * the stored GP parameterization at locations Xcand
 * averaging over reference locations Xref:
 * ds2 = s2 - s2', where the s2s are at Xref and the
 * s2' incorporates Xcand, and everything is averaged
 * over Xref.
 */

void alcGP_R(double *X_in, int *n_in, int *col_in, double *Ki_in, double *d_in,
             double *g_in, int *ncand_in, double *Xcand_in, int *nref_in,
             double *Xref_in, double *phi_in, int *verb_in, double *alc_out)
{
  double **Xref, **Xcand, **X, **Ki;
  
  /* make matrix bones */
  Xref = new_matrix_bones(Xref_in, *nref_in, *col_in);
  Xcand = new_matrix_bones(Xcand_in, *ncand_in, *col_in);
  Ki = new_matrix_bones(Ki_in, *n_in, *n_in);
  X = new_matrix_bones(X_in, *n_in, *col_in);
  
  alcGP(*col_in, X, Ki, *d_in, *g_in, *n_in, *phi_in,
        *ncand_in, Xcand, *nref_in, Xref, *verb_in, alc_out);
  
  /* clean up */
  free(Xref);
  free(Xcand);
  free(Ki);
  free(X);
}

