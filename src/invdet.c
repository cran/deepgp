#include <R_ext/Utils.h>
#include <R.h>
#include <Rmath.h>
#include <assert.h>
#include <Rinternals.h>
#include "misc.h"

char uplo = 'U';

#define dposv dposv_
extern void dposv(char *, size_t *, size_t *, double*, size_t*, double*, size_t*, long*);

/* inv_det_R:
 * 
 * calculates the inverse and determinant of a matrix
 */

void inv_det_R(int *n, double *M, double *Mi, double *ldet)
{

  long info;
  unsigned int i;

  size_t n64;
  n64 = *n;
  dposv(&uplo,&n64,&n64,M,&n64,Mi,&n64,&info);

  *ldet = 0;
  for(i=0; i<*n; i++) *ldet += log(M[i + (*n) * i]);
  *ldet = 2 * (*ldet);
}

/*
 * distance_R:
 *
 * function for calculating the distance matrix between
 * the rows of X1 and X2, with output in D_out -- using
 * a built-in R interface
 */

void distance_R(double *X1_in, int *n1_in, double *X2_in,
                int *n2_in, int *m_in, double *D_out)
{
  double **X1, **X2, **D;

  /* make matrix bones */
  X1 = new_matrix_bones(X1_in, *n1_in, *m_in);
  X2 = new_matrix_bones(X2_in, *n2_in, *m_in);
  D = new_matrix_bones(D_out, *n1_in, *n2_in);

  unsigned int i,j,k;

  /* for each row of X1 and X2 */
  for(i=0; i<*n1_in; i++) {
    for(j=0; j<*n2_in; j++) {
      /* sum the squared entries */
      D[i][j] = 0.0;
      for(k=0; k<*m_in; k++) {
        D[i][j] += sq(X1[i][k] - X2[j][k]);
      }
    }
  }

  /* clean up */
  free(X1);
  free(X2);
  free(D);
}

/*
 * distance_symm_R:
 *
 * function for calculating the distance matrix between
 * the rows of X1 and itself, with output in the symmetric
 * D_out matrix -- using a built-in R interface
 */

void distance_symm_R(double *X_in, int *n_in, int *m_in, double *D_out)
{
  int n, m, i, j, k;
  double **X, **D;

  /* copy integers */
  n = *n_in;
  m = *m_in;

  /* make matrix bones */
  X = new_matrix_bones(X_in, n, m);
  D = new_matrix_bones(D_out, n, n);

  /* for each row of X and itself */
  for(i=0; i<n; i++) {
    D[i][i] = 0.0;
    for(j=i+1; j<n; j++) {
      D[i][j] = 0.0;
      for(k=0; k<m; k++)
        D[i][j] += sq(X[i][k] - X[j][k]);
      D[j][i] = D[i][j];
    }
  }

  /* clean up */
  free(X);
  free(D);
}



