#include <Rmath.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "general.h"

double pnorm(double x, double mu, double sigma, int lower_tail,
             int give_log);

/*
 * er_function:
 * 
 * calculates Gaussian error function
 */ 

double er_function(double x) 
{
  double m;

  m = 2. * pnorm(x * sqrt(2.), 0.0,1.0,1,0) - 1.;
  return m;
}

/*
 * Wij:
 * 
 * calculates Wij matrix
 */

void Wij(double ** W, double ** X1, int n1, double ** X2, int n2, int col,
         double theta, double a, double b) 
{
  double comp1, comp2, y1, y2;
  int i, j, k;

  for (i = 0; i<n1; i++) {
    for (j=0; j<n2; j++) {
      W[i][j] = 1.0; /* initialize product */
      for (k=0; k<col; k++) {
        y1 = X1[i][k];
        y2 = X2[j][k];
        comp1 = (2. * b - y1 - y2)/sqrt(2. * theta);
        comp2 = (2. * a - y1 - y2)/sqrt(2. * theta);
        W[i][j] *= sqrt(M_PI * theta/8.) * exp(-(y1 - y2)*(y1 - y2)/(2. * theta)) * (er_function(comp1) - er_function(comp2));
      }
    }
  }
}

/*
 * Wij_R:
 * 
 * R interface to the C-side function that returns the
 * Wij matrix
 */

void Wij_R(double *X1_in, int *n1_in, double *X2_in, int *n2_in, int *col_in,
            double *theta_in, double *a_in, double *b_in, double *W_out)
{
  double **X1, **X2, **W;

  /* new matrix bones */
  X1 = new_matrix_bones(X1_in, *n1_in, *col_in);
  X2 = new_matrix_bones(X2_in, *n2_in, *col_in);
  W = new_matrix_bones(W_out, *n1_in, *n2_in);

  /* call C function */
  Wij(W, X1, *n1_in, X2, *n2_in, *col_in, *theta_in, *a_in, *b_in);

  free(X1);
  free(X2);
}
