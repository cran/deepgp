#include <assert.h>
#include <math.h>
#include <stdlib.h>

/*
 * new_vector:
 *
 * allocates a new double array of size n1
 */

double* new_vector(unsigned int n)
{
  double *v;
  if(n == 0) return NULL;
  v = (double*) malloc(sizeof(double) * n);
  return v;
}

/*
 * new_matrix:
 *
 * create a new n1 x n2 matrix which is allocated like
 * and n1*n2 array, but can be referenced as a 2-d array
 */

double** new_matrix(unsigned int n1, unsigned int n2)
{
  int i;
  double **m;
  
  if(n1 == 0 || n2 == 0) return NULL;
  m = (double**) malloc(sizeof(double*) * n1);
  m[0] = (double*) malloc(sizeof(double) * (n1*n2));
  for(i=1; i<n1; i++) m[i] = m[i-1] + n2;
  
  return m;
}

/*
 * new_matrix_bones:
 *
 * create a double ** Matrix from a double * vector
 * should be freed with the free command, rather than
 * delete_matrix
 */

double ** new_matrix_bones(double *v, unsigned int n1, unsigned int n2)
{
  double **M;
  int i;
  M = (double **) malloc(sizeof(double*) * n1);
  M[0] = v;
  for(i=1; i<n1; i++) M[i] = M[i-1] + n2;
  return(M);
}

/*
 * delete_matrix:
 *
 * delete a matrix allocated as above
 */

void delete_matrix(double** m)
{
  if(m == NULL) return;
  assert(*m);
  free(*m);
  assert(m);
  free(m);
}

/*
 * sq:
 *
 * calculate the square of x
 */

double sq(double x)
{
  return x*x;
}


