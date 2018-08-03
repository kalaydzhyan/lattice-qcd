#ifndef _DIAGONALIZE_H_
#define _DIAGONALIZE_H_

#include <stdlib.h>
#include <math.h>

/*
  Performs diagonalization of given real symmetric NxN matrix A, finds its eigenvectors and
  eigenvalues. A is stored as: element (i,j) = *( *(A+i) + j). Eigenvalues are returned via
  vector D of length N. On output A is replaced with matrix Q of eigenvectors. Relation betweeen
  matrices:
                 diag( D ) = Q^T * A * Q

  Return values:  0 for success and -1 for errors (iteration count is too large,
  memory allocation failed).
*/
extern int RealSymmetricDiagonalize( double **A, double *D, unsigned int N );

#define MATRIX_ELEMENT(i,j,N)   ((i) + (j) * (N))

extern int Real_Symmetric_Diagonalize( double *A, double *D, unsigned int N );

#endif
