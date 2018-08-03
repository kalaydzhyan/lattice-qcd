#ifndef _SVD_H_
#define _SVD_H_

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <complex.h>

#include <geometry.h>

#define SVD_MAX_ITER  30

/* Performs SVD decomposition 
       A = U W V
   of arbitrary complex MxN matrix A. It was promised that routine works always ;)
   In practice - it is NOT: for M < N expect problems (I'm too lazy to correct this).
      A    - C-ordered (A[i][j] = A[i*N + j] ) matrix to consider;
      M, N - number of rows and columns of A correspondingly;
      W    - N-dim double array, on output contains singular (non-negative!) values
             of A;
      V    - unitary NxN matrix V exactly as in SVD formula above.
   On output A is replaced by MxN matrix U
*/
extern void cmp_svd( double complex *A, uint M, uint N, double *W, double complex *V );


/*
  Householder reduction of hermitian matrix A[N][N] (first arg - &A[0][0]) to
  tridiagonal form A'. On output A is replaced by unitary matrix Q effecting
  the transformation. The relation between matrices:
                     A' = Q^+ * A * Q
  D returns diagonal elements and E returns super-diagonal ones with E[0] = 0.
*/
extern void cmp_householder( double complex *A, uint N, double *D, double *E );
extern void cmp_householder_noeigenvectors( double complex *A, uint N, double *D, double *E );


/*
  QL with implicit shifts, allows to diagonalize tridiagonal real symmetric NxN input
  matrix A, stored in C-order. Finds eigenvalues and eigenvectors.
    D (dim D = N) -- diagonal elements of A, on output contains eigenvalues.
    E (dim E = N) -- sub/super-diagonal elements of A with E[0] - undefined,
                     on output is destroyed.
  Matrix of eigenvectors Q is returned via A with:
      diag( D ) = A^T * [original matrix] * A
*/
extern void cmp_QLShifts( double complex *A, uint N, double *D, double *E );
extern void cmp_QLShifts_noeigenvectors( double complex *A, uint N, double *D, double *E );

/*
  Wrapper for two above routines - finds eigensystem of complex hermitian matrix A
  stored in C-order. On output A is replaced with eigenvectors, D containes eigenvalues.
  Relation between matrices on output:
        diag[ D ] = A^+ * [original matrix] * A
*/
extern void cmp_Eigensystem( double complex *A, uint N, double *D );
extern void cmp_Eigensystem_noeigenvectors( double complex *A, uint N, double *D );

extern double complex cmp_determinant( double complex *A, uint N );
extern double         dbl_determinant( double         *A, uint N );
#endif
