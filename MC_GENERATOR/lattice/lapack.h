#ifndef _LAPACK_H_
#define _LAPACK_H_

#include <complex.h>
#include <geometry.h>


/* *************************************************************************** */
/* ****  Computes eigenvalues of general complex NxN matrix (C-ordered)  ***** */
/* ****  NOTE: original matrix A is destroyed on exit                    ***** */
/* *************************************************************************** */
extern void lapack_cmp_Eigenvalues( double complex *A, uint N , double complex *W );



/* *************************************************************************** */
/* Performs SVD decomposition  A = U W V of arbitrary complex MxN matrix A.
      A    - C-ordered (A[i][j] = A[i*N + j] ) matrix to consider;
      M, N - number of rows and columns of A correspondingly;
      W    - N-dim double array, on output contains singular (non-negative!)
             values of A;
      V    - unitary NxN matrix V exactly as in SVD formula above.
   On output A is replaced by MxN matrix U */
/* *************************************************************************** */
extern void lapack_cmp_svd( double complex *A, uint M, uint N, double *W, double complex *V );


/* *************************************************************************** */
/* ****  Computes eigenvalues and/or eigenvectors of              ************ */
/* ****  complex Hermitian NxN matrix A (C-ordered)               ************ */
/* *************************************************************************** */
extern void lapack_cmp_Eigensystem( double complex *A, uint N, double *D  );
extern void lapack_cmp_Eigensystem_noeigenvectors( double complex *A, uint N, double *D  );


/* *************************************************************************** */
/* Computes an LU factorization of a general M-by-N matrix A  using partial pivoting
   with row interchanges.  Factorization has the form
          A = P * L * U
   where P is a permutation matrix, L is lower triangular with unit diagonal  elements
   (lower  trapezoidal if m > n), and U is upper triangular (upper trapezoidal if m < n).*/
/* *************************************************************************** */
extern void lapack_cmp_LU( double complex *A, uint M, uint N, uint *perm );
extern void lapack_dbl_LU( double         *A, uint M, uint N, uint *perm );


/* *************************************************************************** */
/* ***** Determinants via above LU decomposition ***************************** */
/* *************************************************************************** */
extern double complex lapack_cmp_determinant( double complex *A, uint N );
extern double         lapack_dbl_determinant( double         *A, uint N );



/* *************************************************************************** */
/* ******* Variant of QR factorization, see man page for zungqr ************** */
/* *************************************************************************** */
extern void lapack_zungqr( double complex *A, uint M, uint N );


/* *************************************************************************** */
/* ** Left/Right eigenvectors/eigenvalues of general complex square matrix *** */
/* ** Original matrix A is destroyed on exit. For definition of left/right *** */
/* ** eigenvectors see manpage for zgeev     ********************************* */
/* *************************************************************************** */
extern void lapack_zgeevl( uint N,  double complex *A, double complex *W, double complex *VL );
extern void lapack_zgeevr( uint N,  double complex *A, double complex *W, double complex *VR );
extern void lapack_zge( uint N,  double complex *A, double complex *W );


/* *************************************************************************** */
/* ****  Computes eigenvalues and eigenvectors of    ************************* */
/* ****  real symmetric NxN matrix A (C-ordered)     ************************* */
/* ****  Original A is either destroyed or replaced  ************************* */
/* ****  with eigenvectors (in columns)              ************************* */
/* *************************************************************************** */
extern void lapack_dbl_Eigensystem( double *A, uint N, double *D  );
extern void lapack_dbl_Eigensystem_noeigenvectors( double *A, uint N, double *D  );
#endif
