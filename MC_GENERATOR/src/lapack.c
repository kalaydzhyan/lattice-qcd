#include <lapack.h>

#define UNDERSC(a)  a ## _
#define UNDERSC2(a) a ## _                                                                                         

#define ELEM(i,j,N)   ((i) * (N) + (j))

/* *************************************************************************** */
/* *****  FORTRAN interface declarations               *********************** */
/* *****  See corresponding man pages [lapack distro]  *********************** */
/* *************************************************************************** */
extern void UNDERSC(zgehrd)( int *N, int *ILO, int *IHI,
			     double complex *A, int *LDA, double complex *TAU,
			     double complex *WORK, int *LWORK, int *INFO );

extern void UNDERSC(zhseqr)( char *JOB, char *COMPZ, int *N, int *ILO, int *IHI,
			     double complex *H, int *LDH,
			     double complex *W,
			     double complex *Z, int *LDZ,
			     double complex *WORK, int *LWORK, int *INFO );

extern void UNDERSC(zgesvd)( char *JOBU, char *JOBVT, int *M, int *N,
			     double complex *A, int *LDA,
			     double *S,
			     double complex *U, int *LDU,
			     double complex *VT, int *LDVT,
			     double complex *WORK, int *LWORK, double *RWORK,
			     int *INFO );

extern void UNDERSC(zheev)(  char *JOBZ, char *UPLO, int *N,
			     double complex *A, int *LDA, double *W,
			     double complex *WORK, int *LWORK, double *RWORK,
			     int *INFO );

extern void UNDERSC(zgetrf)( int *M, int *N, double complex *A, int *LDA, int *IPIV,
			     int *INFO);

extern void UNDERSC(dgetrf)( int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);

extern void UNDERSC(zgeqrf)( int *M, int *N, double complex *A, int *LDA,
			     double complex *TAU, double complex *WORK, int *LWORK, int *INFO );

extern void UNDERSC(zungqr)( int *M, int *N, int *K, double complex *A, int *LDA,
			     double complex *TAU, double complex *WORK, int *LWORK, int *INFO );

extern void UNDERSC(zgeev)( char *JOBVL, char *JOBVR, int *N,
			    double complex *A, int *LDA,
			    double complex *W,
			    double complex *VL, int *LDVL,
			    double complex *VR, int *LDVR,
			    double complex *WORK, int *LWORK, double *RWORK,
			    int *INFO );
extern void UNDERSC(dsyev)( char *JOBZ, char *UPLO, int *N,
			    double *A, int *LDA,
			    double *W, double *WORK, int *LWORK, int *INFO );
/* *************************************************************************** */
/* *****  FORTRAN interface routines  **************************************** */
/* *************************************************************************** */
static void zgehrd( int N, int ILO, int IHI,
		    double complex *A, int LDA, double complex *TAU,
		    double complex *WORK, int LWORK, int *INFO ){
  UNDERSC(zgehrd)( &N, &ILO, &IHI, A, &LDA, TAU, WORK, &LWORK, INFO ); 
};

static void zhseqr( char JOB, char COMPZ, int N, int ILO, int IHI,
		    double complex *H, int LDH,
		    double complex *W,
		    double complex *Z, int LDZ,
		    double complex *WORK, int LWORK, int *INFO ){
  UNDERSC(zhseqr)( &JOB, &COMPZ, &N, &ILO, &IHI, H, &LDH, W, Z, &LDZ, WORK, &LWORK, INFO );
};

static void zgesvd( char JOBU, char JOBVT, int M, int N,
		    double complex *A, int LDA, double *S,
		    double complex *U, int LDU, double complex *VT, int LDVT,
		    double complex *WORK, int LWORK, double *RWORK,
		    int *INFO ){
  UNDERSC(zgesvd)( &JOBU, &JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT,
		   WORK, &LWORK, RWORK, INFO );
};

static void zheev( char JOBZ, char UPLO, int N, double complex *A, int LDA, double *W,
		   double complex *WORK, int LWORK, double *RWORK, int *INFO ){
  UNDERSC(zheev)( &JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, RWORK, INFO );
};

static void zgetrf( int M, int N, double complex *A, int LDA, int *IPIV,
		    int *INFO ){
  UNDERSC(zgetrf)( &M, &N, A, &LDA, IPIV, INFO );
};

static void dgetrf( int M, int N, double *A, int LDA, int *IPIV, int *INFO ){
  UNDERSC(dgetrf)( &M, &N, A, &LDA, IPIV, INFO );
};

static void zgeqrf( int M, int N, double complex *A, int LDA,
		    double complex *TAU, double complex *WORK, int LWORK, int *INFO ){
  UNDERSC(zgeqrf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, INFO );
};

static void zungqr( int M, int N, int K, double complex *A, int LDA,
		    double complex *TAU, double complex *WORK, int LWORK, int *INFO ){
  UNDERSC(zungqr)( &M, &N, &K, A, &LDA, TAU, WORK, &LWORK, INFO );
};

static void zgeev( char JOBVL, char JOBVR, int N,
		   double complex *A, int LDA,
		   double complex *W,
		   double complex *VL, int LDVL,
		   double complex *VR, int LDVR,
		   double complex *WORK, int LWORK, double *RWORK,
		   int *INFO ){
  UNDERSC(zgeev)( &JOBVL, &JOBVR, &N, A, &LDA, W, VL, &LDVL, VR, &LDVR, WORK, &LWORK, RWORK, INFO );
};

static void dsyev( char JOBZ, char UPLO, int N, double *A, int LDA,
		   double *W, double *WORK, int LWORK, int *INFO ){
  UNDERSC(dsyev)( &JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, INFO );
};
/* *************************************************************************** */
/* ****  Computes eigenvalues of general complex NxN matrix (C-ordered)  ***** */
/* ****  NOTE: original matrix A is destroyed on exit                    ***** */
/* *************************************************************************** */
static uint    N_Hessenberg = 0;
static uint lwork_Hessenberg = 0;
static double complex *work_Hessenberg = NULL;
static double complex *tau_Hessenberg = NULL;

void lapack_cmp_Eigenvalues( double complex *A, uint N , double complex *W ){
  double complex z;
  int i, k;
  if( N != N_Hessenberg ){
    N_Hessenberg = N;
    if( work_Hessenberg ){
      free( work_Hessenberg );
      work_Hessenberg = NULL;
    };
    if( tau_Hessenberg ){
      free( tau_Hessenberg );
      tau_Hessenberg = NULL;
    };
  };
  if( !tau_Hessenberg
      && !(tau_Hessenberg=(double complex *)calloc((N-1), sizeof(double complex)) ))
    panic("Memory");
  if( !work_Hessenberg ){
    zgehrd( N, 1, N, A, N, tau_Hessenberg, &z, -1, &k );
    if( k ) panic("Query of workspace size failed");
    lwork_Hessenberg = (uint) rint( creal(z) );
    if( !(work_Hessenberg=(double complex *)calloc( lwork_Hessenberg, sizeof(double complex))))
      panic("Memory");
  };
  for( i = 0; i < N-1; i++ )
    for( k = i+1; k < N; k++ ){
      z = A[ELEM(i,k,N)];
      A[ELEM(i,k,N)] = A[ELEM(k,i,N)];
      A[ELEM(k,i,N)] = z;
    };
  zgehrd( N, 1, N, A, N, tau_Hessenberg, work_Hessenberg, lwork_Hessenberg, &k );
  if( k ) panic2("Hessenberg reduction failed, lapack return code %d", k );
  zhseqr( 'E', 'N', N, 1, N, A, N, W, &z, 1, work_Hessenberg, lwork_Hessenberg, &k);
  if( k ) panic2("Failed to compute all eigenvalues, lapack returned %d", k );
};
/* *************************************************************************** */
/*
       ZGESVD computes the singular value decomposition (SVD) of a complex  M-
       by-N matrix A, optionally computing the left and/or right singular vec-
       tors. The SVD is written
            A = U * SIGMA * conjugate-transpose(V)

       where SIGMA is an M-by-N matrix which is zero except for  its  min(m,n)
       diagonal  elements,  U  is an M-by-M unitary matrix, and V is an N-by-N
       unitary matrix.  The diagonal elements of SIGMA are the singular values
       of  A;  they  are real and non-negative, and are returned in descending
       order.  The first min(m,n) columns of U and V are the  left  and  right
       singular vectors of A.

       Note that the routine returns V**H, not V.
*/
/* *************************************************************************** */
static uint           M_svd = 0;
static uint           N_svd = 0;
static uint           lwork_svd = 0;
static double complex *U_svd = NULL;
static double complex *work_svd = NULL;
static double        *rwork_svd = NULL;

void lapack_cmp_svd( double complex *A, uint M, uint N, double *W, double complex *V ){
  int i, k, idx, idx2;
  double complex z;
  if( M != M_svd || N != N_svd ){
    M_svd = M;
    N_svd = N;
    lwork_svd = 0;
    if(  work_svd ){  free( work_svd);  work_svd = NULL; };
    if( rwork_svd ){  free(rwork_svd); rwork_svd = NULL; };
    if(     U_svd ){  free(    U_svd);     U_svd = NULL; };
  };
  if( !rwork_svd ){
    i = ( M_svd < N_svd ) ? M_svd : N_svd;
    i *= 5;
    if( !(rwork_svd = (double *) calloc( i, sizeof(double)) ) ) panic("Memory");
  };
  if( !U_svd
      && !(U_svd = (double complex *) calloc( M_svd * N_svd, sizeof(double complex)) ) )
    panic("Memory");
  if( !work_svd ){
    zgesvd( 'O', 'A', M_svd, N_svd, U_svd, M_svd, W, NULL, 1, V, N_svd, &z, -1, rwork_svd, &i );
    if( i ) panic("Query of workspace size failed");
    lwork_svd = (uint) rint( creal(z) );
    if( !(work_svd = (double complex *) calloc( lwork_svd, sizeof(double complex))) )
      panic("Memory");
  };
  for( i = 0; i < M_svd; i++ ){
    idx = ELEM(i, 0, N_svd );
    for( k = 0; k < N_svd; k++ )
      U_svd[ ELEM( k, i, M_svd ) ] = A[ idx + k ];
  };
  bzero( W, N_svd * sizeof(double) );

  zgesvd( 'O', 'A', M_svd, N_svd, U_svd, M_svd, W, NULL, 1, V, N_svd, work_svd, lwork_svd, rwork_svd, &i );
  if( i ) panic2("SVD failed, lapack return code %d", i );

  for( i = 0; i < M_svd; i++ ){
    idx = ELEM(i, 0, N_svd );
    for( k = 0; k < N_svd; k++ )
      A[ idx + k ] = U_svd[ ELEM( k, i, M_svd ) ];
  };
  for( i = 0; i < N_svd - 1; i++ )
    for( k = i+1; k < N_svd; k++ ){
      idx = ELEM(i, k, N_svd );
      idx2= ELEM(k, i, N_svd );
      z = V[idx]; V[idx] = V[idx2];  V[idx2] = z;
    };
};
/* *************************************************************************** */
/* ****  Computes eigenvalues and (optionally) eigenvectors of    ************ */
/* ****  complex Hermitian NxN matrix A (C-ordered)               ************ */
/* *************************************************************************** */
static uint N_eigen = 0;
static uint lwork_eigen = 0;
static double complex *work_eigen = NULL;
static double *rwork_eigen = NULL;

static void lapack_cmp_Eigensystem_real( double complex *A, uint N, double *D , int flag ){
  int i, k, idx, idx2;
  double complex z;
  if( N_eigen != N ){
    N_eigen = N;
    if(  work_eigen ){ free( work_eigen);  work_eigen = NULL; };
    if( rwork_eigen ){ free(rwork_eigen); rwork_eigen = NULL; };
  };
  if( !rwork_eigen
      && !(rwork_eigen = (double *) calloc( 3 * N_eigen - 2, sizeof(double))) )
    panic("Memory");
  if( !work_eigen ){
    zheev( 'V', 'U', N_eigen, A, N_eigen, D, &z, -1, rwork_eigen, &k );
    if( k ) panic("Query of workspace size failed");
    lwork_eigen = (uint) rint( creal(z) );
    if( !(work_eigen = (double complex *) calloc( lwork_eigen, sizeof(double complex))) )
      panic("Memory");
  };
  zheev( (flag) ? 'V' : 'N', 'U', N_eigen, A, N_eigen, D, work_eigen, lwork_eigen, rwork_eigen, &k);
  if( k ) panic2("ZHEEV failed, lapack return code %d", k );

  if( flag )

    for( i = 0; i < N_eigen; i++ ){
      idx  = ELEM( i, i, N_eigen );
      A[idx] = conj( A[idx] );
      for( k = i+1; k < N_eigen; k++ ){
	idx  = ELEM( i, k, N_eigen );
	idx2 = ELEM( k, i, N_eigen );
	z = A[idx]; A[idx] = conj(A[idx2]); A[idx2] = conj(z);
      };
    };
};

void lapack_cmp_Eigensystem( double complex *A, uint N, double *D  ){
  lapack_cmp_Eigensystem_real( A, N, D, 1 );
};
void lapack_cmp_Eigensystem_noeigenvectors( double complex *A, uint N, double *D  ){
  if( N == 2 ){
    double x = creal(A[0]), y = creal(A[3]), z = 2.0 * cabs(A[1]);
#ifdef DEBUG
    if( fabs(cimag(A[0])) > FLT_EPSILON || fabs(cimag(A[3])) > FLT_EPSILON
	|| cabs( A[2] - conj(A[1]) ) > FLT_EPSILON ) panic("Non hermitian input");
#endif
    z *= z;
    D[0] = x - y;
    D[0] *= D[0];
    z = 0.5 * sqrt( z + D[0]);
    D[0] = D[1] = 0.5 * (x + y);
    D[0] += z;
    D[1] -= z;
    return;
  };
  lapack_cmp_Eigensystem_real( A, N, D, 0 );
};
/* *************************************************************************** */
/* Computes an LU factorization of a general M-by-N matrix A  using partial pivoting
   with row interchanges.  The factorization has the form
          A = P * L * U
   where P is a permutation matrix, L is lower triangular with unit diagonal  elements
   (lower  trapezoidal if m > n), and U is upper triangular (upper trapezoidal if m < n).*/
/* *************************************************************************** */
void lapack_cmp_LU( double complex *A, uint M, uint N, uint *perm ){
  int i;
  zgetrf( M, N, A, M, perm, &i );
  if( i ) panic2("ZGETRF failed, lapack return code %d", i );
};
void lapack_dbl_LU( double *A, uint M, uint N, uint *perm ){
  int i;
  dgetrf( M, N, A, M, perm, &i );
  if( i ) panic2("DGETRF failed, lapack return code %d", i );
};
/* *************************************************************************** */
/* ***** Complex and real determinants via above LU factorization   ********** */
/* *************************************************************************** */
static uint N_det = 0;
static uint *perm_det = NULL;

double complex lapack_cmp_determinant( double complex *A, uint N ){
  int i, sign = 1;
  double complex z = 1.0;
  if( N == 1 )      return A[0];
  else if( N == 2 ) return A[0] * A[3] - A[2] * A[1];
  else if( N == 3 ) return
    A[0] * A[4] * A[8] + A[1] * A[5] * A[6] + A[2] * A[3] * A[7] -
    A[6] * A[4] * A[2] - A[7] * A[5] * A[0] - A[8] * A[3] * A[1];
  if( N_det != N ){
    N_det = N;
    if( perm_det ){  free(perm_det); perm_det = NULL; };
  };
  if( !perm_det 
      && !(perm_det = (uint *) calloc( N_det, sizeof(uint) )) ) panic("Memory");
  lapack_cmp_LU( A, N_det, N_det, perm_det );
  for( i = 0; i < N_det; i++ ){
    z *= A[ELEM(i,i,N_det)];
    if( (i+1) != perm_det[i] ) sign *= -1;
  };
  return z * sign;
};
double lapack_dbl_determinant( double *A, uint N ){
  int i, sign = 1;
  double z = 1.0;
  if( N == 1 )      return A[0];
  else if( N == 2 ) return A[0] * A[3] - A[2] * A[1];
  else if( N == 3 ) return
    A[0] * A[4] * A[8] + A[1] * A[5] * A[6] + A[2] * A[3] * A[7] -
    A[6] * A[4] * A[2] - A[7] * A[5] * A[0] - A[8] * A[3] * A[1];
  if( N_det != N ){
    N_det = N;
    if( perm_det ){  free(perm_det); perm_det = NULL; };
  };
  if( !perm_det 
      && !(perm_det = (uint *) calloc( N_det, sizeof(uint) )) ) panic("Memory");
  lapack_dbl_LU( A, N_det, N_det, perm_det );
  for( i = 0; i < N_det; i++ ){
    z *= A[ELEM(i,i,N_det)];
    if( (i+1) != perm_det[i] ) sign *= -1;
  };
  return z * sign;
};
/* *************************************************************************** */
/* *************************************************************************** */
static uint M_qr = 0;
static uint N_qr = 0;
static double complex *tau_qr = NULL;
static double complex *work_qr = NULL;
static uint lwork_qr = 0;
static double complex *U_qr = NULL;

void lapack_zungqr( double complex *A, uint M, uint N ){
  int i, k, idx;
  double complex z;
  if( M < N ) panic("Lapack does not like M<N");
  if( M_qr != M || N_qr != N ){
    M_qr = M;
    N_qr = N;
    if(  tau_qr ){ free(  tau_qr );   tau_qr = NULL; };
    if( work_qr ){ free( work_qr );  work_qr = NULL; };
    if(    U_qr ){ free(    U_qr );     U_qr = NULL; };
  };
  if( !tau_qr
      && !(tau_qr = (double complex *) calloc( (M_qr < N_qr) ? M_qr:N_qr,sizeof(double complex))) )
    panic("Memory");
  if( !U_qr
      && !(U_qr = (double complex *) calloc(M_qr*N_qr,sizeof(double complex))) )
    panic("Memory");
  if( !work_qr ){
    zgeqrf( M_qr, N_qr, U_qr, M_qr, tau_qr, &z, -1, &i );
    if( i ) panic("Query of workspace size failed");
    lwork_qr = (uint) rint( creal(z) );
    if( !(work_qr = (double complex *) calloc( lwork_qr, sizeof(double complex))) )
      panic("Memory");
  };
  for( i = 0; i < M_qr; i++ ){
    idx = ELEM(i, 0, N_qr );
    for( k = 0; k < N_qr; k++ ) U_qr[ ELEM( k, i, M_qr ) ] = A[ idx + k ];
  };
  zgeqrf( M_qr, N_qr, U_qr, M_qr, tau_qr, work_qr, lwork_qr, &i );
  if( i ) panic2("ZGEQRF failed, lapack return code %d", i );
  zungqr( M_qr, N_qr, N_qr, U_qr, M_qr, tau_qr, work_qr, lwork_qr, &i );
  if( i ) panic2("ZUNGQR failed, lapack return code %d", i );
  for( i = 0; i < M_qr; i++ ){
    idx = ELEM(i, 0, N_qr );
    for( k = 0; k < N_qr; k++ ) A[ idx + k ] = U_qr[ ELEM( k, i, M_qr ) ];
  };
};
/* *************************************************************************** */
/* ** Left/Right eigenvectors/eigenvalues of general complex square matrix *** */
/* *************************************************************************** */
static uint N_zgeev = 0;
static double complex *work_zgeev = NULL;
static uint           lwork_zgeev = 0;
static double         *rwork_zgeev = NULL;

static void lapack_zgeev( uchar vl, uchar vr, uint N, 
			  double complex *A, double complex *W,
			  double complex *VL, double complex *VR ){
  int i, k, idx_i, idx_k;
  double complex z;
  if( (vl && !VL) || (vr && !VR) ) panic("Illegal input");
  if( N_zgeev != N ){
    N_zgeev = N;
    if(  work_zgeev ){ free( work_zgeev);  work_zgeev = NULL; };
    if( rwork_zgeev ){ free(rwork_zgeev); rwork_zgeev = NULL; };
  };
  if( !rwork_zgeev 
      && !(rwork_zgeev = (double *) calloc(2*N_zgeev, sizeof(double))) )
    panic("Memory");
  if( !work_zgeev ){
    zgeev( 'N', 'N', N_zgeev, A, N_zgeev, W, NULL, 1, NULL, 1, &z, -1, rwork_zgeev, &i );
    if( i ) panic2("Query of workspace size failed, status %d", i );
    lwork_zgeev = (uint) rint( creal(z) );
    if( !(work_zgeev = (double complex *) calloc( lwork_zgeev, sizeof(double complex))) )
      panic("Memory");
  };
  for( i = 0; i < N_zgeev - 1; i++ )
    for( k = i+1; k < N_zgeev; k++ ){
      idx_i = ELEM( i, k, N_zgeev );
      idx_k = ELEM( k, i, N_zgeev );
      z = A[idx_i];   A[idx_i] = A[idx_k];   A[idx_k] = z;
    };
  zgeev( (vl)? 'V':'N',  (vr)? 'V':'N', N_zgeev, A, N_zgeev, W, VL, N_zgeev, VR, N_zgeev,
	 work_zgeev, lwork_zgeev, rwork_zgeev, &i );
  if( i ) panic2("ZGEEV failed, lapack status is %d", i );

  if( vl || vr )
    for( i = 0; i < N_zgeev - 1; i++ )
      for( k = i+1; k < N_zgeev; k++ ){
	idx_i = ELEM( i, k, N_zgeev );
	idx_k = ELEM( k, i, N_zgeev );
	if( vl ){
	  z = VL[idx_i];   VL[idx_i] = VL[idx_k];   VL[idx_k] = z;
	};
	if( vr ){
	  z = VR[idx_i];   VR[idx_i] = VR[idx_k];   VR[idx_k] = z;
	};
      };
};

void lapack_zgeevl( uint N,  double complex *A, double complex *W, double complex *VL ){
  lapack_zgeev( 1, 0, N, A, W, VL, NULL );
};
void lapack_zgeevr( uint N,  double complex *A, double complex *W, double complex *VR ){
  lapack_zgeev( 0, 1, N, A, W, NULL, VR );
};
void lapack_zge( uint N,  double complex *A, double complex *W ){
  lapack_zgeev( 0, 0, N, A, W, NULL, NULL );
};



/* *************************************************************************** */
/* ****  Computes eigenvalues and (optionally) eigenvectors of    ************ */
/* ****  real symmetric NxN matrix A (C-ordered)                  ************ */
/* *************************************************************************** */
static uint N_dsyev_eigen = 0;
static uint lwork_dsyev_eigen = 0;
static double *work_dsyev_eigen = NULL;

static void lapack_dbl_Eigensystem_real( double *A, uint N, double *D , int flag ){
  int i, k, idx, idx2;
  double z;
  if( N_dsyev_eigen != N ){
    N_dsyev_eigen = N;
    if(  work_dsyev_eigen ){ free( work_dsyev_eigen);  work_dsyev_eigen = NULL; };
  };
  if( !work_dsyev_eigen ){
    dsyev( 'V', 'U', N_dsyev_eigen, A, N_dsyev_eigen, D, &z, -1, &k );
    if( k ) panic("Query of workspace size failed");
    lwork_dsyev_eigen = (uint) rint(z);
    if( !(work_dsyev_eigen = (double *) calloc( lwork_dsyev_eigen, sizeof(double))) )
      panic("Memory");
  };
  dsyev( (flag) ? 'V' : 'N', 'U', N_dsyev_eigen, A, N_dsyev_eigen, D, work_dsyev_eigen, lwork_dsyev_eigen, &k);
  if( k ) panic2("DSYEV failed, lapack return code %d", k );
  if( flag )
    for( i = 0; i < N_dsyev_eigen; i++ )
      for( k = i+1; k < N_dsyev_eigen; k++ ){
	idx  = ELEM( i, k, N_dsyev_eigen );
	idx2 = ELEM( k, i, N_dsyev_eigen );
	z = A[idx]; A[idx] = A[idx2]; A[idx2] = z;
      };
};
void lapack_dbl_Eigensystem( double *A, uint N, double *D  ){
  lapack_dbl_Eigensystem_real( A, N, D, 1 );
};
void lapack_dbl_Eigensystem_noeigenvectors( double *A, uint N, double *D  ){
  if( N == 2 ){
    double x = A[0], y = A[3], z = 2.0 * A[1];
#ifdef DEBUG
    if( fabs( A[2] - A[1] ) > FLT_EPSILON ) panic("Non symmetric input");
#endif
    z *= z;
    D[0] = x - y;
    D[0] *= D[0];
    z = 0.5 * sqrt( z + D[0]);
    D[0] = D[1] = 0.5 * (x + y);
    D[0] += z;
    D[1] -= z;
    return;
  };
  lapack_dbl_Eigensystem_real( A, N, D, 0 );
};
