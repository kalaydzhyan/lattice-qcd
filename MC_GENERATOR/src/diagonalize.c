#include <diagonalize.h>
/*
  Householder reduction of real symmetric NxN matrix A, stored as
  element (i,j) = *( *(A+i) + j)
  to tridiagonal form A'. On output A is replaced by the ortogonal
  matrix Q effecting the pransformation. The relation between matrices:
                     A' = Q^T * A * Q
  D returns diagonal elements and E returns off-diagonal ones with *E=0.
*/
static inline void householder( double **A, double *D, double *E, unsigned int N ){
  int i, j, l, k;
  double h, hh, scale, f, g;
  if( N > 1 ){
    for( i = N-1 ; i > 0; i-- ){
      l = i-1;
      scale = h = 0;
      if( l > 0 ){
	for( k = 0; k < i ; k++ ) scale += fabs(*(*(A+i)+k));
	if( scale == 0.0 ) *(E+i) = *(*(A+i)+l);
	else{
	  for( k = 0; k < i ; k++ ){
	    *(*(A+i)+k) /= scale;
	    h += (*(*(A+i)+k)) * (*(*(A+i)+k));
	  };
	  g = sqrt(h);
	  if( (f = *(*(A+i)+l) ) > 0 ) g*= -1;
	  *(E+i) = scale * g;
	  h -= f*g;
	  *(*(A+i)+l) = f - g;
	  for( f = 0, j = 0; j < i ; j++ ){
	    *(*(A+j)+i) = *(*(A+i)+j) / h;
	    for( g = 0, k = 0 ; k <= j ; k++ ) g+= (*(*(A+j)+k)) * (*(*(A+i)+k));
	    if( l > j ) for( k = j+1 ; k < i ; k++ ) g+= (*(*(A+k)+j)) * (*(*(A+i)+k));
	    *(E+j) = g/h;
	    f += (*(E+j)) * (*(*(A+i)+j));
	  };/* j */
	  hh = f/( h + h );
	  for( j = 0; j < i ; j++ ){
	    f = *(*(A+i)+j);
	    g = *(E+j) - hh * f;
	    *(E+j) = g;
	    for( k = 0; k <= j ; k++ ) *(*(A+j)+k) -= f * (*(E+k)) + g * (*(*(A+i)+k));
	  }; /* j */
	}; /* scale == 0 */
      }else  /* l > 0 */
	*(E+i) = *(*(A+i)+l);
      *(D+i) = h;
    }; /* i */
  }; /* N > 1 */
  /* ------------------------------------------------ */
  *D = *E = 0;
  for( i = 0 ; i < N ; i++ ){
    if( *(D+i) != 0.0  )
      for( j = 0 ; j < i ; j++ ){
	for( g = 0, k = 0 ; k < i ; k++ ) g += (*(*(A+i)+k)) * (*(*(A+k)+j));
	for( k = 0 ; k < i ; k++ ) *(*(A+k)+j) -= g * (*(*(A+k)+i));
      };
    *(D+i) = *(*(A+i)+i);
    *(*(A+i)+i) = 1;
    if( i ) for( j = 0 ; j < i ; j++ ) *(*(A+i)+j) = *(*(A+j)+i) = 0;
  }; /* i */
};
/* ------------------------------------------------------------------ */
static inline void HouseHolder( double *A, double *D, double *E, unsigned int N ){
  int i, j, l, k;
  double h, hh, scale, f, g;
  if( N > 1 ){
    for( i = N-1 ; i > 0; i-- ){
      l = i-1;
      scale = h = 0;
      if( l > 0 ){
	for( k = 0; k < i ; k++ ) scale += fabs(A[MATRIX_ELEMENT(i,k,N)]);
	if( scale == 0.0 ) E[i] = A[MATRIX_ELEMENT(i,l,N)];
	else{
	  for( k = 0; k < i ; k++ ){
	    A[MATRIX_ELEMENT(i,k,N)] /= scale;
	    h += A[MATRIX_ELEMENT(i,k,N)] * A[MATRIX_ELEMENT(i,k,N)];
	  };
	  g = sqrt(h);
	  if( (f = A[MATRIX_ELEMENT(i,l,N)] ) > 0 ) g*= -1;
	  E[i] = scale * g;
	  h -= f*g;
	  A[MATRIX_ELEMENT(i,l,N)] = f - g;
	  for( f = 0, j = 0; j < i ; j++ ){
	    A[MATRIX_ELEMENT(j,i,N)] = A[MATRIX_ELEMENT(i,j,N)] / h;
	    for( g = 0, k = 0 ; k <= j ; k++ ) g+= A[MATRIX_ELEMENT(j,k,N)] * A[MATRIX_ELEMENT(i,k,N)];
	    if( l > j ) for( k = j+1 ; k < i ; k++ ) g+= A[MATRIX_ELEMENT(k,j,N)] * A[MATRIX_ELEMENT(i,k,N)];
	    E[j] = g/h;
	    f += E[j] * A[MATRIX_ELEMENT(i,j,N)];
	  };/* j */
	  hh = f/( h + h );
	  for( j = 0; j < i ; j++ ){
	    f = A[MATRIX_ELEMENT(i,j,N)];
	    g = E[j] - hh * f;
	    E[j] = g;
	    for( k = 0; k <= j ; k++ ) 
	      A[MATRIX_ELEMENT(j,k,N)] -= f * E[k] + g * A[MATRIX_ELEMENT(i,k,N)];
	  }; /* j */
	}; /* scale == 0 */
      }else  /* l > 0 */
	E[i] = A[MATRIX_ELEMENT(i,l,N)];
      D[i] = h;
    }; /* i */
  }; /* N > 1 */
  /* ------------------------------------------------ */
  D[0] = E[0] = 0.0;
  for( i = 0 ; i < N ; i++ ){
    if( D[i] != 0.0  )
      for( j = 0 ; j < i ; j++ ){
	for( g = 0, k = 0 ; k < i ; k++ ) g += A[MATRIX_ELEMENT(i,k,N)] * A[MATRIX_ELEMENT(k,j,N)];
	for( k = 0 ; k < i ; k++ ) A[MATRIX_ELEMENT(k,j,N)] -= g * A[MATRIX_ELEMENT(k,i,N)];
      };
    D[i] = A[MATRIX_ELEMENT(i,i,N)];
    A[MATRIX_ELEMENT(i,i,N)] = 1;
    if( i ) for( j = 0 ; j < i ; j++ ) A[MATRIX_ELEMENT(i,j,N)] = A[MATRIX_ELEMENT(j,i,N)] = 0.0;
  }; /* i */
};
/* =================================================================================== */
/*
  QL with implicit shifts, allows to diagonalize tridiagonal real symmetric NxN input
  matrix A, stored as: element (i,j) = *( *(A+i) + j). Finds eigenvalues and eigenvectors.
  Input vector D of length N is the diagonal elements of A, on output contains eigenvalues.
  Input vector E of length N is off-diagonal elements with *E arbitrary, on output it is
  destroyed. Matrix of eigenvectors Q is returned via A. Relation between matrices:
      diag( D ) = Q^T * A * Q
  Return values: 0 for success and -1 if iteration count is too large.
*/
static inline int QLshifts( double **A, double *D, double *E, unsigned int N ){
  int i,l,k,iter,m;
  double dd, g, s,c,p, f,b, r;
  if( N > 1 ){
    for( i = 1; i < N; i++ ) *(E+i-1) = *(E+i);
    *(E+N-1) = 0;
    for( l = 0; l < N; l++ ){
      iter = 0;
    label1: for( m = l; m < N-1 ; m++){
	dd = fabs(*(D+m)) + fabs(*(D+m+1));
	if( (fabs(*(E+m))+dd) == dd  ) goto label2;
      };
      m=N-1;
    label2: if( m != l ){
	if( iter++ == 30 ) return -1;
	g = (*(D+l+1) - *(D+l))/(2*(*(E+l)));
	g = ( g > 0 ) ? (g + sqrt(g*g+1)) : (g - sqrt(g*g+1));
	g = *(D+m) - *(D+l) + (*(E+l))/g;
	s = c = 1.0;
	p = 0.0;
	for( i = m-1 ; i >= l; i--){
	  f = s * (*(E+i));
	  b = c * (*(E+i));
	  if( fabs(f) >= fabs(g) ){
	    c = g/f;
	    r = sqrt(c*c+1);
	    s = 1.0/r;
	    *(E+i+1) = f*r;
	    c *= s;
	  }else{
	    s = f/g;
	    r = sqrt(s*s+1);
	    c = 1.0/r;
	    *(E+i+1) = g*r;
	    s *= c;
	  };
	  g = *(D+i+1) - p;
	  r = (*(D+i) - g) * s + 2 * c * b;
	  p = r * s;
	  *(D+i+1) = g + p;
	  g = c * r - b;
	  // --------------------
	  for( k = 0; k < N; k++ ){
	    f = *(*(A+k)+i+1);
	    *(*(A+k)+i+1) = s * (*(*(A+k)+i)) + c*f;
	    *(*(A+k)+i)   = c * (*(*(A+k)+i)) - s*f;
	  };
	  // --------------------
	}; /* i */
	*(D+l) -= p;
	*(E+l) = g;
	*(E+m) = 0;
	goto label1;
      }; /* m != l */
    }; /* l */
  }; /* N > 1 */
  return 0;
};
/* --------------------------------------------------------- */
static inline int QL_Shifts( double *A, double *D, double *E, unsigned int N ){
  int i,l,k,iter,m;
  double dd, g, s,c,p, f,b, r;
  if( N > 1 ){
    for( i = 1; i < N; i++ ) E[i-1] = E[i];
    E[N-1] = 0;
    for( l = 0; l < N; l++ ){
      iter = 0;
    label1: for( m = l; m < N-1 ; m++){
	dd = fabs(D[m]) + fabs(D[m+1]);
	if( (fabs(E[m])+dd) == dd  ) goto label2;
      };
      m=N-1;
    label2: if( m != l ){
	if( iter++ == 30 ) return -1;
	g = (D[l+1] - D[l])/(2*E[l] );
	g = ( g > 0 ) ? (g + sqrt(g*g+1)) : (g - sqrt(g*g+1));
	g = D[m] - D[l] + E[l]/g;
	s = c = 1.0;
	p = 0.0;
	for( i = m-1 ; i >= l; i--){
	  f = s * E[i];
	  b = c * E[i];
	  if( fabs(f) >= fabs(g) ){
	    c = g/f;
	    r = sqrt(c*c+1);
	    s = 1.0/r;
	    E[i+1] = f*r;
	    c *= s;
	  }else{
	    s = f/g;
	    r = sqrt(s*s+1);
	    c = 1.0/r;
	    E[i+1] = g*r;
	    s *= c;
	  };
	  g = D[i+1] - p;
	  r = (D[i] - g) * s + 2 * c * b;
	  p = r * s;
	  D[i+1] = g + p;
	  g = c * r - b;
	  // --------------------
	  for( k = 0; k < N; k++ ){
	    f = A[MATRIX_ELEMENT(k,i+1,N)];
	    A[MATRIX_ELEMENT(k,i+1,N)] = s * A[MATRIX_ELEMENT(k,i,N)] + c*f;
	    A[MATRIX_ELEMENT(k,i  ,N)] = c * A[MATRIX_ELEMENT(k,i,N)] - s*f;
	  };
	  // --------------------
	}; /* i */
	D[l] -= p;
	E[l] = g;
	E[m] = 0;
	goto label1;
      }; /* m != l */
    }; /* l */
  }; /* N > 1 */
  return 0;
};
/* =================================================================================== */
/*
  Performs diagonalization of given real symmetric NxN matrix A, finds its eigenvectors and
  eigenvalues. A is stored as: element (i,j) = *( *(A+i) + j). Eigenvalues are returned via
  vector D of length N. On output A is replaced with matrix Q of eigenvectors. Relation betweeen
  matrices:
                 diag( D ) = Q^T * A * Q

  Return values:  0 for success and -1 for errors (iteration count is too large,
  memory allocation failed).
*/
int RealSymmetricDiagonalize( double **A, double *D, unsigned int N ){
  int ret;
  double *E = (double *) calloc( N, sizeof(double));
  if( !E ) return -1;
  householder( A, D, E, N );
  ret = QLshifts( A, D, E, N );
  free(E);
  return ret;
};
int Real_Symmetric_Diagonalize( double *A, double *D, unsigned int N ){
  int ret;
  double *E = (double *) calloc( N, sizeof(double));
  if( !E ) return -1;
  HouseHolder( A, D, E, N );
  ret = QL_Shifts( A, D, E, N );
  free(E);
  return ret;
};
