#include <svd.h>

/*  Routine is likely to be used frequently, so we want to avoid
    frequent malloc/free of internal vars */
static uint N_internal = 0;
static double *rv = NULL;
static double complex *theta = NULL;
static double complex *chi = NULL;
#ifdef DEBUG
static uint M_internal = 0;
static double complex *save = NULL;
#endif

#define   M_EL(i,j,N)    ((i) * (N) + (j))
#define   SVD_SIGN(a,b)  (((b) >= 0.0) ? fabs(a) : -fabs(a)) 

/* -------------------------------------------------------------------------- */
static inline double pythag( double a, double b ){
  double sq, mod_a = fabs(a), mod_b = fabs(b);
  if( mod_a > mod_b ){
    sq = mod_b/mod_a;
    return mod_a * sqrt( 1.0 + sq * sq );
  };
  if( mod_b < DBL_EPSILON ) return 0.0;
  sq = mod_a / mod_b;
  return mod_b * sqrt( 1.0 + sq * sq );
};
/* -------------------------------------------------------------------------- */
void cmp_svd( double complex *A, uint M, uint N, double *W, double complex *V ){
  int i, j, k, l=0, flag, its, jj, nm=0, idx, idx_d;
  double c, f, g, h, s, x, y, z, anorm, scale;
  double complex zeta, xi;
  if( !M || !N || !A || !W || !V ) panic("NULL input");
  if( N > N_internal ){
    N_internal = N;
    if(  rv   ){ free( rv  );    rv = NULL; };
    if( theta ){ free(theta); theta = NULL; };
    if(  chi  ){ free( chi );   chi = NULL; };
  };
#ifdef DEBUG
  if( N > N_internal || M > M_internal ){
    N_internal = N;
    M_internal = M;
    if(  save ){ free( save ); save = NULL; };
  };
#endif
  if( !rv
      && !(rv    = (double *)         calloc( N, sizeof(double)         )) ) panic("Memory");
  if( !theta
      && !(theta = (double complex *) calloc( N, sizeof(double complex) )) ) panic("Memory");
  if( !chi
      && !(chi   = (double complex *) calloc( N, sizeof(double complex) )) ) panic("Memory");
#ifdef DEBUG
  if( !save
      && !(save   = (double complex *) calloc( N * M, sizeof(double complex) )) ) panic("Memory");
  bcopy( A, save, N * M * sizeof(double complex));
#endif
  g = scale = anorm = 0.0;
  for( i = 0; i < N; i++ ){
    l = i + 1;
    idx_d = M_EL(i,i,N);
    theta[i] = chi[i] = 1.0;
    rv[i] = scale * g;
    g = s = scale = 0.0;
    if( i < M ){  /* Housholder reduction from the left to remove
		     everything below diagonal */
      for( k = i; k < M; k++ ) scale += cabs( A[M_EL(k,i,N)] );
      if( scale != 0.0 ){
	for( k = i; k < M; k++ ){
	  idx = M_EL(k,i,N);
	  A[idx] /= scale;
	  f = cabs( A[idx] );
	  s += f * f;
	};
	f = cabs( A[idx_d] );
	if( f < DBL_EPSILON ) theta[i] = -1.0;
	else                  theta[i] = - A[idx_d] / f;
	g = sqrt(s);
	h = f * g + s;
	A[idx_d] = f + g;
	for( j = l; j < N; j++ ){ /* Reduce matrix */
	  idx = M_EL(i,j,N);
	  zeta = - conj(theta[i]) * A[idx_d] * A[idx];
	  for( k = l; k < M; k++ ) zeta += conj(A[M_EL(k,i,N)]) * A[M_EL(k,j,N)];
	  zeta /= h;
	  A[idx] += zeta * theta[i] * A[idx_d];
	  for( k = l; k < M; k++ ) A[M_EL(k,j,N)] -= zeta * A[M_EL(k,i,N)];
	};
	for( k = i; k < M; k++ ) A[M_EL(k,i,N)] *= scale;
      };
    };
    W[i] = scale * g;
    g = s = scale = 0.0;
    if( i < M && i != (N - 1) ){ /* Housholder reduction from the right to remove
				    everything to the right of super-diagonal */
      idx = M_EL(i,0,N);
      idx_d = idx + l;
      for( k = l; k < N; k++ ) scale += cabs( A[idx + k] );
      if( scale != 0.0 ){
	for( k = l; k < N; k++ ){
	  j = idx + k;
	  A[j] /= scale;
	  f = cabs( A[j] );
	  s += f * f;
	};
	f = cabs( A[idx_d] );
	if( f < DBL_EPSILON ) chi[i] = -1.0;
	else                  chi[i] = - A[idx_d] / f;
	g = sqrt(s);
	h = f * g + s;
	A[idx_d] = f + g;
	for( j = l; j < M; j++ ){
	  idx = M_EL(j,l,N);
	  zeta = - conj( chi[i] ) * A[idx_d] * A[idx];
	  for( k = l + 1; k < N; k++ ) zeta += A[M_EL(j,k,N)] * conj( A[M_EL(i,k,N)] );
	  zeta /= h;
	  A[idx] +=  zeta * A[idx_d] * chi[i];
	  for( k = l + 1; k < N; k++ ) A[M_EL(j,k,N)] -= zeta * A[M_EL(i,k,N)];
	};
	for( k = l; k < N; k++ ) A[M_EL(i,k,N)] *= scale;
      };
    };
    if( (f = fabs(W[i]) + fabs(rv[i])) > anorm ) anorm = f;
  };
  /* Accumulating right-hand transformations */
  for( i = N - 1; i >= 0 ; i-- ){
    if( i < N - 1 ){
      if( g != 0.0 ){
	idx = M_EL(i,l,N);
	V[M_EL(l,i,N)] = - conj(chi[i])/g;
	for( j = l+1; j < N; j++ ) /* :-( ; besides, A[i][l] is real */
	  V[M_EL(j,i,N)] = ( conj(A[M_EL(i,j,N)]) / A[idx] ) / g;
	for( j = l; j < N; j++ ){
	  zeta = - chi[i] * A[idx] * V[M_EL(l,j,N)];
	  for( k = l+1; k < N; k++ ) zeta += A[M_EL(i,k,N)] * V[M_EL(k,j,N)];
	  for( k = l; k < N; k++ ) V[M_EL(k,j,N)] -= zeta * V[M_EL(k,i,N)];
	};
      };
      for( j = l; j < N; j++ ) V[M_EL(i,j,N)] = V[M_EL(j,i,N)] = 0.0;
    };
    V[M_EL(i,i,N)] = 1.0;
    g = rv[i];
    l = i;
  };
  for( zeta = 1.0, j = 1; j < N; j++ ){ /* Take care of complex phases */
    zeta *= theta[j-1] * conj(chi[j-1]);
    for( i = 0; i < N; i++ ) V[M_EL(i,j,N)] *= zeta;
    chi[j-1] = zeta; /* Numbering shifted by one! */
  };
  /* Accumulating left-hand transformations */
  i = (N < M) ? N - 1 : M - 1;
  for( ; i >= 0; i-- ){
    idx_d = M_EL(i,i,N);
    l = i+1;
    g = W[i];
    idx = M_EL(i,0,N);
    for( j = l; j < N; j++ ) A[idx+j] = 0.0;
    if( g != 0.0 ){
      g = 1.0/g;
      for( j = l; j < N; j++ ){
	zeta = 0.0;
	for( k = l; k < M; k++ ){
	  idx = M_EL(k,0,N);
	  zeta += A[idx+j] * conj(A[idx+i]);
	};
	zeta = (zeta/A[idx_d]) * g;
	A[M_EL(i,j,N)] = zeta * A[idx_d] * theta[i];
	for( k = l; k < M; k++ ){
	  idx = M_EL(k,0,N);
	  A[idx+j] -= zeta * A[idx+i];
	};
      };
      A[idx_d] *= -g;
      for( j = l; j < M; j++ ) A[M_EL(j,i,N)] *= conj(theta[i]) * g;
    }else{
      for( j = i; j < M; j++ ) A[M_EL(j,i,N)] = 0.0;
    };
    A[idx_d] += 1.0;
  };
  for( j = 0; j < N; j++ ) /* Take care of complex phases */
    for( i = 0; i < M; i++ ){
      idx = M_EL(i,j,N);
      A[idx] *= theta[j];
      if(j) A[idx] *= chi[j-1];
    };
#ifdef DEBUG
  for( i = 0; i < M; i++ )
    for( j = 0; j < N; j++ ){
      zeta = 0.0;
      for( k = 0; k < N; k++ ){
	zeta += A[M_EL(i,k,N)] * conj(V[M_EL(j,k,N)]) * W[k];
	if( k != N - 1 ) zeta += A[M_EL(i,k,N)] * conj(V[M_EL(j,k+1,N)]) * rv[k+1];
      };
      if( cabs(zeta - save[M_EL(i,j,N)]) > FLT_EPSILON )
	panic2("Bidiagonal already failed, %.10f", cabs(zeta - save[M_EL(i,j,N)]) );
    };
#endif
  /* Problem was reduced to diagonalization of REAL bidiagonal matrix.
     Code below is stolen from NR */
  for( k = N - 1; k >= 0 ; k-- ){
    for( its = 0; its < SVD_MAX_ITER; its++ ){
      flag = 1;
      for( l = k; l >= 0; l-- ){
	nm = l - 1;
	if( fabs(rv[l]) < DBL_EPSILON ){
	  flag = 0;
	  break;
	};
	if( fabs(W[nm]) < DBL_EPSILON ) break;
      };
      if( flag ){
	c = 0.0;
	s = 1.0;
	for( i = l; i <= k; i++ ){
	  f = s * rv[i];
	  rv[i] *= c;
	  if( fabs(f) < DBL_EPSILON ) break;
	  g = W[i];
	  h = pythag( f, g );
	  W[i] = h;
	  h = 1.0/h;
	  c = g * h;
	  s = -f * h;
	  for( j = 0; j < M; j++ ){
	    idx_d = M_EL(j,nm,N);
	    idx   = M_EL(j,i,N);
	    xi   = A[idx_d];
	    zeta = A[idx];
	    A[idx_d] = xi * c + zeta * s;
	    A[idx]   = zeta * c - xi * s;
	  };
	}; /* i cycle */
      }; /* if( flag ) */
      z = W[k];
      if( l == k ){
	if( z < 0.0 ){
	  W[k] = -z;
	  for( j = 0; j < N; j++ ) V[M_EL(j,k,N)] *= -1.0;
	};
	break;
      };
      if( (its + 1) == SVD_MAX_ITER )
	panic2("SVD didn't converge in %d iterations", SVD_MAX_ITER );
      x = W[l];
      nm = k - 1;
      y = W[nm];
      g = rv[nm];
      h = rv[k];
      f = ( (y-z)*(y+z) + (g-h)*(g+h) )/(2.0 * h * y);
      g = pythag( f, 1.0 );
      f = ((x-z)*(x+z) + h * (( y/(f+SVD_SIGN(g,f)) ) - h))/x;
      c = s = 1.0;
      for( j = l; j <= nm; j++ ){
	i = j + 1;
	g = rv[i];
	y = W[i];
	h = s * g;
	g *= c;
	z = pythag(f,h);
	rv[j] = z;
	c = f/z;
	s = h/z;
	f = x * c + g * s;
	g = g * c - x * s;
	h = y * s;
	y *= c;
	for( jj = 0; jj < N; jj++ ){
	  idx_d = M_EL( jj, j, N );
	  idx   = M_EL( jj, i, N );
	  xi   = V[idx_d];
	  zeta = V[idx];
	  V[idx_d] = xi * c + zeta * s;
	  V[idx]   = zeta * c - xi * s;
	};
	z = pythag( f, h );
	W[j] = z;
	if( z != 0.0 ){
	  z = 1.0 / z;
	  c = f * z;
	  s = h * z;
	};
	f = c * g + s * y;
	x = c * y - s * g;
	for( jj = 0; jj < M; jj++ ){
	  idx_d = M_EL( jj, j, N );
	  idx   = M_EL( jj, i, N );
	  xi   = A[idx_d];
	  zeta = A[idx];
	  A[idx_d] = xi * c + zeta * s;
	  A[idx]   = zeta * c - xi * s;
	};
      }; /* j cycle */
      rv[l] = 0.0;
      rv[k] = f;
      W[k] = x;
    }; /* its cycle */
  }; /* k cycle */
#ifdef DEBUG
  for( i = 0; i < N; i++ )
    for( j = 0; j < N; j++ ){
      zeta = ( i==j ) ? -1.0 : 0.0;
      for( k = 0; k < M; k++ ) zeta += conj(A[M_EL(k,i,N)]) * A[M_EL(k,j,N)];
      if( cabs(zeta) > FLT_EPSILON ) panic("A is not unitary on output");
    };
  for( i = 0; i < N; i++ )
    for( j = 0; j < N; j++ ){
      zeta = ( i==j ) ? -1.0 : 0.0;
      for( k = 0; k < N; k++ ) zeta += conj(V[M_EL(k,i,N)]) * V[M_EL(k,j,N)];
      if( cabs(zeta) > FLT_EPSILON ) panic("V is not unitary on output");
    };
  for( i = 0; i < M; i++ )
    for( j = 0; j < N; j++ ){
      zeta = 0.0;
      for( k = 0; k < N; k++ )
	zeta += A[M_EL(i,k,N)] * conj(V[M_EL(j,k,N)]) * W[k];
      if( cabs(zeta - save[M_EL(i,j,N)]) > FLT_EPSILON )
	panic2("Diagonalization failed, %.10f", cabs(zeta - save[M_EL(i,j,N)]) );
    };
#endif
  for( i = 0; i < N; i++ ){
    idx = M_EL( i, 0, N );
    for( j = i+1; j < N; j++ ){
      idx_d = M_EL( j, i, N );
      zeta = V[idx+j];
      V[idx+j] = conj(V[idx_d]);
      V[idx_d] = conj(zeta);
    };
    V[idx+i] = conj(V[idx+i]);
  };
};
/* ********************************************************************************* */
/*
  Householder reduction of hermitian matrix A[N][N] (first arg - &A[0][0]) to
  tridiagonal form A'. On output A is replaced by unitary matrix Q effecting
  the transformation. The relation between matrices:
                     A' = Q^+ * A * Q
  D returns diagonal elements and E returns super-diagonal ones with E[0] = 0.
  flag - if zero we skip eigenvectors calculation.
*/

#define   ELEM(i,j,N)   ((i) * (N) + (j))

static uint N_hs = 0;
static double complex *hs_E = NULL;

static void cmp_householder_real( double complex *A, uint N, double *D, double *E, int flag ){
  int i, j, k, l, idx, idx_d;
  double scale, h, g, f;
  double complex z, zeta;
#ifdef DEBUG
  double complex *orig = NULL;
  for( i = 0; i < N; i++ )
    for( j = i; j < N; j++ )
      if( cabs( A[ELEM(i,j,N)] - conj(A[ELEM(j,i,N)])) > DBL_EPSILON ) panic("Non hermitian input");
  if( flag ){
    if( !(orig = (double complex *) calloc( N * N, sizeof(double complex) )) ) panic("Memory");
    bcopy( A, &orig[0], N * N * sizeof(double complex) );
  };
#endif
  if( N > N_hs ){
    N_hs = N;
    if( hs_E ){ free( hs_E ); hs_E = NULL; };
  };
  if( !hs_E 
      && !( hs_E = (double complex *) calloc( N, sizeof(double complex)) )) panic("Memory");
  for( i = N - 1; i > 0; i-- ){
    l = i - 1;
    h = scale = 0.0;
    idx_d = ELEM(i,0,N);
    if( l ){
      for( k = 0; k < i; k++ ) scale += cabs( A[idx_d+k] );
      if( scale < DBL_EPSILON ){
	hs_E[i] = conj(A[idx_d + l]); /* super-diagonal */
      }else{
	for( k = 0; k < i; k++ ){
	  idx = idx_d + k;
	  A[idx] /= scale;
	  f = cabs( A[idx] );
	  h += f * f;
	};
	g = sqrt(h);
	idx = idx_d + l;
	z = A[idx];
	f = cabs(z);
	if( f < DBL_EPSILON )   z  = 1.0;
	else                    z /= f;
	hs_E[i] = - scale * g * conj(z); /* super-diagonal */
	h += g * f;
	A[idx] = z * (g + f); /* Store u^+ in i-th row of A */
	zeta = 0.0;
 	for( j = 0; j < i; j++ ){
	  idx = ELEM(j,0,N);
	  if( flag ) A[idx+i] = conj(A[idx_d + j])/h; /* Store u/H in i-th column of A */
	  z = 0.0;    /* form j-th element of u^+ A; note that upper triangle was just spoiled */
	  for( k = 0; k <= j ; k++ )  z += A[idx_d + k] * conj( A[  idx + k  ] );
	  for( k = j+1; k < i; k++ )  z += A[idx_d + k] *       A[ELEM(k,j,N)] ;
	  hs_E[j] = z/h; /* form p^+ in temporaly unused place */
	  zeta += hs_E[j] * conj(A[idx_d + j]); /* this is p^+ u */
	};
	zeta /= h + h;
 	for( j = 0; j < i; j++ ){
	  idx = ELEM(j,0,N);
	  z = conj(A[idx_d+j]); /* j-th component of u (not u^+) */
	  hs_E[j] = conj(hs_E[j]) - z * zeta; /* j-th component of q (not q^+) */
	  for( k = 0; k <= j ; k++ )
	    A[idx+k] -= hs_E[j] * A[idx_d + k] + z * conj(hs_E[k]);
	};	
      }; /* if scale */
    }else
      hs_E[i] = conj( A[idx_d + l] ); /* super-diagonal */
    D[i] = h;
  }; /* i */
  /* ----------------------------------------------------- 
     Matrix had been reduced to complex hermitian tridiagonal form.
     Collect the corresponding unitary transformation. */
  E[0] = D[0] = 0.0;
  if( flag ){
    for( i = 0; i < N; i++ ){
      idx_d = ELEM(i,0,N);
      idx = idx_d + i;
      if( D[i] != 0.0 )
	for( j = 0; j < i; j++ ){
	  z = 0.0;
	  for( k = 0; k < i; k++ ) z += A[idx_d + k] * A[ELEM(k,j,N)];
	  for( k = 0; k < i; k++ ) A[ELEM(k,j,N)] -= A[ELEM(k,i,N)] * z;
	};
#ifdef DEBUG
      if( fabs(cimag(A[idx])) > DBL_EPSILON ) panic("Imaginary diagonal element");
#endif
      D[i] = creal( A[idx] );
      A[idx] = 1.0;
      for( j = 0; j < i; j++ ) A[idx_d + j] = A[ELEM(j,i,N)] = 0.0;
    }; /* i */
  }else{
    for( i = 0; i < N; i++ ) D[i] = creal( A[ELEM(i,i,N)] );
  };
#ifdef DEBUG
  if( flag ){
    for( i = 0; i < N; i++ )
      for( j = 0; j < N; j++ ){
	z = ( i==j ) ? -1.0 : 0.0;
	for( k = 0; k < N; k++ ) z += conj(A[ELEM(k,i,N)]) * A[ELEM(k,j,N)];
	if( cabs(z) > FLT_EPSILON ) panic("A is not unitary");
      };
    for( i = 0; i < N; i++ )
      for( j = 0; j < N; j++ ){
	z = 0.0;
	for( k = 0; k < N; k++ ){
	  z += A[ELEM(i,k,N)] * conj(A[ELEM(j,k,N)]) * D[k];
	  if( k != N - 1 ){
	    z += A[ELEM(i,k,N)] * conj(A[ELEM(j,k+1,N)]) * hs_E[k+1];
	    z += A[ELEM(i,k+1,N)] * conj(A[ELEM(j,k,N)]) * conj(hs_E[k+1]);
	  };
	};
	if( cabs(orig[ELEM(i,j,N)] - z) > FLT_EPSILON ){
	  fprintf(stderr, "----------- [%d,%d] ---------------\n", i, j);
	  fprintf(stderr, "Orig: (%f, %f)\n", creal(orig[ELEM(i,j,N)]), cimag(orig[ELEM(i,j,N)]) );
	  fprintf(stderr, "Got : (%f, %f)\n", creal(z), cimag(z) );
	  panic("Complex three-diagonalization FAILED");
	};
      };
  };
#endif
  /* ----------------------------------------------------------
     Remove remaining complex entries in sub/super diagonals */
  for( z = 1.0, j = 1; j < N; j++ ){
    E[j] = cabs(hs_E[j]);
    if( flag ){
      if( E[j] != 0.0 ) z *= conj( hs_E[j] )/E[j];
      for( i = 0; i < N; i++ ) A[ELEM(i,j,N)] *= z;
    };
  };
#ifdef DEBUG
  if( flag ){
    for( i = 0; i < N; i++ )
      for( j = 0; j < N; j++ ){
	z = ( i==j ) ? -1.0 : 0.0;
	for( k = 0; k < N; k++ ) z += conj(A[ELEM(k,i,N)]) * A[ELEM(k,j,N)];
	if( cabs(z) > FLT_EPSILON ) panic("A is not unitary");
      };
    for( i = 0; i < N; i++ )
      for( j = 0; j < N; j++ ){
	z = 0.0;
	for( k = 0; k < N; k++ ){
	  z += A[ELEM(i,k,N)] * conj(A[ELEM(j,k,N)]) * D[k];
	  if( k != N - 1 ){
	    z += A[ELEM(i,k,N)] * conj(A[ELEM(j,k+1,N)]) * E[k+1];
	    z += A[ELEM(i,k+1,N)] * conj(A[ELEM(j,k,N)]) * E[k+1];
	  };
	};
	if( cabs(orig[ELEM(i,j,N)] - z) > FLT_EPSILON ){
	  fprintf(stderr, "----------- [%d,%d] ---------------\n", i, j);
	  fprintf(stderr, "Orig: (%f, %f)\n", creal(orig[ELEM(i,j,N)]), cimag(orig[ELEM(i,j,N)]) );
	  fprintf(stderr, "Got : (%f, %f)\n", creal(z), cimag(z) );
	  panic("Three-diagonalization FAILED");
	};
      };
    free(orig);
  };
#endif
};
/* ********************************************************************************* */
void cmp_householder( double complex *A, uint N, double *D, double *E ){
  cmp_householder_real( A, N, D, E, 1 );
};
/* ********************************************************************************* */
void cmp_householder_noeigenvectors( double complex *A, uint N, double *D, double *E ){
  cmp_householder_real( A, N, D, E, 0 );
};
/* ********************************************************************************* */
/*
  QL with implicit shifts, allows to diagonalize tridiagonal real symmetric NxN input
  matrix A, stored in C-order. Finds eigenvalues and eigenvectors.
    D (dim D = N) -- diagonal elements of A, on output contains eigenvalues.
    E (dim E = N) -- sub/super-diagonal elements of A with E[0] - undefined,
                     on output is destroyed.
  Matrix of eigenvectors Q is returned via A with:
      diag( D ) = A^T * [original matrix] * A
  flag - if zero we skip eigenvectors calculation
*/
static void cmp_QLShifts_real( double complex *A, uint N, double *D, double *E, int flag ){
  int i, k, m, l, iter, idx;
  double dd, g, s,c,p, f,b, r;
  double complex z;
#ifdef DEBUG
  double *rotation = NULL, *d = NULL, *e = NULL;
  if( !(d = (double *)calloc( N , sizeof(double))) ) panic("Memory");
  if( !(e = (double *)calloc( N , sizeof(double))) ) panic("Memory");
  if( !(rotation = (double *)calloc( N * N, sizeof(double))) ) panic("Memory");
  bcopy( &D[0], &d[0], N * sizeof(double));
  bcopy( &E[0], &e[0], N * sizeof(double));
  bzero( rotation, N * N * sizeof(double));
  for( i = 0; i < N; i++ ) rotation[ELEM(i,i,N)] = 1.0;
#endif
  for( i = 1; i < N; i++ ) E[i-1] = E[i];
  E[N-1] = 0.0;
  for( l = 0; l < N; l++ ){
    iter = 0;
    do{
      for( m = l; m < N - 1; m++ ){
	dd = fabs(D[m]) + fabs(D[m+1]);
	if( fabs(E[m]) < DBL_EPSILON ) break;
      };
      if( m != l ){
	if( iter++ == SVD_MAX_ITER ) panic("Too many iterations in QLShifts");
	g = 0.5 * (D[l+1] - D[l])/E[l];
	r = pythag( g, 1.0 );
	g = D[m] - D[l] + E[l]/(g + SVD_SIGN(r,g));
	s = c = 1.0;
	p = 0.0;
	for( i = m-1; i >= l; i-- ){
	  f = s * E[i];
	  b = c * E[i];
	  E[i+1] = (r = pythag( f, g ));
	  if( r < DBL_EPSILON ){
	    D[i+1] -= p;
	    E[m] = 0.0;
	    break;
	  };
	  s = f/r;
	  c = g/r;
	  g = D[i+1] - p;
	  r = (D[i] - g) * s + 2.0 * c * b;
	  D[i+1] = g + ( p = s * r );
	  g = c * r - b;
	  if( flag ){
	    for( k = 0; k < N; k++ ){
	      idx = ELEM( k, i, N );
	      z = A[idx+1];
	      A[idx+1] = s * A[idx] + c * z;
	      A[idx]   = c * A[idx] - s * z;
	    };
	  };
#ifdef DEBUG
	  for( k = 0; k < N; k++ ){
	    idx = ELEM( k, i, N );
	    f = rotation[idx+1];
	    rotation[idx+1] = s * rotation[idx] + c * f;
	    rotation[idx]   = c * rotation[idx] - s * f;
	  };
#endif
	}; /* i */
	if( r < DBL_EPSILON && i >= l ) continue;
	D[l] -= p;
	E[l] = g;
	E[m] = 0.0;
      };
    }while( m != l );
  }; /* l */
#ifdef DEBUG
  for( i = 0; i < N; i++ )
    for( k = 0; k < N; k++ ){
      g = 0.0;
      if( i == k ) g = -d[i];
      if( k == i + 1 ) g = -e[k];
      if( k == i - 1 ) g = -e[i];
      for( m = 0; m < N; m++ )
	g += rotation[ELEM(i,m,N)] * D[m] * rotation[ELEM(k,m,N)];
      if( fabs(g) > FLT_EPSILON ) panic("QLShifts failed");
    };
  free( rotation);
  free(d);
  free(e);
#endif
};
/* ********************************************************************************* */
void cmp_QLShifts( double complex *A, uint N, double *D, double *E ){
  cmp_QLShifts_real( A, N, D, E, 1 );
};
/* ********************************************************************************* */
void cmp_QLShifts_noeigenvectors( double complex *A, uint N, double *D, double *E ){
  cmp_QLShifts_real( A, N, D, E, 0 );
};
/* ********************************************************************************* */
/*
  Wrapper for two above routines - finds eigensystem of complex hermitian matrix A
  stored in C-order. On output A is replaced with eigenvectors, D containes eigenvalues.
  Relation between matrices on output:
        diag[ D ] = A^+ * [original matrix] * A
  flag - if zero we skip eigenvectors calculation
*/
static uint    N_eigen = 0;
static double *E_eigen = NULL;

static void cmp_Eigensystem_real( double complex *A, uint N, double *D , int flag ){
#ifdef DEBUG
  int i,j,k;
  double complex z, *save = NULL;
  if( flag ){
    if( !(save = (double complex *)calloc( N * N, sizeof(double complex))) )
      panic("Memory");
    bcopy( A, save, N * N * sizeof(double complex) );
  };
#endif
  if( N > N_eigen ){
    N_eigen = N;
    if( E_eigen ){
      free( E_eigen );
      E_eigen = NULL;
    };
  };
  if( !E_eigen
      && !( E_eigen = (double *)calloc( N, sizeof(double)) ) ) panic("Memory");

  cmp_householder_real( A, N, D, E_eigen, flag );
  cmp_QLShifts_real(    A, N, D, E_eigen, flag );
#ifdef DEBUG
  if( flag ){
    for( i = 0; i < N; i++ ){
      for( j = 0; j < N; j++ ){
	z = - save[ELEM(i,j,N)];
	for( k = 0; k < N; k++ )
	  z += A[ELEM(i,k,N)] * D[k] * conj(A[ELEM(j,k,N)]);
	if( cabs(z) > FLT_EPSILON ) panic("Complex diagonalization failed");
      };
    };
    free(save);
  };
#endif
};
/* ********************************************************************************* */
void cmp_Eigensystem( double complex *A, uint N, double *D ){
  cmp_Eigensystem_real( A, N, D, 1 );
};
/* ********************************************************************************* */
void cmp_Eigensystem_noeigenvectors( double complex *A, uint N, double *D ){
  if( N == 2 ){
    double x = creal(A[0]), y = creal(A[3]), z = 2.0 * cabs(A[1]);
#ifdef DEBUG
    if( fabs(cimag(A[0])) > FLT_EPSILON
	|| fabs(cimag(A[3])) > FLT_EPSILON
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
  cmp_Eigensystem_real( A, N, D, 0 );
};
/* ************************************************************************* */
static void cmp_recursive_determinant( double complex *A, uint N, uint start,
				       double complex *value, int *parity ){
  int i, j, idx_start, idx;
  double complex z;
  double help, max;
  if( !start ){
    *value = 0.0;
    *parity = 1;
  };
  idx_start = ELEM( start, 0, N );

  if( start == (N - 1) ){
    *value = A[idx_start + start];
  }else{
    max = cabs(A[idx_start + start]);
    j = start;
    for( i = start+1; i < N; i++ )
      if( (help = cabs(A[ELEM(i,start,N)])) > max ){
	max = help;
	j = i;
      };
    if( max < DBL_EPSILON ) panic2("Singular matrix: %.14f", max);
    if( j != start ){
      idx = ELEM( j, 0, N);
      for( i = start; i < N; i++ ){
	z = A[idx_start + i];
	A[idx_start + i] = A[idx + i];
	A[idx + i] = z;
      };
      *parity = -(*parity);
    };
    for( i = start + 1; i < N; i++ ){
      idx = ELEM( i, 0, N);
      for( j = start + 1; j < N; j++ ){
	z = A[idx + start] * A[idx_start + j];
	A[idx+j] *= A[idx_start + start];
	A[idx+j] -= z;
      };
    };
    cmp_recursive_determinant(A, N, start + 1, value, parity);
    j = N - start - 2;
    for( i = 0; i < j; i++ )   *value /= A[idx_start + start];
  };
};
/* ************************************************************************* */
double complex cmp_determinant( double complex *A, uint N ){
  double complex ret = 0.0;
  int parity = 1;
  if( N == 1 ){
    return A[0];
  }else if( N == 2 ){
    return A[0] * A[3] - A[2] * A[1];
  }else if( N == 3 ){
    return A[0] * A[4] * A[8] + A[1] * A[5] * A[6] + A[2] * A[3] * A[7] -
      A[6] * A[4] * A[2] - A[7] * A[5] * A[0] - A[8] * A[3] * A[1];
  };
  cmp_recursive_determinant( A, N, 0, &ret, &parity);
  return parity * ret;
};
/* ************************************************************************* */
static void dbl_recursive_determinant( double *A, uint N, uint start,
				       double *value, int *parity ){
  int i, j, idx_start, idx;
  double z;
  double help, max;
  if( !start ){
    *value = 0.0;
    *parity = 1;
  };
  idx_start = ELEM( start, 0, N );

  if( start == (N - 1) ){
    *value = A[idx_start + start];
  }else{
    max = cabs(A[idx_start + start]);
    j = start;
    for( i = start+1; i < N; i++ )
      if( (help = fabs(A[ELEM(i,start,N)])) > max ){
	max = help;
	j = i;
      };
    if( max < DBL_EPSILON ) panic2("Singular matrix: %.14f", max);
    if( j != start ){
      idx = ELEM( j, 0, N);
      for( i = start; i < N; i++ ){
	z = A[idx_start + i];
	A[idx_start + i] = A[idx + i];
	A[idx + i] = z;
      };
      *parity = -(*parity);
    };
    for( i = start + 1; i < N; i++ ){
      idx = ELEM( i, 0, N);
      for( j = start + 1; j < N; j++ ){
	z = A[idx + start] * A[idx_start + j];
	A[idx+j] *= A[idx_start + start];
	A[idx+j] -= z;
      };
    };
    dbl_recursive_determinant(A, N, start + 1, value, parity);
    j = N - start - 2;
    for( i = 0; i < j; i++ )   *value /= A[idx_start + start];
  };
};
/* ************************************************************************* */
double dbl_determinant( double *A, uint N ){
  double ret = 0.0;
  int parity = 1;
  if( N == 1 ){
    return A[0];
  }else if( N == 2 ){
    return A[0] * A[3] - A[2] * A[1];
  }else if( N == 3 ){
    return A[0] * A[4] * A[8] + A[1] * A[5] * A[6] + A[2] * A[3] * A[7] -
      A[6] * A[4] * A[2] - A[7] * A[5] * A[0] - A[8] * A[3] * A[1];
  };
  dbl_recursive_determinant( A, N, 0, &ret, &parity);
  return parity * ret;
};
