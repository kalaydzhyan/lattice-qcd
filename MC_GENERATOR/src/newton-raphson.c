#include <newton-raphson.h>

extern double fmax(double x, double y);
extern double fmin(double x, double y);

#define LOGGING 0

/* Keep track of current problem dimension to reduce malloc/free overhead */
static uint NR_N = 0;

static uint   *NR_indx = NULL;
static double *NR_fjac = NULL;
static double *NR_g = NULL;
static double *NR_p = NULL;
static double *NR_xold = NULL;
static double *NR_fvec = NULL;
static double *NR_fvec_tmp = NULL;
static double *NR_high = NULL;
static double *NR_low = NULL;
static int     NR_error = 0;

/* ------------------------------------------------------------------------------------ */
/* Returns 1/2 \vec{F} \vec{F} at x, where \vec{F} is vector of user-supplied functions */
static double fminimum( double *x, void *data,
			int (*callback)( double *, void *, double *) ){
  uint i;
  double ret = 0.0;
  if( (NR_error = (*callback)( x, data, NR_fvec )) ) return 0.0;
  for( i = 0; i < NR_N; i++ )
    ret += NR_fvec[i] * NR_fvec[i];
  return 0.5 * ret;
};
/* ------------------------------------------------------------------------------------ */
/* Given an N-dimensional point x_old, the value of the function f_old and gradient g there,
   and a direction p, finds a new point x along p from x_old where the function 'fminimum' above
   has decreased "sufficiently". The new function value is returned in f. p is normally the Newton
   direction. Return status is zero on a normal exit. Otherwise, returns 1 indicating
   that x is too close to x_old. In a minimization algorithm, this usually signals convergence
   and can be ignored. However, in a zero-finding algorithm the caller should check whether the
   convergence is spurious. */
static int line_search( double *x_old, double f_old, double *g, double *p,
			double *x, double *f,
			void *data, int (*callback)( double *, void *, double * ) ){
  int i;
  double slope, temp, test, tmplam;
  double a, b, alam, alam2 = 0, alamin, disc, f2 = 0, rhs1, rhs2;
  for( a = 0.0, slope = 1.0, i = 0 ; i < NR_N; i++ ){
    b = fabs( p[i] );
    a = fmax( a, b );
    if( b > DBL_EPSILON ){
      test = x_old[i] + p[i];
      if( test > NR_high[i] ){
	temp = (NR_high[i] - x_old[i]) / p[i];
      }else if( test < NR_low[i] ){
	temp = (NR_low[i] - x_old[i]) / p[i];
      }else{
	temp = 1.0;
      };
      slope = fmin( slope, temp );
    };
  };
  if( fabs(slope) * a < FLT_EPSILON ) return 1;
#ifdef DEBUG
  if( slope < 0.0 ) panic2("Negative scale in line_search, %.10f", slope );
#endif
  for( i = 0; i < NR_N; i++ ) p[i] *= slope;

  for( slope = 0.0, i = 0; i < NR_N; i++ )  slope += g[i] * p[i];
  if( slope >= 0.0 )
    panic("Roundoff problem in line_search");

  for( test = 0.0, i = 0; i < NR_N; i++ ){  /* Compute lambda_min. */
    temp = fabs(p[i])/fmax( fabs(x_old[i]), 1.0 );
    if (temp > test) test = temp;
  };
  alamin = NR_TOLX / test;
  alam = 1.0;              /*Always try full Newton step first. */
  for( ; ; ){
    for( i = 0; i < NR_N; i++ ){
      x[i] = x_old[i] + alam * p[i];
      if( x[i] > NR_high[i] ) x[i] = NR_high[i];
      if( x[i] < NR_low[i]  ) x[i] = NR_low[i];
    };
    *f = fminimum( x, data, callback );
    if( NR_error ) return NR_error;

    if( alam < alamin ){
      /* Convergence on \delta x. For zero finding, the calling program should verify the convergence. */
      return 1;
    }else if( *f <= (f_old + NR_ALF * alam * slope) ){ /* Sufficient function decrease.*/
      return 0; 
    }else { /* Backtrack. */
      if( fabs(alam - 1.0) < DBL_EPSILON ){ /* First time. */
	tmplam = -slope/(2.0 * ( *f- f_old - slope )); 
      }else{ /* Subsequent backtracks. */
	rhs1 = *f - f_old - alam  * slope;
	rhs2 = f2 - f_old - alam2 * slope;
	a = ( rhs1/(alam*alam) - rhs2/(alam2*alam2) )/(alam-alam2);
	b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	if( fabs(a) < DBL_EPSILON ){
	  tmplam = -slope/(2.0*b);
	}else {
	  disc = b*b - 3.0 * a * slope;
	  if( disc < 0.0 ){
	    tmplam = 0.5 * alam;
	  }else if( b <= 0.0 ){
	    tmplam = ( -b + sqrt(disc) )/(3.0*a);
	  }else
	    tmplam = -slope/(b+sqrt(disc));
	};
	if (tmplam > 0.5 * alam )
	  tmplam = 0.5 * alam; /* lambda <= 0.5 * lambda_1 */
      };
    };
    alam2 = alam;
    f2 = *f;
    alam = fmax( tmplam, 0.1 * alam );  /* lambda >= 0.1 * lambda_1 */
  }; /* next iteration */
};
/* ------------------------------------------------------------------------------------ */
/* Computes forward-difference approximation to Jacobian. On input, x is the point at which
   the Jacobian is to be evaluated, fvec is the vector of function values at that point,
   'callback' are user-supplied routines that returns the functions at x.
   On output, df is the Jacobian array. */
static void fdjac( double *x, double *f, double *df,
		   void *data, int (*callback)( double *, void *, double * )){
  uint i, j;
  double h, temp;
  for( j = 0; j < NR_N; j++ ){
    temp = x[j];
    h = FLT_EPSILON * fabs(temp);
    if( h < DBL_EPSILON ) h = FLT_EPSILON;
    x[j] = temp + h;
    h = x[j]  - temp;
    if( (NR_error = (*callback)( x, data, NR_fvec_tmp )) ) return;
    x[j] = temp;
    for( i = 0; i < NR_N; i++ )
      df[i * NR_N + j] = ( NR_fvec_tmp[i] - f[i] )/h;
  };
};
/* ------------------------------------------------------------------------------------ */
/* For a given a matrix A replaces it by the LU decomposition of a rowwise permutation of itself.
   On output A is arranged as in equation (2.3.14) above.
   indx is an output vector that records the row permutation effected by the partial
   pivoting; Normal return status is \pm 1 which is sign det(A). For singular input matrix return
   status is 0. 'det' on output is optional determinant of input matrix - if you don't need it
   just pass NULL here. This routine is used in combination with lubksb to solve linear equations
   or invert a matrix. */
static uint N_LU = 0;
static double *vv_LU = NULL; /* stores the implicit scaling of each row. */

int ludcmp( double *A, uint n, uint *indx, double *det ){
  uint i, imax = 0, j, k, idx;
  double big, dum, sum, temp, dt = 1.0; /* No row interchanges yet. */
  if( N_LU != n ){
    N_LU = n;
    if( vv_LU ){ free(vv_LU); vv_LU = NULL; };
  };
  if( !vv_LU 
      && !(vv_LU = (double *)calloc( N_LU, sizeof(double))) ) panic("Memory");
  for( i = 0; i < N_LU; i++ ){ /* Loop over rows to get the implicit scaling information. */
    idx = i * N_LU;
    for( big = 0.0, j = 0; j < N_LU; j++ )
      if( (temp = fabs(A[idx + j])) > big) big = temp;
    if( big < FLT_EPSILON ){ /* Singular matrix */
      if( det ) *det = 0.0;
      return 0; 
    };
    vv_LU[i] = 1.0/big;
  };
  for( j = 0; j < N_LU; j++ ){ /*This is the loop over columns of Crout's method. */
    for( i = 0; i < j; i++ ){ /* This is equation (2.3.12) except for i = j. */
      idx = i * N_LU;
      sum = A[idx + j];
      for( k = 0; k < i; k++ ) sum -= A[idx + k] * A[k * N_LU + j];
      A[idx+j] = sum;
    };
    for( big = 0.0, i = j; i < N_LU; i++ ){ /* Initialize for the search for largest pivot element. */
      /* This is i = j of equation (2.3.12) and i = j+1, ..., N of equation (2.3.13). */
      idx = i * N_LU;
      sum = A[idx+j];
      for( k = 0; k < j; k++ )   sum -= A[idx + k] * A[k * N_LU + j];
      A[idx + j] = sum;
      if( (dum = vv_LU[i] * fabs(sum)) >= big ){
	/* Is the figure of merit for the pivot better than the best so far? */
	big = dum;
	imax = i;
      };
    };
    if( j != imax ){ /* Do we need to interchange rows? */
      idx = imax * N_LU;
      i   = j * N_LU;
      for( k = 0; k < N_LU; k++ ){ /* Yes, do so... */
	dum = A[idx + k];
	A[idx + k] = A[i + k];
	A[i + k] = dum;
      };
      dt *= -1.0;              /* ... and change parity */
      vv_LU[imax] = vv_LU[j]; /* Also interchange the scale factor. */
    };
    indx[j] = imax;
    if (j < N_LU - 1 ){ /* Now, finally, divide by the pivot element. */
      dum = 1.0/(A[j * (N_LU + 1)]);
      for( i = j+1; i < N_LU; i++ ) A[i * N_LU + j] *= dum;
    };
  }; /* Go back for the next column in the reduction. */
  for( i = 0; i < N_LU; i++ ) dt *= A[ i * (N_LU + 1) ];
  if( det ) *det = dt;
  if( fabs(dt) < FLT_EPSILON ) return 0;
  return (dt > 0.0) ? 1 : -1 ;
};
/* ------------------------------------------------------------------------------------ */
/* Solves the set of n linear equations A X = B. Here A is input, not as the matrix A but rather
   as its LU decomposition, determined by the routine ludcmp. indx is input
   as the permutation vector returned by ludcmp. b is input as the right-hand side vector
   B, and returns with the solution vector X. a, n, and indx are not modified by this routine
   and can be left in place for successive calls with different right-hand sides b. This routine
   takes into account the possibility that b will begin with many zero elements, so it is
   efficient for use in matrix inversion. */
void lubksb( double  *A, uint n, uint *indx, double *b ){
  int i, ip, j, idx, ii = -1;
  double sum;
  for( i = 0; i < n; i++ ){
    /* When ii is set to a positive value, it will become the index of the first nonvanishing
       element of b. We now do the forward substitution, equation (2.3.6). The only new wrinkle
       is to unscramble the permutation as we go. */
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    idx = i * n;
    if( ii >= 0 ){
      for( j = ii; j < i; j++ ) sum -= A[idx+j] * b[j];
    }else if( fabs(sum) > DBL_EPSILON ){
      /* A nonzero element was encountered, so from now on we will have to do the sums in the loop above. */
      ii=i;
    };
    b[i] = sum;
  };
  for( i = n - 1; i >= 0; i-- ){  /* Now we do the backsubstitution, equation (2.3.7). */
    idx = i * n;
    sum = b[i];
    for( j = i + 1; j < n; j++ ) sum -= A[idx+j] * b[j];
    b[i] = sum/A[idx+i];  /* Store a component of the solution vector X. */
  };
};
/* ------------------------------------------------------------------------------------ */
/* Given an initial guess x for a root in n dimensions, find the root by a globally convergent
   Newton-Raphson method. The vector of functions \vec{F} to be zeroed is given by the user-supplied
   callback. Return status is \pm 1 (sign of the Jacobian at root) on a normal return and 0 if the routine
   has presumably converged to a local minimum of the function \vec{F} \vec{F}. In this case
   try restarting from a different initial guess. */
int newton_raphson( uint n, double *x,
		    double *low, double *high,
		    void *data, int (*callback)( double *, void *, double *),
		    int *status ){
  int i, j, its, sign;
  double den, f, fold, sum, temp, test;
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if( NR_N != n ){
    NR_N = n;
    if( NR_indx ){ free(NR_indx);  NR_indx = NULL; };
    if( NR_fjac ){ free(NR_fjac);  NR_fjac = NULL; };
    if( NR_g ){ free(NR_g);  NR_g = NULL; };
    if( NR_p ){ free(NR_p);  NR_p = NULL; };
    if( NR_xold ){ free(NR_xold);  NR_xold = NULL; };
    if( NR_fvec ){ free(NR_fvec);  NR_fvec = NULL; };
    if( NR_fvec_tmp ){ free(NR_fvec_tmp);  NR_fvec_tmp = NULL; };
    if( NR_high ){ free(NR_high);  NR_high = NULL; };
    if( NR_low  ){ free(NR_low);   NR_low = NULL; };
  };
  if( !NR_indx && !(NR_indx = (uint *)calloc( NR_N, sizeof(uint))) ) panic("Memory");
  if( !NR_fjac && !(NR_fjac = (double *)calloc( NR_N*NR_N, sizeof(double))) ) panic("Memory");
  if( !NR_g    && !(NR_g = (double *)calloc( NR_N, sizeof(double))) ) panic("Memory");
  if( !NR_p    && !(NR_p = (double *)calloc( NR_N, sizeof(double))) ) panic("Memory");
  if( !NR_xold && !(NR_xold = (double *)calloc( NR_N, sizeof(double))) ) panic("Memory");
  if( !NR_fvec && !(NR_fvec = (double *)calloc( NR_N, sizeof(double))) ) panic("Memory");
  if( !NR_fvec_tmp && !(NR_fvec_tmp = (double *)calloc( NR_N, sizeof(double))) ) panic("Memory");
  if( !NR_high && !(NR_high = (double *)calloc( NR_N, sizeof(double))) ) panic("Memory");
  if( !NR_low  && !(NR_low  = (double *)calloc( NR_N, sizeof(double))) ) panic("Memory");
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#ifdef DEBUG
  for( i = 0; i < NR_N; i++ )
    if( low[i] >= high[i] || x[i] < low[i] || x[i] > high[i] )
      panic("Limits or initial guess are illegal");
#endif
  bcopy( low, NR_low, NR_N * sizeof(double));
  bcopy( high, NR_high, NR_N * sizeof(double));
  NR_error = 0;
  if( status ) *status = 0;
  f = fminimum( x, data, callback );     /* fvec is also computed by this call */
  if( NR_error ){
    if( status ) *status = NR_error;
    return 0;
  };
  for( test = 0.0, i = 0; i < NR_N; i++ ){ /* Test for initial guess being a root */
    den = fabs(NR_fvec[i]);
    if( den > test ) test = den;
  };
  if( test < NR_TOLF ){
    fdjac( x, NR_fvec, NR_fjac, data, callback );
    if( NR_error ){
      if( status ) *status = NR_error;
      return 0;
    };
    return ( lapack_dbl_determinant( NR_fjac, NR_N) > 0.0 ) ? 1 : -1;
  };

  for( its = 0; its < NR_MAXITS; its++ ){  /* Start of iteration loop. */
    /* If analytic Jacobian is available, you can replace fdjac with your own routine */
    fdjac( x, NR_fvec, NR_fjac, data, callback );
    if( NR_error ){
      if( status ) *status = NR_error;
      return 0;
    };
    for( i = 0; i < NR_N; i++ ){  /* Compute grad f for the line search. */
      for( sum = 0.0, j = 0; j < NR_N; j++ ) sum += NR_fjac[j * NR_N + i] * NR_fvec[j];
      NR_g[i] = sum;
    };
#if LOGGING
    printf("-------- ITERATION %d -------------\n", its );
    printf("Coordinates:\n");
    for( i = 0; i < NR_N; i++ ) printf("\t %f", x[i] );
    printf("\n");
    printf("Functions:\n");
    for( i = 0; i < NR_N; i++ ) printf("\t %f", NR_fvec[i] );
    printf("\n");
    printf("Gradient:\n");
    for( i = 0; i < NR_N; i++ ) printf("\t %f", NR_g[i] );
    printf("\n");
    printf("Jacobian: (det = %.14f)\n", lapack_dbl_determinant( NR_fjac, NR_N) );
    for( i = 0; i < NR_N; i++ ){
      for( j = 0; j < NR_N; j++ )
	printf("\t %f", NR_fjac[i*NR_N+j] );
      printf("\n");
    };
#endif
    for( i = 0; i < NR_N; i++ ) NR_xold[i] = x[i]; /* Store x and f. */
    fold = f; 
    for( i = 0; i < NR_N; i++ ) NR_p[i] = -NR_fvec[i]; /* Right-hand side for linear equations. */
    if( !(sign = ludcmp( NR_fjac, NR_N, NR_indx, NULL )) ){ /* Solve linear equations by LU decomposition. */
#if LOGGING
      printf("**** Degenerate Jacobian ****\n" );
#endif
      return 0;
    };
    lubksb( NR_fjac, NR_N, NR_indx, NR_p );
#if LOGGING
    printf("Newton direction:\n");
    for( i = 0; i < NR_N; i++ ) printf("\t %f", NR_p[i] );
    printf("\n");
#endif
    /* line_search returns new x and f. It also calculates fvec at the new x when it calls fminimum. */
    j = line_search( NR_xold, fold, NR_g, NR_p, x, &f, data, callback );
    if( NR_error ){
      if( status ) *status = NR_error;
      return 0;
    };
#if LOGGING
    printf("Coordinates [after line search]:\n");
    for( i = 0; i < NR_N; i++ ) printf("\t %f", x[i] );
    printf("\n");
#endif
    for( test = 0.0, i = 0; i < NR_N; i++ ){ /* Test for convergence on function values. */
      den = fabs( NR_fvec[i] );
      if( den > test ) test = den;
    };
    if( test < NR_TOLF ) return sign;

    if( j ){ /* Check for gradient of f zero, i.e., spurious convergence. */
      test = 0.0;
      den = fmax( f, 0.5 * NR_N );
      for( i = 0; i < NR_N; i++ ){
	temp = fabs(NR_g[i]) * fmax( fabs(x[i]), 1.0 )/den;
	if (temp > test) test = temp;
      };
      if( test < NR_TOLMIN ){
#if LOGGING
	printf("'Spurious convergence', gradient %f\n", test );
#endif
	return 0;
      };
    };
    for( test = 0.0, i = 0; i < NR_N; i++ ){ /* Test for convergence on \delta x. */
      temp = fabs( x[i] - NR_xold[i] ) / fmax( fabs(x[i]), 1.0 );
      if( temp > test ) test = temp;
    };
    if( test < NR_TOLX ){
#if LOGGING
      printf("'Spurious convergence', delta x  %f\n", test );
#endif
      return 0;
    };
  };
  return 0;
};
/* ------------------------------------------------------------------------------------ */
/* Golden section search of maximum (sign=1) / minimum (sign=-1) of given function
   in the interval [a,c]. Triple a,b,c MUST bracket extremum. Returns function's
   extremal value, 'double *t' points to extremum location on exit. */
double find_extremum( double a, double b, double c, double *t, int sign,
		      void *data, double (*callback)( double , void * ) ){
  double f1, f2, x0, x1, x2, x3;
  x0 = a;
  x3 = c;
  if( fabs(c - b) > fabs(b - a) ){
    x1 = b; x2 = b + NR_C * (c - b);
  }else{
    x2 = b; x1 = b - NR_C * (b - a);
  };
  f1 = sign * (*callback)( x1, data );
  f2 = sign * (*callback)( x2, data );
  while( fabs(x3 - x0) > NR_PRES * (fabs(x1) + fabs(x2)) ){
    if( f1 < f2 ){
      x0 = x1; x1 = x2; x2 = NR_R * x1 + NR_C * x3;
      f1 = f2; 
      f2 = sign * (*callback)( x2, data );
    }else{
      x3 = x2; x2 = x1; x1 = NR_R * x2 + NR_C * x0;
      f2 = f1;
      f1 = sign * (*callback)( x1, data );
    };
  };
  if( f2 < f1 ){
    *t = x1;
    return f1;
  };
  *t = x2;
  return f2;
};
/* ------------------------------------------------------------------------------------ */
double find_root( double start, double F_start, double end,   double F_end, double *t,
		  void *data, double (*callback)( double, void * ) ){
  double F = 0.0;
  if( F_start * F_end > 0.0 )
    panic3("Root is not located in [%f;%f]", start, end );
  if( start > end ){
    F =   start;   start =   end;   end = F;
    F = F_start; F_start = F_end; F_end = F;
  };
  do{
    *t = 0.5 * ( start + end );
    F = (*callback)( *t , data );
    if( F * F_start > 0.0 )   start = *t;
    else                      end   = *t;
  }while( end - start > NR_PRES );
  return F;
};
