#include <SU_N-utils.h>

/* ************************************************************************* */
double complex SU_N_determinant( SU_N A ){
#if (SU_N_RANK == 2)
  return A.U[0][0] * A.U[1][1] - A.U[1][0] * A.U[0][1];
#elif (SU_N_RANK == 3)
  return A.U[0][0] * A.U[1][1] * A.U[2][2] + A.U[0][1] * A.U[1][2] * A.U[2][0] +
    A.U[0][2] * A.U[1][0] * A.U[2][1] - A.U[2][0] * A.U[1][1] * A.U[0][2] -
    A.U[2][1] * A.U[1][2] * A.U[0][0] - A.U[2][2] * A.U[1][0] * A.U[0][1];
#else
  return lapack_cmp_determinant( &(A.U[0][0]), SU_N_RANK );
#endif
};
/* ************************************************************************* */
SU_N SU_N_conj( SU_N A ){
  int i, j;
  double complex z;
  for( i = 0; i < SU_N_RANK; i++ ){
    A.U[i][i] = conj( A.U[i][i] );
    for( j = i+1; j < SU_N_RANK; j++ ){
      z = conj( A.U[i][j] );
      A.U[i][j] = conj( A.U[j][i] );
      A.U[j][i] = z;
    };
  };
  return A;
};
/* ************************************************************************* */
double SU_N_norm( SU_N A ){
  int i;
  double ret = 0.0;
  for( i = 0; i < SU_N_RANK; i++ ) ret += creal( A.U[i][i] );
  return ret/((double) SU_N_RANK);
};
/* ************************************************************************* */
SU_N SU_N_sum( SU_N A, SU_N B ){
  int i, j;
  for( i = 0; i < SU_N_RANK; i++ )
    for( j = 0; j < SU_N_RANK; j++ )
      A.U[i][j] += B.U[i][j];
  return A;
};
/* ************************************************************************* */
SU_N SU_N_mult( SU_N A, SU_N B ){
  int i, j, k;
  SU_N ret;
  bzero( &ret, sizeof(SU_N) );
  for( i = 0; i < SU_N_RANK; i++ )
    for( j = 0; j < SU_N_RANK; j++ )
      for( k = 0; k < SU_N_RANK; k++ )
	ret.U[i][j] += A.U[i][k] * B.U[k][j];
  return ret;
};
/* ************************************************************************* */
int SU_N_projection_U_N( SU_N *A ){
  uchar i;
  double D[ SU_N_RANK ];
  SU_N V;
  lapack_cmp_svd( &(A->U[0][0]), SU_N_RANK, SU_N_RANK, &D[0], &(V.U[0][0]) );
  for( i = 0; i < SU_N_RANK; i++ ){
#ifdef DEBUG
    if( D[i] < 0.0 ) panic2("Negative singular value: %.10f", D[i] );
#endif
    if( D[i] < FLT_EPSILON ){
      error2("Almost degenerate matrix, singular value %.10f", D[i] );
      return 1;
    };
  };
  *A = SU_N_mult( *A, V );
  return 0;
};
/* ************************************************************************** */
static void SU_N_random_real( SU_N *A, uchar flag ){
  int i, j;
  do{
    for( i = 0; i < SU_N_RANK; i++ )
      for( j = 0; j < SU_N_RANK; j++ )
	A->U[i][j] = (2.0 * RND() - 1.0) + I * (2.0 * RND() - 1.0);
    i = ( flag ) ? SU_N_normalize_naive(A) : SU_N_projection_U_N(A) ;
  }while(i);
};
/* ************************************************************************** */
void SU_N_random( SU_N *A ){
  SU_N_random_real( A, 1 );
};
/* ************************************************************************** */
void SU_N_random_U_N( SU_N *A ){
  SU_N_random_real( A, 0 );
};
/* ************************************************************************* */
int SU_N_normalize_naive( SU_N *A ){
  uchar i, j;
  double complex z;
  double phi;
  if( SU_N_projection_U_N( A ) ) return 1;
  z = SU_N_determinant( *A );
  phi = atan2( cimag(z), creal(z) )/((double) SU_N_RANK);
  z = cos(phi) - I * sin(phi);
  for( i = 0; i < SU_N_RANK; i++ )
    for( j = 0; j < SU_N_RANK; j++ )
      A->U[i][j] *= z;
  return 0;
};
/* ************************************************************************* */
/* Projection of general NxN complex matrix X  down to 'nearest'  SU(N) via
   maximization w.r.t. g \in SU(N) of F = Re Tr[ X^+ g]. If svd for X is
         X = U diag[d_i] V
   then we take g in the form    g = U diag[e^{i\phi_i}] V
   and extremize F w.r.t. \phi_i (preserving det g = 1)
*/
static int angles_extremum( double *x, void *data, double *val ){
  uchar i;
  double total = 0.0;
  double *d = (double *) data;
  for( i = 0; i < SU_N_RANK - 1; i++ ) total += x[i];
  total += d[SU_N_RANK];
  total = d[SU_N_RANK - 1] * sin( mod_pi(total) );
  for( i = 0; i < SU_N_RANK - 1; i++ ) val[i] = d[i] * sin(x[i]) + total;
  return 0;
};
/* ---------------------------------------------------------------------------- */
static int jacobian_OK( double *x, double *d ){
  int i, j;
  double A[ (SU_N_RANK - 1) * (SU_N_RANK - 1) ], D[SU_N_RANK - 1], total = 0.0;
  for( i = 0; i < SU_N_RANK - 1; i++ ) total += x[i];
  total += d[SU_N_RANK];
  total = - d[SU_N_RANK - 1] * cos( mod_pi(total) );
  for( i = 0; i < SU_N_RANK - 2; i++ ){
    for( j = i+1; j < SU_N_RANK - 1; j++ )
      A[ i * (SU_N_RANK - 1) + j ] = A[ j * (SU_N_RANK - 1) + i ] = total;
    A[ i * SU_N_RANK ] = total - d[i] * cos(x[i]);
  };
  A[SU_N_RANK * (SU_N_RANK - 2)] = total - d[SU_N_RANK - 2] * cos(x[SU_N_RANK - 2]);
  lapack_dbl_Eigensystem_noeigenvectors( A, SU_N_RANK - 1, D  );
  for( i = 0; i < SU_N_RANK - 1; i++ )
    if( D[i] > 0.0 ) return 0;
  return 1;
};
/* ---------------------------------------------------------------------------- */
int SU_N_normalize( SU_N *A ){
  int i, j;
  double D[ SU_N_RANK + 1];
  double phi[SU_N_RANK - 1], low[SU_N_RANK - 1], high[SU_N_RANK - 1];
  double complex z;
  SU_N V;
  z = SU_N_determinant( *A );
  D[SU_N_RANK] = atan2( cimag(z), creal(z) );
  lapack_cmp_svd( &(A->U[0][0]), SU_N_RANK, SU_N_RANK, &D[0], &(V.U[0][0]) );
  for( i = 0; i < SU_N_RANK; i++ )
    if( D[i] < FLT_EPSILON ) return 1;
  for( i = 0; i < SU_N_RANK - 1; i++ ){
    low[i] = - M_PI;
    high[i] =  M_PI;
  };
  do{
    for( i = 0; i < SU_N_RANK - 1; i++ ) phi[i] = M_PI * (2.0 * RND() - 1.0);
    i = newton_raphson( SU_N_RANK - 1, phi, low, high, (void *) D, &angles_extremum, NULL );
  }while( !i || !jacobian_OK( phi, D ) );
  for( i = 0; i < SU_N_RANK; i++ ){
    if( i == SU_N_RANK - 1 ){
      for( D[0] = 0.0, j = 0; j < SU_N_RANK - 1; j++ ) D[0] += phi[j];
      D[0] = mod_pi(D[SU_N_RANK] + D[0]);
      z = cos(D[0]) - I * sin(D[0]);
    }else{
      z = cos(phi[i]) + I * sin(phi[i]);
    };
    for( j = 0; j < SU_N_RANK; j++ ) V.U[i][j] *= z;
  };
  *A = SU_N_mult( *A, V );
  return 0;
};
