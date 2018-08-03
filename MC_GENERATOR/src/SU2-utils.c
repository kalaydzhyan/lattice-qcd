#include <SU2-utils.h>

#ifdef RAVEN
double complex conj( double complex z ){
  return creal(z) - I * cimag(z);
};
#endif

SU2 SU2_conj( SU2 v ){
  v.alpha = conj(v.alpha);
  v.beta = -v.beta;
  return v;
};
double SU2_det( SU2 v ){
  return creal(v.alpha * conj(v.alpha) + v.beta * conj(v.beta));
};
double SU2_normalize( SU2 *v ){
  double tmp = sqrt( SU2_det(*v) );
  v->alpha /= tmp;
  v->beta /= tmp;
  return tmp;
};
SU2 SU2_mult( SU2 v1, SU2 v2 ){
  SU2 ret = {
    v1.alpha * v2.alpha - v1.beta * conj(v2.beta),
    v1.alpha * v2.beta + v1.beta * conj(v2.alpha)
  };
  return ret;
};
double cnorm( double complex z ){
  double ret = cabs(z);
  return ret * ret;
};
SU2 SU2_random_matrix( double delta ){
  SU2 ret;
  uchar k;
  double x[3], alpha, help = 0;
  if( delta < 0 || delta > 2 ) delta = 2;
  for( k = 0; k < 3; k++ ){
    x[k] = 2.0*RND()-1.0;
    help += x[k]*x[k];
  };
  alpha = (2.0*RND()-1.0) * acos( 1.0 - delta );
  help = sin(alpha)/sqrt(help);
  ret.alpha = cos(alpha) + I * help * x[0];
  ret.beta = help * ( x[1] + I * x[2] );
  return ret;
};
/* =========================================== */
static inline double solid_angle_2( double complex z1, double complex z2 ){
  z1 = 1.0 + z1*conj(z2);
  return atan2( cimag(z1), creal(z1) );
};
double SCS_solid_angle( double complex z1, double complex z2, double complex z3 ){
  return mod_pi( solid_angle_2(z1,z2) + solid_angle_2(z2,z3) + solid_angle_2(z3,z1) );
};
double SCS_solid_angle4( double complex z0, double complex z1,
			 double complex z2, double complex z3 ){
  return mod_pi( SCS_solid_angle(z0,z1,z2) + SCS_solid_angle(z0,z2,z3) );
};
double complex SCS_eigenvector( SU2 total, int sign ){
  double complex ret = I;
  double a = cimag(total.alpha), b = cnorm(total.beta);
  if( sign > 0 )
    ret *= (a - sqrt(a * a + b) );
  else
    ret *= (a + sqrt(a * a + b) );
  return ret/conj( total.beta );
};
double SCS_step_phase( SU2  V, double complex z ){
  double complex zeta = V.alpha - conj( V.beta ) * z;
  return atan2( cimag(zeta), creal(zeta) );
};
double complex SCS_step_shift( SU2 V, double complex z ){
  return (V.beta + conj(V.alpha) * z) / (V.alpha - conj(V.beta) * z);
};
double SCS_scalar_product( double complex z1, double complex z2 ){
  return 1 - 2* cnorm(z1-z2)/((1+cnorm(z1)) * (1+cnorm(z2)));
};
/* =========================================== */
#define DIM     3
#define EPSILON 0.00000001
int SU2_from_SO3( double *matrix, SU2 *U ){ 
  int i,k;
  double help, tmp, eigen[DIM];
  double O[DIM * DIM];
  for( i = 0 ; i < DIM; i++ )
    for( k = i ; k < DIM; k++ )
      O[MATRIX_ELEMENT(i,k,DIM)] = O[MATRIX_ELEMENT(k,i,DIM)] =
	0.5 * ( matrix[MATRIX_ELEMENT(i,k,DIM)] + matrix[MATRIX_ELEMENT(k,i,DIM)] );
  /* find eigenvector with eigenvalue +/- 1 */
  if( Real_Symmetric_Diagonalize( &O[0], &eigen[0], DIM ) ) return -1;
  for( k = 0; k < DIM && fabs(eigen[k]-1.0) > EPSILON; k++ );
  if( k == DIM ){
    error("Eigenvalue 1 NOT found!");
    return -1;
  };
  /* sin 2alpha */
  help  = O[MATRIX_ELEMENT(0,k,DIM)] * O[MATRIX_ELEMENT(1,k,DIM)] * O[MATRIX_ELEMENT(2,k,DIM)] *
    ( matrix[MATRIX_ELEMENT(1,1,DIM)] - matrix[MATRIX_ELEMENT(0,0,DIM)] );
  help += O[MATRIX_ELEMENT(2,k,DIM)] * (
					O[MATRIX_ELEMENT(0,k,DIM)] * O[MATRIX_ELEMENT(0,k,DIM)] * matrix[MATRIX_ELEMENT(0,1,DIM)] - 
					O[MATRIX_ELEMENT(1,k,DIM)] * O[MATRIX_ELEMENT(1,k,DIM)] * matrix[MATRIX_ELEMENT(1,0,DIM)]
					);
  help += (1.0 - O[MATRIX_ELEMENT(2,k,DIM)] * O[MATRIX_ELEMENT(2,k,DIM)]) *
    (O[MATRIX_ELEMENT(1,k,DIM)] * matrix[MATRIX_ELEMENT(2,0,DIM)] - O[MATRIX_ELEMENT(0,k,DIM)] * matrix[MATRIX_ELEMENT(2,1,DIM)] );
  if( help*help > 1 ) return -1;
  /* cos 2alpha*/
  tmp  = O[MATRIX_ELEMENT(0,k,DIM)] * O[MATRIX_ELEMENT(0,k,DIM)] * matrix[MATRIX_ELEMENT(1,1,DIM)];
  tmp += O[MATRIX_ELEMENT(1,k,DIM)] * O[MATRIX_ELEMENT(1,k,DIM)] * matrix[MATRIX_ELEMENT(0,0,DIM)];
  tmp -= O[MATRIX_ELEMENT(0,k,DIM)] * O[MATRIX_ELEMENT(1,k,DIM)] * (matrix[MATRIX_ELEMENT(0,1,DIM)] + matrix[MATRIX_ELEMENT(1,0,DIM)]);
  if( tmp*tmp > 1 ) return -1;
  help = 0.5 * atan2( help, tmp );
  tmp = sin(help);
  (*U).alpha = cos(help) + I * (tmp * O[MATRIX_ELEMENT(2,k,DIM)] );
  (*U).beta = tmp * O[MATRIX_ELEMENT(1,k,DIM)] + I * ( tmp * O[MATRIX_ELEMENT(0,k,DIM)] );
  return 0;
};

void SO3_from_SU2( double *matrix, SU2 *U ){ 
  double x = creal( (*U).alpha );
  double y = cimag( (*U).alpha );
  double z = creal( (*U).beta );
  double w = cimag( (*U).beta );
  matrix[MATRIX_ELEMENT(0,0,DIM)] = x*x - y*y - z*z + w*w;
  matrix[MATRIX_ELEMENT(1,1,DIM)] = x*x - y*y + z*z - w*w;
  matrix[MATRIX_ELEMENT(2,2,DIM)] = x*x + y*y - z*z - w*w;
  matrix[MATRIX_ELEMENT(0,1,DIM)] = 2*(z*w+x*y);
  matrix[MATRIX_ELEMENT(1,0,DIM)] = 2*(z*w-x*y);
  matrix[MATRIX_ELEMENT(0,2,DIM)] = 2*(y*w-x*z);
  matrix[MATRIX_ELEMENT(2,0,DIM)] = 2*(y*w+x*z);
  matrix[MATRIX_ELEMENT(1,2,DIM)] = 2*(y*z+x*w);
  matrix[MATRIX_ELEMENT(2,1,DIM)] = 2*(y*z-x*w);
};
