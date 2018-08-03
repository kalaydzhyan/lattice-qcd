#ifndef _SU2_UTILS_H_
#define _SU2_UTILS_H_

#include <math.h>
#include <complex.h>
#include <geometry.h>

typedef struct{
  /* This is just the first row of SU2 matrix, |alpha|^2+|beta|^2=1 */
  double complex alpha, beta;
} SU2;

/* returns just |z|^2 */
extern double cnorm( double complex z );
#ifdef RAVEN /* ;) */
extern double complex conj( double complex z );
#endif

/* General SU2 routines */
extern SU2    SU2_conj( SU2 v );
extern double SU2_det( SU2 v );
extern double SU2_normalize( SU2 *v );
extern SU2    SU2_mult( SU2 v1, SU2 v2 );

/* Generates random SU2 matrix g such that 1-1/2Tr g < delta
   (for delta = 2 you get indeed random g) */
extern SU2    SU2_random_matrix( double delta );

/* Stuff related to spin coherent states (SCS) */
extern double         SCS_solid_angle(  double complex z1, double complex z2, double complex z3 );
extern double         SCS_solid_angle4( double complex z0, double complex z1,
                    double complex z2, double complex z3 );
extern double complex SCS_eigenvector(  SU2 total, int sign );
extern double         SCS_step_phase( SU2 V, double complex z );
extern double complex SCS_step_shift( SU2 V, double complex z );
extern double         SCS_scalar_product( double complex z1, double complex z2 );

#include <diagonalize.h>
/* Routines to get SO(3) representation of SU(2) matrix and back;
   Correspondence is given by:
     O \in SO(3), U \in SU(2)
     O^{ab} = 1/2 Tr[ \sigma^a U \sigma^b U^+ ]
   The mapping SO(3)->SU(2) is two-fold and we always choose the
   representative which is closest to 1. */
extern int  SU2_from_SO3( double *matrix, SU2 *U );
extern void SO3_from_SU2( double *matrix, SU2 *U );
#endif
