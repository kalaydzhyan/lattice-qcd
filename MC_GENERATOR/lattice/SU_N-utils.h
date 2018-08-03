#ifndef _SU_N_UTILS_H_
#define _SU_N_UTILS_H_

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <complex.h>


#include <geometry.h>
#include <newton-raphson.h>
#include <lapack.h>

#define SU_N_RANK  3

typedef struct{
  double complex U[SU_N_RANK][SU_N_RANK];
} SU_N;

/* ---------------------------------------------------------- */
/* Hermitian conjugation (any complex A) */
extern SU_N SU_N_conj( SU_N  A );

/* Complex determinant (any complex _non-degenerate_ A) */
extern double complex SU_N_determinant( SU_N A );

/* Projects arbitrary non-degenerate complex A to 'nearest'
   U(N) matrix (determinant is not corrected). Return status
   is 0 on success, otherwise is 1 (degenerate input matrix) */
extern int SU_N_projection_U_N( SU_N *A );

/* Projects non-degenerate complex A to 'nearest' SU(N) element.
   Return status is 0 on success, otherwise is 1 (degenerate input).
   First implementation does not care about [U(1)]^{N-1} subgroup
   and corrects determinant naively */
extern int SU_N_normalize_naive( SU_N *A );

/* Second implementation searches max_g Re Tr[g A^+] for g \in SU(N).
   Slow but honest use of newton-raphson method */
extern int SU_N_normalize(SU_N *A );

/* Matrix multiplication (any complex A, B) */
extern SU_N SU_N_mult( SU_N  A, SU_N B );

/* Matrix addition (any complex A, B) */
extern SU_N SU_N_sum(  SU_N  A, SU_N B );

/* Produces random SU(N) or U(N) matrix */
extern void SU_N_random( SU_N *A );
extern void SU_N_random_U_N( SU_N *A );

/* Returns simply  1/N Re Tr A  for arbitrary complex A */
extern double SU_N_norm( SU_N A );

#endif
