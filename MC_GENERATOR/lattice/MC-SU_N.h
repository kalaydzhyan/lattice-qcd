#ifndef _MC_SU_N_H_
#define _MC_SU_N_H_

#include <geometry.h>
#include <MC-SU2.h>
#include <SU_N-utils.h>

/* ============================================================
  Everywhere SU_N means just NxN complex matrix,
  which might be SU(N) or U(N) in applications.
  ============================================================= */

/* SU(N)/U(N) gauge fields with standart links/sites enumeration */
extern SU_N *F_N;

/* ------ Initialization of SU(N)/U(N) gauge data ------ */
/* Lists available read/write drivers */
extern void SU_N_list_drivers();
extern void  U_N_list_drivers();

/* Geometry and gauge fields initialization. If !filename we use 'hot start' */
extern int  SU_N_init( const char *filename, uchar T, uchar L, uchar D );
extern int   U_N_init( const char *filename, uchar T, uchar L, uchar D );

/* Dumps SU_N gauge data */
extern int  SU_N_dump( const char *filename );
extern int   U_N_dump( const char *filename );

/* Frees SU_N allocated memory (geometry untouched) */
extern void SU_N_free();
extern void  U_N_free();


/* ----------- Monte Carlo ---------------- */
/* SU(N): thermalization sweep with one full pseudo-heatbath plus
          SU(N) overrelaxation at each link;
   U(N) : thermalization sweep with Metropolis and U(N)
          overrelaxation at each link;
   All routines return averaged plaquette value */

extern double MC_U_N();
extern double MC_SU_N();

#define MC_SU_N_COOLING_DELTA   0.05
extern double MC_SU_N_cool();
extern double MC_U_N_cool();

extern double U_N_mono( char **current );

/* ------ Misc. routines ------ */
/* Mean plaquette for given SU(N)/U(N) configuration */
extern double SU_N_mean_plaq();
extern double  U_N_mean_plaq();

/* Random gauge transform */
extern void SU_N_random_gauge();
extern void  U_N_random_gauge();

/* Landau gauge fixing (not perfect !) */
#define LANDAU_GAUGE_PRECISION  0.00001
extern double SU_N_Landau_gauge();
extern double  U_N_Landau_gauge();

#endif
