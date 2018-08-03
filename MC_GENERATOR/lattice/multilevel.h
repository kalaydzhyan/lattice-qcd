#ifndef _MULTILEVEL_H_
#define _MULTILEVEL_H_

#include <strings.h>
#include <MC-SU2.h>
#include <spinor.h>

/* Do not forget to change N_average/N_update when changing N_LEVELS */
#define           N_LEVELS   4
#define LOWEST_LEVEL_WIDTH   2

/* ------------------------------------------------- */
extern void ML_noise();

extern double *ML_correlator( uchar *size, char *filename );
extern double *ML_killing_correlator( uchar *size, char *filename, double gamma_s, double gamma_c );
#endif
