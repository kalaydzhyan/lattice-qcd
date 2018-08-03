#ifndef _MC_KILLING_H_
#define _MC_KILLING_H_

#include <spinor.h>

/* Performs MC_KILLING_SWEEPS Monte Carlo sweeps with Metropolis algorithm and action:
   D=3:  gamma_s \sum_{sites}   |q_{site}|    + gamma_c \sum_{cubes} |q_{cube}|

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   YYYYYYYYYYYYYYYYYYYY Should be properly commented YYYYYYYYYYYYYYYYYYYY
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   D=4: gamma_s \sum_{diamond} |q_{diamond}| + gamma_c \sum_{cubes} |q_{cube}|
   Gauge action is implicitely assumed. Half of MC_KILLING_TRY are overrelaxation attempts.*/

#define MC_KILLING_N_TRY     6
extern int MC_killing( uint N_sweeps, double gamma_s, double gamma_c );
#endif
