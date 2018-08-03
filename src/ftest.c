#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include "defs.h"
#include "types.h"
#include "mt19937ar.h"
#include "lattice.h"
#include "arpack.h"
#include "dirac.h"
#include "panic.h"
#include "timing.h"
#include "linal.h"
#include "polacc.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef SU2
#include "su2math.h"
#endif

#ifdef SU3
#include "su3math.h"
#endif

int main (int argc, char * argv[]) /*{{{*/
{
#ifdef EFIELD
	int E = 1;
	int H = 0;
#endif

#ifdef EFIELD
	init_efield(E,H);
	t_real FE = 2*M_PI*(double)E/(double)(LS*LT);
	t_real FH = 2*M_PI*(double)H/(double)(LS*LS);
	fprintf(stdout, "Constant external electric and magnetic fields with |E| = %4.6lf, |H| = %4.6lf in lattice units are ON !!! \n",FE,FH);
	fprintf(stdout, "E is directed along X, H is directed along Z.\n");
	fprintf(stdout, "DIM: %i, LS: %i, LT: %i\n", DIM, LS, LT);
	fflush(stdout);
#endif

 lat_mov_init();

#ifdef EFIELD

 int s[DIM];
 int idx, idx_moved, mu, nu;
 t_real rplaq = 1;

 FORALL_4D(s)
 {
  idx = INDEX_X_4D(s);
  for(mu = 0; mu<3; mu++)
   for(nu = mu+1; nu<4; nu++)
   {
    t_complex plaq = ev[idx][mu]*conj(ev[idx][nu]);
    idx_moved = (*lat_mov)[idx][mu][FWD];
    plaq *= ev[idx_moved][nu];
    idx_moved = (*lat_mov)[idx][nu][FWD];
    plaq *= conj(ev[idx_moved][mu]);
    if(creal(plaq)!=rplaq)
    {
     fprintf(stdout,"NOT EQUAL: rplaq = %2.6lf, plaq = %2.6lf\n", rplaq, creal(plaq));
     rplaq = creal(plaq);
     fflush(stdout);
    };
   };
 };

#endif

#ifdef EFIELD
 fprintf(stdout,"\n\n");
 for(int q = 0; q < 30; q++)
 {
  s[0] = rand()%LS; s[1] = rand()%LS; s[2] = rand()%LS; s[3] = rand()%LT;
 
  fprintf(stdout,"x: [%i %i %i %i] \n", s[0], s[1], s[2], s[3]);
 
  for(mu = 0; mu<3; mu++)
   for(nu = mu+1; nu<4; nu++)
   {
    t_complex plaq = ev[idx][mu]*conj(ev[idx][nu]);
    idx_moved = (*lat_mov)[idx][mu][FWD];
    plaq *= ev[idx_moved][nu];
    idx_moved = (*lat_mov)[idx][nu][FWD];
    plaq *= conj(ev[idx_moved][mu]);
    fprintf(stdout,"[%i %i]: %2.6lf\n", mu, nu, creal(plaq)); 
    rplaq = creal(plaq);
    fflush(stdout);
   };
   
 };  
#endif

 return EXIT_SUCCESS;
}
