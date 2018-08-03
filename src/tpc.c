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
#include "linal.h"

#ifdef SU2
#include "su2math.h"
#endif

#ifdef SU3
#include "su3math.h"
#endif

void printhelp(void) __attribute__ ((noreturn));

void printhelp(void)
{
	fprintf(stdout, "Determines the topological charge ...\n");
	
	fprintf(stdout, "Usage:\t-i evfilename\n");
	exit(EXIT_SUCCESS);
};

int main (int argc, char * argv[])
{
	int c;
	char *iname = NULL;
	t_op_ev evd;
	int i;
	int res;
	int vol;
	t_real csp = 1.0;
	
	while ((c = getopt(argc, argv, "i:")) != -1)
	{
		switch (c) {
		case 'i':
			iname = malloc(sizeof(char)*(strlen(optarg)+1));
			memset(iname, 0, sizeof(char)*(strlen(optarg)+1));
			strcpy(iname, optarg);
			break;
		case '?':
			if (isprint(optopt))
				fprintf(stderr, "Unknown option '-%c'.\n", optopt);
			else
				fprintf(stderr,
						"Unknown option character '\\x%x'.\n", optopt);
			printhelp();
		default:
			printhelp();
		}
	};
	if (!iname)
	 printhelp();

	evd.op = NULL;
	evd.nev = 0;
	evd.which[0] = 0;
	evd.which[1] = 0;
	evd.val_vec = 0;
	evd.tol = 0;
	evd.evals = NULL;
	evd.evecs = NULL;
	evd.data = 0;

	res = read_ev(iname, &evd, SM_SECTION);

	for(vol = 1, i = 0; i < evd.dim-1; i++)
	 vol *= evd.ls;
	vol *= evd.lt;
        vol *= evd.ndirac;
	vol *= evd.ncolors;

	// Topological charge
	
	int x, ic, id, Q = 0;
	t_real lambda;
	
	if(!((res & E_READ_CANTREAD) || (res & E_READ_CANTOPEN) || (res & E_READ_SECNOTHERE)))
         for(i = 0; i < evd.nev; i++)
         {
          lambda = 2.0*440.0*CCALL(cimag)(evd.evals[i])/(csp*(2.0 - CCALL(creal)(evd.evals[i])));
          // For zero mode - calculate the chirality to determine the topological charge
          t_complex r5 = 0.0 + I*0.0;
          for(x = 0; x<VOL; x++)
	   for(ic=0; ic<NCOLORS; ic++)
            for(id=0; id<NDIRAC; id++)
             r5 += CCALL(conj)(*(evd.evecs + i*vol + NCOLORS*NDIRAC*x + NDIRAC*ic + id))*(*(evd.evecs + i*vol + NCOLORS*NDIRAC*x + NDIRAC*ic + id))*(id < 2? 1.0 : -1.0);
          if(CCALL(creal)(r5)>0.9)
           Q += 1;
          if(CCALL(creal)(r5)<-0.9)
           Q -= 1;   
         };
         
	fprintf(stdout, "%i\n", Q);
	fflush(stdout);
	if (evd.evals != NULL) {
		free(evd.evals);
		evd.evals = NULL;
	}
	if (evd.evecs != NULL) {
		free(evd.evecs);
		evd.evecs = NULL;
	}

	return EXIT_SUCCESS;
}
