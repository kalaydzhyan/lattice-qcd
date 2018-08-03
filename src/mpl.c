/*{{{*/
/*!
 * \file mssc.c
 *
 * \brief Calculate spectrum, chirality, magnetization
 *
 *
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \date   Втр Сен 28 01:10:14 MSD 2004
 */
/*}}}*/

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
	fprintf(stdout, "Advanced Reader of EV files\n");
	fprintf(stdout, "Usage:\t-i evfilename -c csp -n mode_num -t time_slice\n");
	exit(EXIT_SUCCESS);
}

/*{{{*/
/*!
 * \brief  main routine
 *
 * \param argc
 * \param argv[]
 *
 * \return int 
 */
/*}}}*/

int main (int argc, char * argv[]) /*{{{*/
{
	int c;
	char *iname = NULL;
	t_op_ev evd;
	int i, k, swapped;
	int res;
	int vol;
	t_real csp = 1.0;
	int T = 0;
	int mode_num = 0;
	
	while ((c = getopt(argc, argv, "i:c:n:t:")) != -1) /*{{{*/
	{
		switch (c) {
		case 'i':
			iname = malloc(sizeof(char)*(strlen(optarg)+1));
			memset(iname, 0, sizeof(char)*(strlen(optarg)+1));
			strcpy(iname, optarg);
			break;
		case 'c':
			csp = (t_real)atof(optarg);
			break;
		case 'n':
			mode_num = atoi(optarg);
			break;
		case 't':
			T = atoi(optarg);
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
	}/*}}}*/
	if (!iname)
	 printhelp();
	//fprintf(stdout,"Input file :\t%s\n", iname);
	//fflush(stdout);

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
	int vol3d = vol; 
	vol *= evd.lt;
        vol *= evd.ndirac;
	vol *= evd.ncolors;

	if((res & E_READ_CANTREAD) || (res & E_READ_CANTOPEN) || (res & E_READ_SECNOTHERE))
	{
	 fprintf(stdout, "CANNOT READ INPUT FILE!!!\n");
	 fflush(stdout);
	 return EXIT_FAILURE;
	};
	
	// Sorting eigenvalues and eigenvectors
	
	t_real * lambda      = (t_real *)malloc(evd.nev*sizeof(t_real));
	int    * lambda_nums = (int    *)malloc(evd.nev*sizeof(int));
	
	
        for(i = 0; i < evd.nev; i++)
        {
         // Eigenvalues and their positions
         fprintf(stdout, "%i %4.6lf %4.6lf\n", i, (double)CCALL(cimag)(evd.evals[i]), (double)CCALL(creal)(evd.evals[i]));
         lambda[i]      = fabs(2.0*440.0*CCALL(cimag)(evd.evals[i])/(csp*(2.0 - CCALL(creal)(evd.evals[i]))));
         lambda_nums[i] = i;
        };
        
        // Sorting eigenvalues
        /*k = evd.nev;
        do{
         swapped = 0;
         for(i = 0; i<k; i++)
          if(lambda[i] >= lambda[i+1])
          {
           t_real tmpr = lambda[i];
           lambda[i] = lambda[i+1];
           lambda[i+1] = tmpr;
           int    tmpi = lambda_nums[i];
           lambda_nums[i] = lambda_nums[i+1];
           lambda_nums[i+1] = tmpi;
           swapped = 1;
          }; 
         k--;
        }while(swapped); 
        
        //for(i=0; i<evd.nev; i++)
        // fprintf(stdout, "%2.6lf %i\n", (double)lambda[i], lambda_nums[i]);
         
        int x, ic, id; 
                
        int evec_num = lambda_nums[mode_num];
      
        for(x=0; x<vol3d; x++)
        {    
         t_real rho = 0;
         for(ic=0; ic<NCOLORS; ic++)
          for(id=0; id<NDIRAC; id++)
          {
           rho += CCALL(conj)(*(evd.evecs + evec_num*vol + NCOLORS*NDIRAC*(vol3d*T + x) + NDIRAC*ic + id))*(*(evd.evecs + evec_num*vol + NCOLORS*NDIRAC*(vol3d*T + x) + NDIRAC*ic + id));
          };
         //fprintf(stdout, "%6.8lf\n", rho); 
         //fflush(stdout);
        };*/
                         
        free(lambda);
        free(lambda_nums); 
         
        
	if(evd.evals != NULL)
	{
	 free(evd.evals);
	 evd.evals = NULL;
	};
	
	if (evd.evecs != NULL)
	{
	 free(evd.evecs);
	 evd.evecs = NULL;
	};

	return EXIT_SUCCESS;
}/*}}}*/
