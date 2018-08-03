/*{{{*/
/*!
 * \file wd_prop.c
 *
 * \brief Calculate the propagator for the Wilson Fermions
 *
 *
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \author Pavel Buividovich, email: gbuividovich@gmail.com (implemented background magnetic field, chemical potential, SU(3) gauge group in 2008 - 2009)
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
#include "timing.h"
#include "linal.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef SU2
#include "su2math.h"
#endif

#ifdef SU3
#include "su3math.h"
#endif

#define IMAX 5000

void printhelp(void) __attribute__ ((noreturn));

/*{{{*/
/*!
 * \brief Print help/usage in case of wrong option
 *
 * \param void
 *
 * \return void 
 *//*}}}*/
void printhelp(void) /*{{{*/
{
	fprintf(stdout, "Spectrum of Wilson Dirac Operator\n");
	fprintf(stdout, "Usage:\n");
	fprintf(stdout, "[-f]\t -- free operator\n");
	fprintf(stdout, "[-I]\t -- instanton\n");
	fprintf(stdout, "[-a]\t -- aperiodic boundary conditions\n");
	fprintf(stdout, "[-i]\t -- info of gauge configuration file\n");
#ifdef EFIELD
	fprintf(stdout, "-E num\t -- external electric field. Default is zero.\n");
	fprintf(stdout, "-H num\t -- external magnetic field. Default is zero.\n");
#endif
	fprintf(stdout, "-r rho\t -- mass parameter rho. rho > 0\n");
	fprintf(stdout, "[-k num]\t -- anisotropy factor, the ratio of kappa_t to kappa_s\n");
	fprintf(stdout, "[-s]\t -- rho = -|rho|\n");
	fprintf(stdout, " -o \t -- outputfilename\n");
	fprintf(stdout, "[-t]\t -- tolerance of BiCGStab. [DEFAULT: %f]\n", 0.00001);

	exit(EXIT_SUCCESS);
}/*}}}*/

/*{{{*/
/*!
 * \brief  main routine
 *
 * \param argc
 * \param argv[]
 *
 * \return int 
 *//*}}}*/
int main (int argc, char * argv[]) /*{{{*/
{
	t_gauge_vector * gv;
	int c;
	int i, ret;
	char *iname = NULL;
	char *oname = NULL;
	t_real tol = 0.00001; //TODO: define
	t_real kappa = -1.0, rho = -1.0;
	t_op_ev evd;
	t_real rhosign = -1.0;
	int freewd = 0;
	t_real action;
#ifdef EFIELD
	int E = 0, H = 0;
#endif
    anisotropy_factor = 1.0;
	int instanton = 0;
	t_real ic[DIM] = {0.0, 0.0, 0.0, 0.0};
	int aperiodic = 0;
	int a, id, x0;
	t_real dist;
    
  	while ((c = getopt(argc, argv, "ReQVr:E:H:o:sfi:t:l:m:d:Iak:")) != -1) /*{{{*/
	{
		switch (c) {
	        case 'k':
    		        anisotropy_factor = atof(optarg);
	        break; 
		case 'a':
			aperiodic = 1;
			break;
		case 'I':
			instanton = 1;
			break;
		case 'o':
			oname = malloc(sizeof(char)*(strlen(optarg)+1));
			memset(oname, 0, sizeof(char)*(strlen(optarg)+1));
			strcpy(oname, optarg);
			break;
		case 'r':
			rho = atof(optarg);
			break;
#ifdef EFIELD
        case 'E':
            E = atoi(optarg);
            break;
        case 'H':
            H = atoi(optarg);
            break;
#endif
		case 's':
			rhosign = 1.0;
			break;
		case 'f':
			freewd = 1;
			break;
		case 'i':
			iname = malloc(sizeof(char)*strlen(optarg));
			strcpy(iname, optarg);
			break;
		case 't':
			tol = atof(optarg);
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
	if (!(instanton || freewd) && (iname == NULL))
		printhelp();
	else if ((instanton || freewd) && iname != NULL)
		printhelp();
	if (!oname)
		printhelp();
	if(rho < 0.0)
		printhelp();
	kappa = rho2kappa(rhosign*rho);

#ifdef _OPENMP
    fprintf(stdout,"Running in parallel mode with %i threads \n", omp_get_max_threads());
#endif

	fprintf(stdout,"Input file :\t%s\n", (freewd ? "Free" : iname));
	fprintf(stdout,"Output file:\t%s\n", oname);
	fprintf(stdout,"Size       :\t%ix%ix%ix%i\n", LS, LS, LS, LT); 
	fprintf(stdout,"Group      :\tSU(%i)\n", NCOLORS);
	fprintf(stdout,"Tolerance  :\t%12.10f\n", tol);
	fprintf(stdout,"rho        :\t%12.10f\n", rho);
	fprintf(stdout,"kappa      :\t%12.10f\n", kappa);
	fprintf(stdout,"anisotropy :\t%12.10f\n", anisotropy_factor);
	fprintf(stdout,"operator   :\t%s\n", ("Wilson Dirac"));
	gv = lat_gauge_create();
	lat_mov_init();
	if(freewd)
	 lat_gauge_identity(gv);
	else
	 if(instanton)
	 {
	  ic[0] = LS/2.0+0.5;
	  ic[1] = LS/2.0+0.5;
	  ic[2]= LS/2.0+0.5;
	  ic[3] = LT/2.0+0.5;
	  lat_gauge_instanton(gv, 2.0, ic);
	 }
	 else
	 {
	  lat_gauge_load(gv, iname);
	  action = lat_gauge_action(gv);
	  fprintf(stdout,"\n\n\t action     :\t%12.10f\n", action);
	  fflush(stdout);
	 };
	
	lat_gauge_check_unit(gv); 
	if(aperiodic) 
	 lat_gauge_aperiodic(gv);
	 
#ifdef EFIELD
	init_efield(E,H);
	t_real FE = 2*M_PI*(double)E/(double)(LS*LT);
	t_real FH = 2*M_PI*(double)H/(double)(LS*LS);
	fprintf(stdout,"Constant external electric and magnetic fields with |E| = %4.6lf, |H| = %4.6lf in lattice units are ON !!! \n",FE,FH);
	fprintf(stdout,"E is directed along X, H is directed along Z.\n");
	fflush(stdout);
#endif
	
	evd.val_vec = 1;
	evd.tol = tol;
	evd.data = &kappa;
#ifdef SCALAR
         evd.op = &(scalar);
#else
	 evd.op = &(wdirac);
#endif
	
#ifdef TIMING
	init_time();
#endif

 t_cds_vector *source, *solution;
 if(posix_memalign((void*)&source,      16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for source (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&solution,    16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for solution      (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void *)&evd.evecs,        16, NCOLORS*NDIRAC*sizeof(t_cds_vector)))
   panic("%s: can't allocate memory for evd.evecs        (%i Kb)", __func__, NCOLORS*NDIRAC*sizeof(t_cds_vector)/1024);

 for(a=0; a<NCOLORS; a++)
  for(id=0; id<NDIRAC; id++)
  {
   i = NDIRAC*a + id;
   fprintf(stdout, "\n\n\t Inverting the Wilson Dirac operator for source %i ... \n\n", i+1);
   set_zero((t_complex *)source);
   // preparing a gaussian source
   x0 = 0;
   *((t_complex *)source + NDIRAC*NCOLORS*x0 + i) = 1.0 + I*0.0;
   dist = vnorm((t_complex *)source);
   ax((t_complex *)source, 1.0/dist);

#ifdef TIMING
   init_time();
#endif
   ret = bicgstab(&evd, gv, source, solution, tol, IMAX);
   if(ret!=0)
   {
    fprintf(stdout, "\n\n\t BiCGStab has not reached the precision %4.4E in allowed number %i of iterations at source %i. Program terminates.\n\n", tol, IMAX, NDIRAC*a + id);
    fprintf(stderr, "\n\n\t BiCGStab has not reached the precision %4.4E in allowed number %i of iterations at source %i. Program terminates.\n\n", tol, IMAX, NDIRAC*a + id);
    return EXIT_FAILURE;
   };
#ifdef TIMING
   print_time();
#endif
   memcpy(&(evd.evecs[i*VOL*NCOLORS*NDIRAC]), (t_complex *)solution, sizeof(t_cds_vector));
  };

  evd.nev   = NCOLORS*NDIRAC;
  evd.evals = (t_complex *)malloc(NCOLORS*NDIRAC*sizeof(t_complex));
  for(i=0; i<NCOLORS*NDIRAC; i++)
   evd.evals[i] = (double)i;
  write_ev(oname, "shadow", &evd, SM_SECTION);

	
#ifdef TIMING
	print_time();
#endif
	
	write_ev(oname, "shadow", &evd, SM_SECTION);  // if LR == 1, this section contains snev Right eigenvalues: D R = Lambda R
    // Printing out the eigenvalues - for control 
    

    if (evd.evals != NULL)
    { 
     free(evd.evals);
     evd.evals = NULL;
    };
    
    if (evd.evecs != NULL)
    {
     free(evd.evecs);
     evd.evecs = NULL;
    };
    
    lat_gauge_destroy(gv);

    return EXIT_SUCCESS;
}/*}}}*/
