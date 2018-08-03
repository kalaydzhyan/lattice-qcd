/*{{{*/
/*!
 * \file wd.c
 *
 * \brief Calculate spectrum of Wilson Fermions
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

#define DEF_POLACC_DEG	    (30)
#define DEF_ARNOLDI_TOL	    (1e-10)
#define DEF_ARNOLDI_LNEV    (1)
#define DEF_NCV_FACTOR_NOM  (5)
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
	fprintf(stdout, "[-h]\t -- calculate spectrum of H = gamma_5 Q, otherwise Q\n");
	fprintf(stdout, "[-F]\t -- NCV_FACTOR_NOM, parameter of the Arnoldi algorythm, default is 5\n");
#ifdef EFIELD
	fprintf(stdout, "-E num\t -- external electric field. Default is zero.\n");
	fprintf(stdout, "-H num\t -- external magnetic field. Default is zero.\n");
#endif
#ifdef MU
	fprintf(stdout, "-M mu\t -- chemical potential mu, default is zero.\n");
	fprintf(stdout, "[-R]\t -- calculate left and right eigenvectors\n");
	fprintf(stdout, "[-e]\t -- calculate eigenvalues with smallest real parts, not smallest absolute value\n");
	fprintf(stdout, "[-Q]\t -- calculate the spectrum of the squared \"Hermitean\" (at mu = 0) Wilson Dirac operator\n");
#endif
   	fprintf(stdout, "[-V]\t -- print out the eigenvalues found.\n");
	fprintf(stdout, "-r rho\t -- mass parameter rho. rho > 0\n");
	fprintf(stdout, "[-k num]\t -- anisotropy factor, the ratio of kappa_t to kappa_s\n");
	fprintf(stdout, "[-s]\t -- rho = -|rho|\n");
	fprintf(stdout, " -o \t -- outputfilename\n");
	fprintf(stdout, " -S \t -- # of smallest eigenvalues\n");
	fprintf(stdout, "[-L]\t -- # of largest eigenvalues. [DEFAULT: %i]\n", DEF_ARNOLDI_LNEV);
	fprintf(stdout, "[-t]\t -- tolerance of Arnoldi process. [DEFAULT: %f]\n", DEF_ARNOLDI_TOL);
	fprintf(stdout, "[-p]\t -- turn on polynomial acceleration\n");
	fprintf(stdout, "Polynomial acceleration parameters (if -p parameter is specified):\n");
	fprintf(stdout, "-m \t -- polynomial acceleration m\n");
	fprintf(stdout, "[-l]\t -- polynomial acceleration l\n");
	fprintf(stdout, "[-d]\t -- polynomial acceleration degree. [DEFAULT: %i]\n", DEF_POLACC_DEG);

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
	int i;
	char *iname = NULL;
	char *oname = NULL;
	int snev = -1, lnev = DEF_ARNOLDI_LNEV;
	t_real tol = DEF_ARNOLDI_TOL;
	t_real kappa = -1.0, rho = -1.0;
	t_op_ev evd;
	t_real rhosign = -1.0;
	int hermitianQ = 0; /* calculate either Q or H */
	int freewd = 0;
	int polacc = 0;
	double polacc_l = -1.0;
	double polacc_m = -1.0;
	int polacc_deg = DEF_POLACC_DEG;
	t_polacc_data polacc_data;
	t_real action;
#ifdef EFIELD
	int E = 0, H = 0;
#endif
    anisotropy_factor = 1.0;
#ifdef MU
	Mu = 0.0;
	int LR = 0; // LR = 1 - calculate left and right eigenvectors, LR = 0 - no changes
	int SR = 0; // SR = 1 - calculate eigvecs with smallest Real part of eigvals, 0 = no changes
	int SQ = 0; // SQ = 1 - calculate eigvecs of the square of "Hermitean" (at mu = 0) Wilson Dirac, 0 = no changes
#endif
   	int PV = 0; // printout eigenvalues at the end of calculation
	int instanton = 0;
	t_real ic[DIM] = {0.0, 0.0, 0.0, 0.0};
	int aperiodic = 0;
    int ncv_factor_nom = DEF_NCV_FACTOR_NOM, ncv_factor_denom = 2;
    
    	while ((c = getopt(argc, argv, "pReQVr:E:H:M:o:S:L:shfi:t:l:m:d:IaF:k:")) != -1) /*{{{*/
	{
		switch (c) {
                case 'F':
                        ncv_factor_nom = atoi(optarg);
                        if (ncv_factor_nom <= 0)
                            ncv_factor_nom = DEF_NCV_FACTOR_NOM;
                        break;
	        case 'k':
    		        anisotropy_factor = atof(optarg);
	        break; 
		case 'a':
			aperiodic = 1;
			break;
		case 'I':
			instanton = 1;
			break;
		case 'p':
			polacc = 1;
			break;
		case 'o':
			oname = malloc(sizeof(char)*(strlen(optarg)+1));
			memset(oname, 0, sizeof(char)*(strlen(optarg)+1));
			strcpy(oname, optarg);
			break;
		case 'S':
			snev = atoi(optarg);
			break;
		case 'L':
			lnev = atoi(optarg);
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
#ifdef MU
        case 'M':
            Mu = atof(optarg);
            break;
        case 'R':
            LR = 1;
            break; 
        case 'e':
            SR = 1;
            break; 
        case 'Q':
            SQ = 1;
            break;
#endif
		case 'V':
		        PV = 1;
			break;
		case 'l':
			polacc_l = atof(optarg);
			break;
		case 'm':
			polacc_m = atof(optarg);
			break;
		case 's':
			rhosign = 1.0;
			break;
		case 'h':
			hermitianQ = 1;
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
		case 'd':
			polacc_deg = atoi(optarg);
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
#ifdef MU
    if(polacc==1)
    {
     polacc = 0;
     fprintf(stdout, "\n\nPolacc does not work for finite mu, turned OFF automatically!!!\n\n");
     fflush(stdout);
    };
#endif
	if (!(instanton || freewd) && (iname == NULL))
		printhelp();
	else if ((instanton || freewd) && iname != NULL)
		printhelp();
	if (!oname)
		printhelp();
	if (snev < 0 || rho < 0.0)
		printhelp();
	kappa = rho2kappa(rhosign*rho);
	if (polacc && polacc_l < 0.0)
		polacc_l = 1.0 + 8.0*kappa;
	if (polacc && polacc_m < 0.0) {
		fprintf(stdout, "Polynomial acceleration is ON, but polacc_m is not set\n");
		printhelp();
	}
	if (polacc && !hermitianQ) {
		fprintf(stdout, "Polynomial acceleration is ON, but non-hermitian op is set\n");
		printhelp();
	}

#ifdef _OPENMP
    fprintf(stdout,"Running in parallel mode with %i threads \n", omp_get_max_threads());
#endif

	fprintf(stdout,"Input file :\t%s\n", (freewd ? "Free" : iname));
	fprintf(stdout,"Output file:\t%s\n", oname);
	fprintf(stdout,"Size       :\t%ix%ix%ix%i\n", LS, LS, LS, LT); 
	fprintf(stdout,"Group      :\tSU(%i)\n", NCOLORS);
	fprintf(stdout,"Tolerance  :\t%12.10f\n", tol);
	fprintf(stdout,"SM NEV     :\t%i\n", snev);
#ifdef MU
    if(!LR)
    {
     fprintf(stdout,"LM NEV     :\t%i\n", lnev);
    }
    else
    {
     fprintf(stdout, "Calculating %i LEFT and RIGHT eigenvectors\n", snev);    
    }; 
#else
	fprintf(stdout,"LM NEV     :\t%i\n", lnev);
#endif
	fprintf(stdout,"rho        :\t%12.10f\n", rho);
	fprintf(stdout,"kappa      :\t%12.10f\n", kappa);
	fprintf(stdout,"anisotropy :\t%12.10f\n", anisotropy_factor);
#ifdef MU
	fprintf(stdout,"!!!!!!!!!! Running with chemical potential !!!!!!!!!!!!!!!!!!\n");
	fprintf(stdout,"mu         :\t%12.10f\n", Mu);
#endif
#ifdef MU
    fprintf(stdout,"operator   :\t%s %s\n", (hermitianQ ? "Hermitian Wilson Dirac" : "Wilson Dirac"), (SQ ? "SQUARED" : ""));
    if(SR==1)
     fprintf(stdout,"Calculating eigenvalues with the smallest REAL part!!!\n"); 
#else
	fprintf(stdout,"operator   :\t%s\n", (hermitianQ ? "Hermitian Wilson Dirac" : "Wilson Dirac"));
#endif
	gv = lat_gauge_create();
	lat_mov_init();
	if(freewd)
	 if (polacc)
	 {
	  fprintf(stdout, "Polynomial acceleration can't be used for degenerate spectrum, i.e. free spectrum\n");
	  printhelp();
	 }
	 else
	 {
	  lat_gauge_identity(gv);
	 }
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
	  fprintf(stderr,"action     :\t%12.10f\n", action);
	  fflush(stderr);
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
	
	fprintf(stdout, "Starting Arnoldi iterations...\n");
	fflush(stdout);

	evd.nev = snev;
	evd.val_vec = 1;
	evd.tol = tol;
	evd.evals = NULL;
	evd.evecs = NULL;
	if (polacc) {
		evd.which[0] = 'L';
		evd.which[1] = 'M';
		fprintf(stdout, "P_m : %f\n", polacc_m);
		fprintf(stdout, "P_l : %f\n", polacc_l);
		fprintf(stdout, "P_deg : %i\n", polacc_deg);
		build_mapping(polacc_m, polacc_l, &polacc_data);
		polacc_data.op = &(h_wdirac);
		polacc_data.deg = polacc_deg;
		polacc_data.data = &kappa;
		evd.op = &pol_acc;
		evd.data = &polacc_data;
	} else {
		evd.which[0] = 'S';
		evd.which[1] = 'M';
#ifdef MU
       if(SR==1)
       {
        evd.which[0] = 'S';
	    evd.which[1] = 'R';        
       };
#endif
		evd.data = &kappa;
		if (!hermitianQ)
#ifdef SCALAR
			evd.op = &(scalar);
#else
			evd.op = &(wdirac);
#endif
		else
		{
		 evd.op = &(h_wdirac);
#ifdef MU
         if(SQ==1)
          evd.op = &(h_wdirac_sq);
#endif
    		};
	}
	
#ifdef TIMING
	init_time();
#endif
#ifdef MU
    fprintf(stdout, "\n\nCalculating %s eigenvalues/vectors ... \n\n", (LR ? "RIGHT" : "SMALLEST"));
#else
    fprintf(stdout, "\n\nCalculating eigenvalues/vectors ... \n\n");
#endif
	arnoldi(gv, &evd, ncv_factor_nom, ncv_factor_denom);
#ifdef TIMING
	print_time();
#endif
	if (polacc) {
		orig_eval_hw(gv, &evd);
	}
	write_ev(oname, "shadow", &evd, SM_SECTION);  // if LR == 1, this section contains snev Right eigenvalues: D R = Lambda R
    // Printing out the eigenvalues - for control 
    
    if(PV)
    {
#ifdef MU
     fprintf(stdout, "\n\n %s eigenvalues: \n\n", (LR ? "RIGHT" : "SMALLEST"));
#else
     fprintf(stdout, "\n\n Smallest eigenvalues: \n\n");
#endif
     for(i=0; i<evd.nev; i++)
      fprintf(stdout, "%i: %8.8E %8.8E\n", i, (double)CCALL(creal)(evd.evals[i]), (double)CCALL(cimag)(evd.evals[i]));
    };    
    
    /* largest eigenvalues/eigenfunctions, or right eigenvalues/eigenfunctions */
    evd.nev = lnev;
    evd.which[0] = 'L';
    evd.which[1] = 'M';
#ifdef MU
    if(SR==1)
    {
     evd.which[0] = 'L';
     evd.which[1] = 'R';
    };
#endif
    evd.val_vec = 1;
    evd.tol = 1e-3;
    evd.data = &kappa;
#ifdef MU
    if(LR==1)
    {
     Mu = -1.0*Mu; // Hermitean conjugate operator        
     evd.nev = snev;
     evd.val_vec = 1; // Store eigenvectors
     evd.tol = tol;
     evd.evals = NULL;
     evd.evecs = NULL;
     evd.data = &kappa;
     evd.which[0] = 'S';
     evd.which[1] = 'M';
     if(SR==1)
     {
      evd.which[0] = 'S';
      evd.which[1] = 'R';        
     };        
    };
#endif
    if (!hermitianQ)
#ifdef SCALAR
     evd.op = &(scalar);
#else
     evd.op = &(wdirac);
#endif
    else
    {
     evd.op = &(h_wdirac);
#ifdef MU
     if(SQ==1)
      evd.op = &(h_wdirac_sq);
#endif
    };
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
#ifdef MU
    fprintf(stdout, "\n\nCalculating %s eigenvalues/vectors ... \n\n", (LR ? "LEFT" : "LARGEST"));
#endif
	arnoldi(gv, &evd, ncv_factor_nom, ncv_factor_denom);
	write_ev(oname, "shadow", &evd, LM_SECTION); // if LR == 1, this section contains snev LEFT eigenvalues/eigenvectors: L* D = Lambda L*
	
// Printing out the eigenvalues - for control 
    if(PV)
    {
#ifdef MU
     fprintf(stdout, "\n\n %s eigenvalues: \n\n", (LR ? "LEFT" : "LARGEST"));
#else
     fprintf(stdout, "\n\n Largest eigenvalues: \n\n");
#endif
     for(i=0; i<evd.nev; i++)
      fprintf(stdout, "%i: %8.8E %8.8E\n", i, (double)CCALL(creal)(evd.evals[i]), (double)CCALL(cimag)(evd.evals[i]));
    };

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
