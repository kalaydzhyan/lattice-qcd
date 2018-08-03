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
	fprintf(stdout, "Usage:\t-i evfilename -c csp\n");
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
	int i;
	int res;
	int vol;
	t_real csp = 1.0;
	
	while ((c = getopt(argc, argv, "i:c:")) != -1) /*{{{*/
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
	vol *= evd.lt;
        vol *= evd.ndirac;
	vol *= evd.ncolors;
	
//  initializing the commutator of gammas
        t_complex d_gamma[DIM][NDIRAC][NDIRAC] = 
        {{
          { 0.0, 0.0,   0.0,-1.0*I},
          { 0.0, 0.0,-1.0*I,   0.0},
          { 0.0,   I,   0.0,   0.0},
          {   I, 0.0,   0.0,   0.0}
         },
         {
          { 0.0, 0.0, 0.0,-1.0},
          { 0.0, 0.0, 1.0, 0.0},
          { 0.0, 1.0, 0.0, 0.0},
          {-1.0, 0.0, 0.0, 0.0}
         },
         {
          { 0.0,   0.0,-1.0*I, 0.0},
          { 0.0,   0.0,   0.0,   I},
          {   I,   0.0,   0.0, 0.0},
          { 0.0,-1.0*I,   0.0, 0.0}
         },
         {
          { 0.0, 0.0, 1.0, 0.0},
          { 0.0, 0.0, 0.0, 1.0},
          { 1.0, 0.0, 0.0, 0.0},
          { 0.0, 1.0, 0.0, 0.0}
         }
        };

        t_complex d_gamma5[NDIRAC][NDIRAC] = 
        {{1.0, 0.0, 0.0, 0.0},
         {0.0, 1.0, 0.0, 0.0},
         {0.0, 0.0,-1.0, 0.0},
         {0.0, 0.0, 0.0,-1.0}
        };

        t_complex d_sigma[DIM][DIM][NDIRAC][NDIRAC];
        int mu, nu, id1, id2, id0;
        
        for(mu=0; mu<DIM; mu++)
         for(nu=0; nu<DIM; nu++)
          for(id1=0; id1<NDIRAC; id1++)
           for(id2=0; id2<NDIRAC; id2++)
           {
            d_sigma[mu][nu][id1][id2] = 0.0 + I*0.0;
            for(id0=0; id0<NDIRAC; id0++)
            {
             d_sigma[mu][nu][id1][id2] += I*d_gamma[mu][id1][id0]*d_gamma[nu][id0][id2];
             d_sigma[mu][nu][id1][id2] -= I*d_gamma[nu][id1][id0]*d_gamma[mu][id0][id2];
            }; 
           };
                
        t_complex d_g5gm[DIM][NDIRAC][NDIRAC];
        for(mu=0; mu<DIM; mu++)
         for(id1=0; id1<NDIRAC; id1++)
          for(id2=0; id2<NDIRAC; id2++)
          {
           d_g5gm[mu][id1][id2] = 0.0 + I*0.0;
           for(id0=0; id0<NDIRAC; id0++)
            d_g5gm[mu][id1][id2] += I*d_gamma5[id1][id0]*d_gamma[mu][id0][id2];
          }; 
                  
        t_complex d_sigma5[DIM][DIM][NDIRAC][NDIRAC];
                
        for(mu=0; mu<DIM; mu++)
         for(nu=0; nu<DIM; nu++)
          for(id1=0; id1<NDIRAC; id1++)
           for(id2=0; id2<NDIRAC; id2++)
           {
            d_sigma5[mu][nu][id1][id2] = 0.0 + I*0.0;
            for(id0=0; id0<NDIRAC; id0++)
            {
             d_sigma5[mu][nu][id1][id2] += d_gamma[mu][id1][id0]*d_g5gm[nu][id0][id2];
             d_sigma5[mu][nu][id1][id2] += d_g5gm[nu][id1][id0]*d_gamma[mu][id0][id2];
            }; 
           };    
        
        t_real tj[DIM], tj5[DIM];
        for(mu=0; mu<DIM; mu++)
        {
         tj[mu]  = 0.0;
         tj5[mu] = 0.0;
        }; 

	if(!((res & E_READ_CANTREAD) || (res & E_READ_CANTOPEN) || (res & E_READ_SECNOTHERE)))
        {
         for(i = 0; i < evd.nev; i++)
         {
          // Eigenvalue
          t_real lambda = fabs(2.0*CCALL(cimag)(evd.evals[i])/(2.0 - CCALL(creal)(evd.evals[i])));
          t_real signed_lambda = 2.0*440.0*CCALL(cimag)(evd.evals[i])/(csp*(2.0 - CCALL(creal)(evd.evals[i])));
          t_real dnorm  = fabs(CCALL(cabs)(evd.evals[i] - 1.0) - 1.0);
      
          fprintf(stdout, "l: %4.4lf, d: %4.6lf, ", (double)lambda, (double)dnorm);
          fflush(stdout);
          
          fprintf(stdout, "%2.6lf ", (double)signed_lambda);
      
          // Chirality
          t_real chirality = 0.0;
      
          int x, ic;
  
          for(x=0; x<VOL; x++)
           for(ic=0; ic<NCOLORS; ic++)
            for(id1=0; id1<NDIRAC; id1++)
             for(id2=0; id2<NDIRAC; id2++)
              chirality += CCALL(creal)(CCALL(conj)(*(evd.evecs + i*vol + NCOLORS*NDIRAC*x + NDIRAC*ic + id1))*d_gamma5[id1][id2]*(*(evd.evecs + i*vol + NCOLORS*NDIRAC*x + NDIRAC*ic + id2)));
              
          fprintf(stdout, "c: %4.4lf\n", (double)chirality);
          
          //fflush(stdout);
      
          // Magnetization
          t_real     asmn[DIM][DIM];
          t_complex  smn;
          t_real     asmn5[DIM][DIM];
          t_complex  smn5;
      
          for(mu=0; mu<DIM; mu++)
           for(nu=0; nu<DIM; nu++)
           {
            asmn[mu][nu]  = 0;
            asmn5[mu][nu] = 0;
           };    
  
          for(x=0; x<VOL; x++)
           for(mu=0; mu<DIM; mu++)
            for(nu=mu+1; nu<DIM; nu++)
            {
             smn  = 0.0 + I*0.0;
             smn5 = 0.0 + I*0.0;
             for(ic=0; ic<NCOLORS; ic++)
              for(id1=0; id1<NDIRAC; id1++)
               for(id2=0; id2<NDIRAC; id2++)
               {
                smn  += CCALL(conj)(*(evd.evecs + i*vol + NCOLORS*NDIRAC*x + NDIRAC*ic + id1))*d_sigma[mu][nu][id1][id2]*(*(evd.evecs + i*vol + NCOLORS*NDIRAC*x + NDIRAC*ic + id2));
                smn5 += CCALL(conj)(*(evd.evecs + i*vol + NCOLORS*NDIRAC*x + NDIRAC*ic + id1))*d_sigma5[mu][nu][id1][id2]*(*(evd.evecs + i*vol + NCOLORS*NDIRAC*x + NDIRAC*ic + id2));
               };
             asmn[mu][nu] += CCALL(creal)(smn);  
             asmn5[mu][nu] += CCALL(creal)(smn5);  
            };
      
            fprintf(stdout, "%2.6lf\n", (double)asmn[0][1]);
          
          
            
          // axial and vector currents
          t_real     aj[DIM];
          t_complex     j;
          t_real     aj5[DIM];
          t_complex     j5;
      
          for(mu=0; mu<DIM; mu++)
          {
           aj[mu]  = 0;
           aj5[mu] = 0;
          }; 
  
          for(x=0; x<VOL; x++)
           for(mu=0; mu<DIM; mu++)
           {    
            j  = 0.0 + I*0.0;
            j5 = 0.0 + I*0.0;
            for(ic=0; ic<NCOLORS; ic++)
             for(id1=0; id1<NDIRAC; id1++)
              for(id2=0; id2<NDIRAC; id2++)
               {
                j  += CCALL(conj)(*(evd.evecs + i*vol + NCOLORS*NDIRAC*x + NDIRAC*ic + id1))*d_gamma[mu][id1][id2]*(*(evd.evecs + i*vol + NCOLORS*NDIRAC*x + NDIRAC*ic + id2));
                j5 += CCALL(conj)(*(evd.evecs + i*vol + NCOLORS*NDIRAC*x + NDIRAC*ic + id1))*d_g5gm[mu][id1][id2]*(*(evd.evecs + i*vol + NCOLORS*NDIRAC*x + NDIRAC*ic + id2));
               };
            aj[mu]  += CCALL(creal)(j);
            aj5[mu] += CCALL(creal)(j5);  
           };
           
          //fprintf(stdout, "%2.6lf %2.6lf ", (double)aj[2], (double)aj5[2]); 
          
          if(lambda>1)
          {
           for(mu=0; mu<DIM; mu++)
           {
            tj[mu]  +=  aj[mu]/signed_lambda;
            tj5[mu] += aj5[mu]/signed_lambda;
           };
          };  
            
          //fprintf(stdout, "\n");        
          //for(mu=0; mu<DIM; mu++)
          // fprintf(stdout, "%i: %2.6lf %2.6lf\n", mu, (double)aj[mu], (double)aj5[mu]);
          //fprintf(stdout, "\n");
                    
          //fprintf(stdout, "\n");
          //for(mu=0; mu<DIM; mu++)
          // for(nu=mu+1; nu<DIM; nu++)
          // {
          //fprintf(stdout, "%i %i: %4.6lf\n", mu, nu, (double)(asmn[mu][nu]));
          // }; 
          
          //fprintf(stdout, "%2.6lf ", (double)chirality);
          //fprintf(stdout, "%2.6lf %2.6lf ", (double)asmn5[0][1], (double)asmn5[2][3]);
          
          //fprintf(stdout, "\n"); 
          //fflush(stdout);
         };
         
        };
        
        //fprintf(stdout, "\n\n\t VECTOR CURRENT: \t\n\n");
        //fprintf(stdout, "\n\n");
        //for(mu=0; mu<DIM; mu++)
        //{
        // fprintf(stdout, "%6.6lf  ", tj[mu]);
        //};
        
        //fprintf(stdout, "\n\n\t AXIAL CURRENT: \t\n\n");
        //for(mu=0; mu<DIM; mu++)
        //{
        // fprintf(stdout, "%6.6lf  ", tj5[mu]);
        //};
        //fprintf(stdout, "\n\n");
        //fflush(stdout);

	if (evd.evals != NULL) {
		free(evd.evals);
		evd.evals = NULL;
	}
	if (evd.evecs != NULL) {
		free(evd.evecs);
		evd.evecs = NULL;
	};	

	return EXIT_SUCCESS;
}/*}}}*/
