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

#undef NOISY_OUTPUT

void printhelp(void) __attribute__ ((noreturn));

void printhelp(void)
{
	fprintf(stdout, "Advanced Reader of EV files\n");
	
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
	
//  initializing the gamma matrices and various commutators
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
        int mu, nu, id1, id2, id0, ic;
        
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
	
	
	// Defining and initializing the observables
	
	int x;  
	
	double* ar5   = (double* )malloc(VOL*sizeof(double));   
	double* aj    = (double* )malloc(VOL*DIM*sizeof(double));
 	double* asmn  = (double* )malloc(VOL*DIM*DIM*sizeof(double));
	
	for(x = 0; x<VOL; x++)
	{
	 ar5[x]  = 0.0; 
	 for(mu = 0; mu<DIM; mu++)
	  aj[DIM*x + mu] = 0.0;
	 for(mu = 0; mu<DIM; mu++)
	  for(nu = 0; nu<DIM; nu++)
	   asmn[x*DIM*DIM + mu*DIM + nu] = 0.0;
	};
	
	t_real* lambda = (t_real *)malloc(evd.nev*sizeof(t_real));

	if(!((res & E_READ_CANTREAD) || (res & E_READ_CANTOPEN) || (res & E_READ_SECNOTHERE)))
         for(i = 0; i < evd.nev; i++)
         {
          lambda[i] = 2.0*CCALL(cimag)(evd.evals[i])/(2.0 - CCALL(creal)(evd.evals[i]));
          
#ifdef NOISY_OUTPUT
          fprintf(stdout, "lambda = %6.6lf + I*%6.6lf\n", (double)CCALL(creal)(lambda[i]), (double)CCALL(cimag)(lambda[i]));
          fflush(stdout);
#endif
         };
        
        // Chirality 
        // Define zero mode + smallest nonzero mode
        double  mlambda = 10000.0;
        int    nmlambda = 0;
        for(i=0; i<evd.nev; i++)
         if((lambda[i] < mlambda) && (lambda[i] > 0.1))
         {
          mlambda  = lambda[i];
          nmlambda = i;
         }; 
#ifdef NOISY_OUTPUT
	fprintf(stdout, "Minimal nonzero lambda: %6.6lf\n", mlambda);
	fflush(stdout);	
#endif
	int id;
        for(i=0; i<evd.nev; i++)
         if(lambda[i] < 0.1 || i==nmlambda)
         {
          t_complex r5 = 0.0 + I*0.0;
          for(x = 0; x<VOL; x++)
          {
           for(ic=0; ic<NCOLORS; ic++)
            for(id=0; id<NDIRAC; id++)
             r5 += CCALL(conj)(*(evd.evecs + i*vol + NCOLORS*NDIRAC*x + NDIRAC*ic + id))*(*(evd.evecs + i*vol + NCOLORS*NDIRAC*x + NDIRAC*ic + id))*(id < 2? 1.0 : -1.0);
           ar5[x] += (double)CCALL(creal)(r5);    
          };     
         };
         
        // End of chirality
        
        // magnetization

        for(i=0; i<evd.nev; i++)
         if(lambda[i] < 0.1 || i==nmlambda)
         {
          t_complex smn;    
	  for(x = 0; x<VOL; x++)
           for(mu=0; mu<DIM; mu++)
            for(nu=mu+1; nu<DIM; nu++)
            {
             smn = 0.0 + I*0.0;
             for(ic=0; ic<NCOLORS; ic++)
              for(id1=0; id1<NDIRAC; id1++)
               for(id2=0; id2<NDIRAC; id2++)
                smn += CCALL(conj)(*(evd.evecs + i*vol + NCOLORS*NDIRAC*x + NDIRAC*ic + id1))*d_sigma[mu][nu][id1][id2]*(*(evd.evecs + i*vol + NCOLORS*NDIRAC*x + NDIRAC*ic + id2));
              asmn[x*DIM*DIM + DIM*mu + nu] += (double)CCALL(creal)(smn);
            };
         };         
         
        if(!((res & E_READ_CANTREAD) || (res & E_READ_CANTOPEN) || (res & E_READ_SECNOTHERE)))
         for(i = 0; i < evd.nev; i++) 
         {
	  // Checking whether a given eigenvalue is nonzero and has pair
	  int j, has_pair = 0;
	  for(j=0; j<evd.nev; j++)
	   if(j!=i && fabs(lambda[j]+lambda[i])<1 && fabs(lambda[i])>0.1 )
	    has_pair = 1;
#ifdef NOISY_OUTPUT
	  if(!has_pair)
	   fprintf(stdout, "lambda = %4.4lf has no pair!\n", lambda[i]);
#endif
	  if(has_pair)
	  {   
           // axial and vector currents
           
           t_complex     j, j5;
	   for(x = 0; x<VOL; x++)
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
              aj[x*DIM + mu]  += (double)CCALL(creal)(j)/lambda[i];
            };    

          };
         };  

        // Processing the data
        
        double ar52   = 0.0; //rho_5^2
        double ar5s23 = 0.0; //rho5 sigma23
        double ajm[DIM]; // j_mu
        double ajm2[DIM]; // j_mu^2
        double ar52jm2[DIM]; // rho_5^2 j_mu^2
        double dmar52[DIM]; // (d_mu rho_5)^2
        double asmn1[DIM][DIM], asmn2[DIM][DIM]; //s_munu
        double dmar52jn2[DIM][DIM]; //(d_mu rh05)^2 j_nu^2
        
        FILE* cfile = fopen("instanton_current.dat","w");
        for(x=0; x<VOL; x++)
        {
         for(mu=0; mu<DIM; mu++)
          fprintf(cfile, "%4.6E  ", aj[x*DIM + mu]);
         fprintf(cfile, "\n"); 
        };  
        fclose(cfile);  
        
        FILE* mfile = fopen("instanton_magnetization.dat","w");
        for(x=0; x<VOL; x++)
        {
         for(mu=0; mu<DIM-1; mu++)
          for(nu=mu+1; nu<DIM; nu++)
           fprintf(mfile, "%4.6E  ", asmn[x*DIM*DIM + DIM*mu + nu]);
         fprintf(mfile, "\n"); 
        };  
        fclose(mfile);  
        
        FILE* pfile = fopen("instanton_chirality.dat","w");
        for(x=0; x<VOL; x++)
        {
         fprintf(pfile, "%4.6E  ", ar5[x]);
         fprintf(pfile, "\n"); 
        };  
        fclose(pfile);  
                
        for(mu=0; mu<DIM; mu++)
         ajm[mu] = 0.0; 
        for(mu=0; mu<DIM; mu++)
         ajm2[mu] = 0.0; 
        for(mu=0; mu<DIM; mu++)
         ar52jm2[mu] = 0.0; 
        for(mu=0; mu<DIM; mu++)
         dmar52[mu] = 0.0; 
        for(mu=0; mu<DIM; mu++)
         for(nu=0; nu<DIM; nu++)
          asmn1[mu][nu] = 0.0; 
        for(mu=0; mu<DIM; mu++)
         for(nu=0; nu<DIM; nu++)
          asmn2[mu][nu] = 0.0; 
        for(mu=0; mu<DIM; mu++)
         for(nu=0; nu<DIM; nu++)
          dmar52jn2[mu][nu] = 0.0; 
          
        lat_mov_init();  
 
        for(x=0; x<VOL; x++)
        {
         ar52 += ar5[x]*ar5[x];
         
         ar5s23 += ar5[x]*asmn[x*DIM*DIM + 2*DIM + 3];
         
         for(mu=0; mu<DIM; mu++)
          ajm[mu] += aj[x*DIM + mu];
         
         for(mu=0; mu<DIM; mu++)
          ajm2[mu] += aj[x*DIM + mu]*aj[x*DIM + mu];
      
         for(mu=0; mu<DIM; mu++)
          ar52jm2[mu] += ar5[x]*ar5[x]*aj[x*DIM + mu]*aj[x*DIM + mu];
          
         for(mu=0; mu<DIM; mu++)
         {
          int mx1 = (*lat_mov)[x][mu][FWD];
          int mx2 = (*lat_mov)[x][mu][BWD];
          dmar52[mu] += (ar5[mx1] - ar5[mx2])*(ar5[mx1] - ar5[mx2]);
         }; 
          
         for(mu=0; mu<DIM; mu++)
          for(nu=mu+1; nu<DIM; nu++)
           asmn1[mu][nu] += asmn[x*DIM*DIM + DIM*mu + nu];  
           
         for(mu=0; mu<DIM; mu++)
          for(nu=mu+1; nu<DIM; nu++)
           asmn2[mu][nu] += asmn[x*DIM*DIM + DIM*mu + nu]*asmn[x*DIM*DIM + DIM*mu + nu];  
           
         for(mu=0; mu<DIM; mu++)
          for(nu=0; nu<DIM; nu++)
          {
           int mx1 = (*lat_mov)[x][mu][FWD];
           int mx2 = (*lat_mov)[x][mu][BWD];
           dmar52jn2[mu][nu] += (ar5[mx1] - ar5[mx2])*(ar5[mx1] - ar5[mx2])*aj[x*DIM + nu]*aj[x*DIM + nu];  
           
          }; 
           
        };
        
#ifdef NOISY_OUTPUT        
	fprintf(stdout,"ar5^2: \n");
        fprintf(stdout, "%6.4lE", ar52/(double)VOL);
	fprintf(stdout,"\n\n");
#else
	fprintf(stdout, "%6.4lE ", ar52/(double)VOL);
#endif         

#ifdef NOISY_OUTPUT        
	fprintf(stdout,"ar5s23: \n");
        fprintf(stdout, "%6.4lE", ar5s23/(double)VOL);
	fprintf(stdout,"\n\n");
#else
	fprintf(stdout, "%6.4lE ", ar5s23/(double)VOL);
#endif         

#ifdef NOISY_OUTPUT
        fprintf(stdout, "ajm: \n"); 
        for(mu=0; mu<DIM; mu++)
         fprintf(stdout, "%i: %6.4lE ", mu, ajm[mu]/(double)VOL);
        fprintf(stdout,"\n\n"); 
#else
	for(mu=0; mu<DIM; mu++)
         fprintf(stdout, "%6.4lE ", ajm[mu]/(double)VOL);
#endif

#ifdef NOISY_OUTPUT
        fprintf(stdout, "ajm^2: \n"); 
        for(mu=0; mu<DIM; mu++)
         fprintf(stdout, "%i: %6.4lE ", mu, ajm2[mu]/(double)VOL);
        fprintf(stdout,"\n\n"); 
#else
	for(mu=0; mu<DIM; mu++)
         fprintf(stdout, "%6.4lE ", ajm2[mu]/(double)VOL);
#endif

#ifdef NOISY_OUTPUT
	fprintf(stdout, "ar52jm2: \n");  
        for(mu=0; mu<DIM; mu++)
         fprintf(stdout, "%i: %6.4lE ", mu, ar52jm2[mu]/(double)VOL);
	fprintf(stdout,"\n\n");
#else
	for(mu=0; mu<DIM; mu++)
         fprintf(stdout, "%6.4lE ", ar52jm2[mu]/(double)VOL);
#endif

#ifdef NOISY_OUTPUT
	fprintf(stdout, "(d_mu ar5)^2: \n");  
        for(mu=0; mu<DIM; mu++)
         fprintf(stdout, "%i: %6.4lE ", mu, dmar52[mu]/(double)VOL);
	fprintf(stdout,"\n\n");
#else
	for(mu=0; mu<DIM; mu++)
         fprintf(stdout, "%6.4lE ", dmar52[mu]/(double)VOL);
#endif

#ifdef NOISY_OUTPUT
	fprintf(stdout, "asmn: \n");          
        for(mu=0; mu<DIM; mu++)
         for(nu=mu+1; nu<DIM; nu++)
          fprintf(stdout, "%i%i: %6.4lE ", mu, nu, asmn1[mu][nu]/(double)VOL);  
	fprintf(stdout,"\n\n");
#else
	for(mu=0; mu<DIM; mu++)
         for(nu=mu+1; nu<DIM; nu++)
          fprintf(stdout, "%6.4lE ", asmn1[mu][nu]/(double)VOL);
#endif        

#ifdef NOISY_OUTPUT                  
	fprintf(stdout, "asmn^2: \n");
        for(mu=0; mu<DIM; mu++)
         for(nu=mu+1; nu<DIM; nu++)
          fprintf(stdout, "%i%i: %6.4lE ", mu, nu, asmn2[mu][nu]/(double)VOL);  
	fprintf(stdout,"\n\n");
#else
	for(mu=0; mu<DIM; mu++)
         for(nu=mu+1; nu<DIM; nu++)
          fprintf(stdout, "%6.4lE ", asmn2[mu][nu]/(double)VOL);
#endif

#ifdef NOISY_OUTPUT                  
	fprintf(stdout, "(d_mu ar5)^2(j_nu)^2: \n");
        for(mu=0; mu<DIM; mu++)
         for(nu=0; nu<DIM; nu++)
          fprintf(stdout, "%i%i: %6.4lE ", mu, nu, dmar52jn2[mu][nu]/(double)VOL);  
	fprintf(stdout,"\n\n");
#else
	for(mu=0; mu<DIM; mu++)
         for(nu=0; nu<DIM; nu++)
          fprintf(stdout, "%6.4lE ", dmar52jn2[mu][nu]/(double)VOL);
#endif

        fprintf(stdout, "\n");
        fflush(stdout);  
        
	free(ar5);
	free(asmn); 
	free(aj);
	free(lambda);

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
