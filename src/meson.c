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
#include "statistics.h"

#ifdef SU2
#include "su2math.h"
#endif

#ifdef SU3
#include "su3math.h"
#endif

#undef  NOISY_OUTPUT
#define MAXINAMES    (150)
#define MAX_ZERO_EV  (0.0001) // Boundary between "exactly zero" and nonzero eigenvalues

#define EV_SIZE      (VOL*NDIRAC*NCOLORS)

void set_czero(t_complex *A, int n)
{
 int i;
 for(i = 0; i<n; i++)
  A[i] = 0.0 + I*0.0;
}

void set_dzero(t_real *A, int n)
{
 int i;
 for(i = 0; i<n; i++)
  A[i] = 0.0;
}

char      outfname_suffix[300];

// Initializing the gamma matrices and various commutators 
// Tensor -> 6 components
int munu(int mu, int nu)
{
 if(mu==0 && nu==1)
  return 0;
 if(mu==0 && nu==2)
  return 1;
 if(mu==0 && nu==3)
  return 2;
 if(mu==1 && nu==2)
  return 3;
 if(mu==1 && nu==3)
  return 4;
 if(mu==2 && nu==3)
  return 5;
 return -1;
}

 t_complex dgamma[DIM][NDIRAC][NDIRAC] = {
  {
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

t_complex dunity[NDIRAC][NDIRAC] = {
   {1.0, 0.0, 0.0, 0.0},
   {0.0, 1.0, 0.0, 0.0},
   {0.0, 0.0, 1.0, 0.0},
   {0.0, 0.0, 0.0, 1.0},
  };
  
t_complex gamma5[NDIRAC][NDIRAC] = {
   {1.0, 0.0, 0.0, 0.0},
   {0.0, 1.0, 0.0, 0.0},
   {0.0, 0.0,-1.0, 0.0},
   {0.0, 0.0, 0.0,-1.0},
  };
  
t_complex sigma[6][NDIRAC][NDIRAC];

void initialize_gamma()
{
 int mu, nu, id1, id2, id0;
 
 for(mu=0; mu<DIM-1; mu++)
  for(nu=mu+1; nu<DIM; nu++)
   for(id1=0; id1<NDIRAC; id1++)
    for(id2=0; id2<NDIRAC; id2++)
    {
     sigma[munu(mu, nu)][id1][id2] = 0.0 + I*0.0;
     for(id0=0; id0<NDIRAC; id0++)
     {
      sigma[munu(mu, nu)][id1][id2] += I*dgamma[mu][id1][id0]*dgamma[nu][id0][id2];
      sigma[munu(mu, nu)][id1][id2] -= I*dgamma[nu][id1][id0]*dgamma[mu][id0][id2];
     }; 
    };
};

// Initializing the OVD structure
void init_evd(t_op_ev *evd)
{
 evd->op = NULL;
 evd->nev = 0;
 evd->which[0] = 0;
 evd->which[1] = 0;
 evd->val_vec = 0;
 evd->tol = 0;
 evd->evals = NULL;
 evd->evecs = NULL;
 evd->data = 0;
}

void free_evd(t_op_ev *evd)
{
 if (evd->evals != NULL)
 {
  free(evd->evals);
  evd->evals = NULL;
 };
 if (evd->evecs != NULL)
 {
  free(evd->evecs);
  evd->evecs = NULL;
 };
}

int main (int argc, char * argv[]) // format of input: meson prefix suffix datafiles
{ 
 int ifile, t, t1, t2, mu, xs, x, id1, id2, id3, id4, a, b, isrc1, isrc2;
 t_complex tmp;
 t_complex  jxjy[DIM][LT];
 t_real    ajxjy[DIM][LT];
 t_real   a2jxjy[DIM][LT];
 t_real    cjxjy[DIM][LT][LT];
 t_real    cc = 0.0, dcc = 0.0;
 
 
 //Split inames
 int nfiles = argc - 3;
 if(nfiles<=0)
  panic("%s: no input files specified!!!\n", __func__);
  
 initialize_gamma(); 
 
 set_dzero(    &(ajxjy[0][0]), DIM*LT);
 set_dzero(   &(a2jxjy[0][0]), DIM*LT);
 set_dzero( &(cjxjy[0][0][0]), DIM*LT*LT);
 
 // Loop over all files
 for(ifile = 0; ifile < nfiles; ifile++)
 {
  t_op_ev evd;
  init_evd(&evd);
  fprintf(stdout, "Processing the file %s...\n", argv[ifile + 3]);
  int res = read_ev(argv[ifile + 3], &evd, SM_SECTION);
  if((res & E_READ_CANTREAD) || (res & E_READ_CANTOPEN) || (res & E_READ_SECNOTHERE))
   panic("%s: cannot read the file %s!!!\n", __func__, argv[ifile + 3]);
  if(evd.nev!=NDIRAC*NCOLORS)
   panic("%s: nev != no. of sources in file %s!!!\n", __func__, argv[ifile + 3]);
  
  set_czero(&(jxjy[0][0]), DIM*LT);
  
  /* Estimating the chiral condensate */
  
  x = 0;
  tmp = 0.0 + I*0.0;
  for(id1=0; id1<NDIRAC; id1++)
   for(a=0; a<NCOLORS; a++)
   {
    isrc1  = NDIRAC*a + id1;
    tmp   += evd.evecs[isrc1*EV_SIZE + x*NCOLORS*NDIRAC + NDIRAC*a + id1];        
   };
  
  cc  += CCALL(creal)(tmp);
  dcc += CCALL(creal)(tmp)*CCALL(creal)(tmp);  
  
  /* Now go the currents */

  for(t=0; t<LT; t++)
   for(mu=0; mu<DIM; mu++) 
    for(xs=0; xs<LS3; xs++)
     for(id1=0; id1<NDIRAC; id1++)
      for(id2=0; id2<NDIRAC; id2++)
       for(id3=0; id3<NDIRAC; id3++)
        for(id4=0; id4<NDIRAC; id4++)
         for(a=0; a<NCOLORS; a++)
          for(b=0; b<NCOLORS; b++)
          {
           x = t*LS3 + xs;        
           isrc1 = NDIRAC*b + id2;
           isrc2 = NDIRAC*b + id3;
           tmp   =             evd.evecs[isrc1*EV_SIZE + x*NCOLORS*NDIRAC + NDIRAC*a + id1];
           tmp  *= CCALL(conj)(evd.evecs[isrc2*EV_SIZE + x*NCOLORS*NDIRAC + NDIRAC*a + id4]);
           tmp  *= dgamma[mu][id2][id3]*dgamma[mu][id4][id1];
           tmp  *= (id3 < 2)? 1.0 : -1.0;
           tmp  *= (id4 < 2)? 1.0 : -1.0; 
           jxjy[mu][t] += tmp;         
          };
    
  for(t=0; t<LT; t++)
   for(mu=0; mu<DIM; mu++)
    ajxjy[mu][t] += CCALL(creal)(jxjy[mu][t]);
    
  for(t=0; t<LT; t++)
   for(mu=0; mu<DIM; mu++)
    a2jxjy[mu][t] += CCALL(creal)(jxjy[mu][t])*CCALL(creal)(jxjy[mu][t]);  
  
  for(t1=0; t1<LT; t1++)
   for(t2=0; t2<LT; t2++)
    for(mu=0; mu<DIM; mu++)
     cjxjy[mu][t1][t2] += CCALL(creal)(jxjy[mu][t1]*jxjy[mu][t2]);
          
  //freeing up the memory
  free_evd(&evd);
 };
 
 char corrfname[300], covfname[300]; // correlator + covariance matrix 
 sprintf(corrfname, "%scorrelator_%s.dat", argv[1], argv[2]);
 sprintf(covfname,  "%scovariance_%s.dat", argv[1], argv[2]);
 
 fprintf(stdout, "Correlator file: %s\n", corrfname);
 fprintf(stdout, "Covariance file: %s\n",  covfname);
 
 // Normalizing and processing
 for(t=0; t<LT; t++)
  for(mu=0; mu<DIM; mu++)
   ajxjy[mu][t] = ajxjy[mu][t]/(t_real)nfiles;
   
 for(t=0; t<LT; t++)
  for(mu=0; mu<DIM; mu++)
   a2jxjy[mu][t] = sqrt(a2jxjy[mu][t]/(t_real)nfiles - ajxjy[mu][t]*ajxjy[mu][t])/sqrt((t_real)nfiles);  
  
 for(t1=0; t1<LT; t1++)
  for(t2=0; t2<LT; t2++)
   for(mu=0; mu<DIM; mu++)
    cjxjy[mu][t1][t2] = cjxjy[mu][t1][t2]/(t_real)nfiles - ajxjy[mu][t1]*ajxjy[mu][t2];
    
 // Output
 FILE* corrfile = fopen(corrfname, "w");
 for(t=0; t<LT; t++)
 {
  fprintf(corrfile, "%i ", t); 
  for(mu=0; mu<DIM; mu++)
   fprintf(corrfile, "%6.6E  %6.6E  ", ajxjy[mu][t], a2jxjy[mu][t]);
  fprintf(corrfile, "\n"); 
 };  
 fclose(corrfile); 
 
 FILE* covfile = fopen(covfname, "w");
 for(mu=0; mu<DIM; mu++)
 {
  for(t1=0; t1<LT; t1++)
  {
   for(t2=0; t2<LT; t2++) 
    fprintf(covfile, "%6.6E ", cjxjy[mu][t1][t2]);
   fprintf(covfile, "\n"); 
  }; 
  fprintf(covfile, "\n\n");  
 };  
 fclose(covfile);
      
 fprintf(stdout, "\n\n");
  for(mu=0; mu<DIM; mu++)
   for(t=0; t<LT; t++) 
    fprintf(stdout, "\t < jxjy[mu = %i][t = %i] >: %6.6E +/- %6.6E\n", mu, t, ajxjy[mu][t], a2jxjy[mu][t]);        
  fprintf(stdout, "\n\n");
   
 cc = cc/(t_real)nfiles;
 dcc = sqrt(dcc/(t_real)nfiles - cc*cc)/sqrt((t_real)nfiles);
 
 t_real ccs = (0.2/0.1)*(0.2/0.1)*(0.2/0.1)/1.4;
 
 t_real cc31 = pow(ccs*cc + ccs*dcc, 0.333);
 t_real cc32 = pow(ccs*cc - ccs*dcc, 0.333);
 
 cc  = 0.5*(cc31 + cc32);
 dcc = 0.5*fabs(cc31 - cc32);
  
 fprintf(stdout, "Chiral condensate (assuming a = 0.1 fm): %4.4E +/- %4.4E GeV\n", cc, dcc);
 
 return EXIT_SUCCESS;
}
