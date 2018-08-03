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

#define NMASSES      (3)

#define NCORRS	     (10)

int x_shift(int position, int delta, int direction)
{
 int pos = position, x,y,z,t;

 t = pos / (LS*LS*LS);
 pos %= (LS*LS*LS);
 x = pos / (LS*LS);
 pos %= (LS*LS);
 y = pos / LS;
 pos %= LS;
 z = pos;

 switch (direction)
        {

         case 0:
         t = (t + delta) % LT;
         break;

         case 1:
         x = (x + delta) % LS;
         break;

         case 2:
         y = (y + delta) % LS;
         break;

         case 3:
         z = (z + delta) % LS;
         break;
        }

 return (t*LS*LS*LS + x*LS*LS + y*LS + z);

}



void set_czero(t_complex *A, int n)
{
 int i;
 for(i = 0; i<n; i++)
  A[i] = 0.0 + 0.0*I;    
}

void set_dzero(double *A, int n)
{
 int i;
 for(i = 0; i<n; i++)
  A[i] = 0.0;
}

// Simulation parameters
double    beta = 0.0;
double    spacing = 0.0;
int       H = 0;
double    Hphys = 0.0;
char      outfname_suffix[300];
int       max_nz_evals = 0;
double    mass_phys = 0.0;
t_complex mass_latt = 0.0 + I*0.0;

void read_sim_params(char* param_file)
{
 FILE *params = fopen(param_file, "r");
 fscanf(params, "%lf", &beta);
 fscanf(params, "%lf", &spacing);
 fscanf(params, "%i",  &H);
 fscanf(params, "%lf", &mass_phys); // Physical mass in MeV - regulator for banks-casher etc.
 fscanf(params, "%i",  &max_nz_evals);
 fscanf(params, "%s", outfname_suffix);
 fclose(params);

 Hphys     = 3*2*3.14159*H*(0.197326*0.197326)/(LS*LS*spacing*spacing); //3 = 1/q, q is the quark charge
 mass_latt = (t_real)(spacing*mass_phys/197.326) + 0.0*I;

 fprintf(stdout, "\n\t RUN PARAMETERS: \t\n\n");
 fprintf(stdout, "\t beta:\t\t %2.4lf\n", beta);
 fprintf(stdout, "\t spacing:\t %2.4lf fm\n", spacing);
 fprintf(stdout, "\t H (lat.u.):\t %i\n", H);
 fprintf(stdout, "\t H (ph.u.):\t %2.4lf GeV^2\n", Hphys);
 fprintf(stdout, "\t Mass (ph.u.):\t %3.2lf MeV\n", mass_phys);
 fprintf(stdout, "\t Mass (lat.u.):\t %3.2lf + I*%3.2lf MeV\n", (double)(CCALL(creal)(mass_latt)), (double)(CCALL(cimag)(mass_latt)) );
 fprintf(stdout, "\t Nonzero ev.:\t %i \n", max_nz_evals);
 fprintf(stdout, "\t Suffix:\t %s\n", outfname_suffix);
 fprintf(stdout, "\n\n");
}

//  Initializing the gamma matrices and various commutators 
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

// Projecting eigenvalues
void project_ev(t_complex* evals, int nev)
{
 int i;
 for(i=0; i<nev; i++)
  evals[i] = 2.0*evals[i]/(2.0 - evals[i]);
}

// Sorting the eigenvalues
int sort_ev(t_complex* revecs, t_complex* evals, int nev)
{	 
 t_cds_vector* tmp;
 if(posix_memalign((void*)&tmp, 16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for tmp (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
#ifdef NOISY_OUTPUT
 fprintf(stdout, "Sorting eigenvalues ... \n");
 fflush(stdout);
#endif
 int i, k, ki;
 for(i=0; i<nev; i++)
 {
  t_real minlambda = 100000000;
  ki = 0;
  for(k=i; k<nev; k++)
   if(CCALL(cabs)(evals[k]) < minlambda)
   {
    minlambda = CCALL(cabs)(evals[k]);
    ki = k;
   };
  k = ki;
  memcpy(                            tmp, &(revecs[i*VOL*NCOLORS*NDIRAC]), VOL*NCOLORS*NDIRAC*sizeof(t_complex));
  memcpy(&(revecs[i*VOL*NCOLORS*NDIRAC]), &(revecs[k*VOL*NCOLORS*NDIRAC]), VOL*NCOLORS*NDIRAC*sizeof(t_complex));
  memcpy(&(revecs[k*VOL*NCOLORS*NDIRAC]),                             tmp, VOL*NCOLORS*NDIRAC*sizeof(t_complex));
  t_complex tmpval = evals[i];
  evals[i] = evals[k];
  evals[k] = tmpval;
 };
 free(tmp);

 int nzero = 0;
 while(CCALL(cabs)(evals[nzero])<MAX_ZERO_EV)
  nzero++;

#ifdef NOISY_OUTPUT
 fprintf(stdout, "DONE: sorting left evals/evecs!\n");
 fflush(stdout);
 for(i=0; i<nev; i++)
  fprintf(stdout, "%i: \t%2.6lf %2.6lf\n", i, CCALL(creal)(evals[i]), CCALL(cimag)(evals[i]));
 fprintf(stdout, "Number of zero modes: %i\n", nzero);
 fflush(stdout);
#endif
 return nzero;
}

// op - a link to the array of dim n_op_idx X NDIRAC X NDIRAC

t_complex R[NDIRAC][NDIRAC][NDIRAC][NDIRAC][NCORRS];

void init_props(t_complex *revecs, t_complex *evals, int nev, t_complex mass)
{
 int i, x, a, b, id1, id2, id3, id4, t;

 t_complex DCS1[NDIRAC][NDIRAC], DCS2[NDIRAC][NDIRAC][NCORRS];
 t_complex CS1[NCOLORS][NCOLORS][NDIRAC][NDIRAC][NCORRS], CS2[NCOLORS][NCOLORS][NDIRAC][NDIRAC][NCORRS];
 
 set_czero(&(R[0][0][0][0][0]), NDIRAC*NDIRAC*NDIRAC*NDIRAC*NCORRS);

 for(x=0; x<VOL; x++)
 {
  set_czero(&(DCS1[0][0]),         NDIRAC*NDIRAC);
  set_czero(&(DCS2[0][0][0]),      NDIRAC*NDIRAC*NCORRS);
  set_czero(&(CS1[0][0][0][0][0]), NCOLORS*NCOLORS*NDIRAC*NDIRAC*NCORRS);
  set_czero(&(CS2[0][0][0][0][0]), NCOLORS*NCOLORS*NDIRAC*NDIRAC*NCORRS);
    
  /* Initializing Disconnected Prop 1 - G(x,x) */    
  for(i=0; i<nev; i++)
   for(a=0; a<NCOLORS; a++)
    for(id1=0; id1<NDIRAC; id1++)
     for(id2=0; id2<NDIRAC; id2++)
      DCS1[id1][id2] += CCALL(conj)(revecs[i*EV_SIZE + NCOLORS*NDIRAC*x + NDIRAC*a + id1])*revecs[i*EV_SIZE + NCOLORS*NDIRAC*x + NDIRAC*a + id2]/(evals[i] + mass);

  /* Initializing Disconnected Prop 2 - G(x + t,x + t) */    
  for(t=0; t<NCORRS; t++)    
   for(i=0; i<nev; i++)
    for(a=0; a<NCOLORS; a++)
     for(id1=0; id1<NDIRAC; id1++)
      for(id2=0; id2<NDIRAC; id2++)
       DCS2[id1][id2][t] += CCALL(conj)(revecs[i*EV_SIZE + NCOLORS*NDIRAC*x_shift(x,t,0) + NDIRAC*a + id1])*revecs[i*EV_SIZE + NCOLORS*NDIRAC*x_shift(x,t,0) + NDIRAC*a + id2]/(evals[i] + mass);    

  /* Initializing Connected Prop 1 - G(x, x+t) */
  for(t=0; t<NCORRS; t++)    
   for(i=0; i<nev; i++)
    for(a=0; a<NCOLORS; a++)
     for(b=0; b<NCOLORS; b++)
      for(id1=0; id1<NDIRAC; id1++)
       for(id2=0; id2<NDIRAC; id2++)
        CS1[a][b][id1][id2][t] += CCALL(conj)(revecs[i*EV_SIZE + NCOLORS*NDIRAC*x + NDIRAC*a + id1])*revecs[i*EV_SIZE + NCOLORS*NDIRAC*x_shift(x,t,0) + NDIRAC*b + id2]/(evals[i] + mass);

  /* Initializing Connected Prop 2 - G(x+t, x) */
  for(t=0; t<NCORRS; t++)
   for(i=0; i<nev; i++)
    for(a=0; a<NCOLORS; a++)
     for(b=0; b<NCOLORS; b++)
      for(id1=0; id1<NDIRAC; id1++)
       for(id2=0; id2<NDIRAC; id2++)
        CS2[a][b][id1][id2][t] += CCALL(conj)(revecs[i*EV_SIZE + NCOLORS*NDIRAC*x_shift(x,t,0) + NDIRAC*a + id1])*revecs[i*EV_SIZE + NCOLORS*NDIRAC*x + NDIRAC*b + id2]/(evals[i] + mass);

  /* Disconnected part of two-point current-current correlator */
  for(t=0; t<NCORRS; t++)
   for(id1=0; id1<NDIRAC; id1++)
    for(id2=0; id2<NDIRAC; id2++)
     for(id3=0; id3<NDIRAC; id3++)
      for(id4=0; id4<NDIRAC; id4++)        
       R[id1][id2][id3][id4][t] += DCS1[id1][id2]*DCS2[id3][id4][t];    
       
  /* Disconnected part of two-point current-current correlator */
  for(t=0; t<NCORRS; t++)
   for(id1=0; id1<NDIRAC; id1++)
    for(id2=0; id2<NDIRAC; id2++)
     for(id3=0; id3<NDIRAC; id3++)
      for(id4=0; id4<NDIRAC; id4++)
       for(a=0; a<NCOLORS; a++)
        for(b=0; b<NCOLORS; b++)
         R[id1][id2][id3][id4][t] -= CS1[a][b][id1][id4][t]*CS2[b][a][id3][id2][t];
  /* Progress indicator */        
  if(x%LS3==0)
  {         
   fprintf(stdout, ".");
   fflush(stdout);
  }; 
 };
 
 for(t=0; t<NCORRS; t++)  
  for(id1=0; id1<NDIRAC; id1++)
   for(id2=0; id2<NDIRAC; id2++)
    for(id3=0; id3<NDIRAC; id3++)
     for(id4=0; id4<NDIRAC; id4++)
      R[id1][id2][id3][id4][t] = R[id1][id2][id3][id4][t]/(t_real)VOL;
}

int main (int argc, char * argv[])
{ 
 int ifile;
 
 t_complex* revecs;
 t_complex* revals;
 int incl_zm = 0; // TODO: input by user
 
 //Split inames
 int nfiles = argc - 1;
 if(!nfiles)
  panic("%s: no input files specified!!!\n", __func__);
  
 initialize_gamma(); 
 read_sim_params("params.in");

 // Loop over all files
 for(ifile = 0; ifile < nfiles; ifile++)
 {
  t_op_ev evd;
  init_evd(&evd);
  
  fprintf(stdout, "Processing the file %s...\n", argv[ifile + 1]);
  int res = read_ev(argv[ifile + 1], &evd, SM_SECTION);
  if((res & E_READ_CANTREAD) || (res & E_READ_CANTOPEN) || (res & E_READ_SECNOTHERE))
   panic("%s: cannot read the file %s!!!\n", __func__, argv[ifile + 1]);
  revecs = (t_complex *)malloc(evd.nev*VOL*NCOLORS*NDIRAC*sizeof(t_complex));
  revals = (t_complex *)malloc(evd.nev*sizeof(t_complex)); 
  memcpy(revecs, &(evd.evecs[0]), evd.nev*VOL*NCOLORS*NDIRAC*sizeof(t_complex));
  memcpy(revals, &(evd.evals[0]), evd.nev*sizeof(t_complex)); 
  
  if(evd.evecs!=NULL)
  {
   free(evd.evecs);
   evd.evecs = NULL;
  };
  
  if(evd.evals!=NULL)
  {
   free(evd.evals);
   evd.evals = NULL;
  }; 
 
  project_ev(&(revals[0]), evd.nev);
  int nzero = sort_ev(&(revecs[0]), &(revals[0]), evd.nev); //TODO: sort_ev for complex spectrum???
  int first_mode = (incl_zm)? 0 : nzero;
  int n_use_evals = 2*max_nz_evals + ((first_mode)? 0 : nzero);
  if(n_use_evals>evd.nev)
   panic("%s: I do not have %i evs, only %i!!!\n", __func__, n_use_evals, evd.nev);

  init_props(&(revecs[first_mode*VOL*NCOLORS*NDIRAC]), &(revals[first_mode]), n_use_evals, mass_latt);
  
  //freeing up the memory
  free_evd(&evd);
  free(revecs);
  free(revals);
 };
 
 return EXIT_SUCCESS;
}
