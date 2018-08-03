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
#ifdef MU
double    MuPhys  = 0.0;
#endif
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
#ifdef MU
 double tmpmu;
 fscanf(params, "%lf", &tmpmu);
 Mu = (t_real)tmpmu;
#endif
 fscanf(params, "%lf", &mass_phys); // Physical mass in MeV - regulator for banks-casher etc.
 fscanf(params, "%i",  &max_nz_evals);
 fscanf(params, "%s", outfname_suffix);
 fclose(params);

 Hphys     = 3*2*3.14159*H*(0.197326*0.197326)/(LS*LS*spacing*spacing); //3 = 1/q, q is the quark charge
 mass_latt = (t_real)(spacing*mass_phys/197.326) + 0.0*I;
#ifdef MU
 MuPhys    = (t_real)(197.326*Mu/spacing);
#endif

 fprintf(stdout, "\n\t RUN PARAMETERS: \t\n\n");
 fprintf(stdout, "\t beta:\t\t %2.4lf\n", beta);
 fprintf(stdout, "\t spacing:\t %2.4lf fm\n", spacing);
 fprintf(stdout, "\t H (lat.u.):\t %i\n", H);
 fprintf(stdout, "\t H (ph.u.):\t %2.4lf GeV^2\n", Hphys);
#ifdef MU
 fprintf(stdout, "\t Mu (lat.u.):\t %2.4lf\n", Mu);
 fprintf(stdout, "\t Mu (ph.u.):\t %2.4lf MeV\n", MuPhys);
#endif
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
#ifdef MU
int sort_ev(t_complex* revecs, t_complex* levecs, t_complex* evals, int nev)
#else
int sort_ev(t_complex* revecs, t_complex* evals, int nev)
#endif
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
#ifdef MU
  memcpy(                            tmp, &(levecs[i*VOL*NCOLORS*NDIRAC]), VOL*NCOLORS*NDIRAC*sizeof(t_complex));
  memcpy(&(levecs[i*VOL*NCOLORS*NDIRAC]), &(levecs[k*VOL*NCOLORS*NDIRAC]), VOL*NCOLORS*NDIRAC*sizeof(t_complex));
  memcpy(&(levecs[k*VOL*NCOLORS*NDIRAC]),                             tmp, VOL*NCOLORS*NDIRAC*sizeof(t_complex));
#endif  
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

t_complex R1[NDIRAC][NDIRAC];
t_complex R2[NDIRAC][NDIRAC][NDIRAC][NDIRAC];

#ifdef MU
void init_props(t_complex *revecs, t_complex *levecs, t_complex *evals, int nev, t_complex mass)
#else
void init_props(t_complex *revecs, t_complex *evals, int nev, t_complex mass)
#endif
{
 int i, x, a, b, id1, id2, id3, id4;

 t_complex S[NCOLORS][NCOLORS][NDIRAC][NDIRAC];
 
 set_czero(&(R1[0][0]), NDIRAC*NDIRAC);
 set_czero(&(R2[0][0][0][0]), NDIRAC*NDIRAC*NDIRAC*NDIRAC);

 for(x=0; x<VOL; x++)
 {
  set_czero(&(S[0][0][0][0]), NDIRAC*NDIRAC*NCOLORS*NCOLORS);
  
  for(i=0; i<nev; i++)
   for(a=0; a<NCOLORS; a++)
    for(b=0; b<NCOLORS; b++)
     for(id1=0; id1<NDIRAC; id1++)
      for(id2=0; id2<NDIRAC; id2++)
#ifdef MU
       S[a][b][id1][id2] += CCALL(conj)(levecs[i*EV_SIZE + NCOLORS*NDIRAC*x + NDIRAC*a + id1])*revecs[i*EV_SIZE + NCOLORS*NDIRAC*x + NDIRAC*b + id2]/(evals[i] + mass);
#else
       S[a][b][id1][id2] += CCALL(conj)(revecs[i*EV_SIZE + NCOLORS*NDIRAC*x + NDIRAC*a + id1])*revecs[i*EV_SIZE + NCOLORS*NDIRAC*x + NDIRAC*b + id2]/(evals[i] + mass);
#endif
    
  for(a=0; a<NCOLORS; a++)
   for(id1=0; id1<NDIRAC; id1++)
    for(id2=0; id2<NDIRAC; id2++)
     R1[id1][id2] += S[a][a][id2][id1];
  
  for(a=0; a<NCOLORS; a++)
   for(b=0; b<NCOLORS; b++)
    for(id1=0; id1<NDIRAC; id1++)
     for(id2=0; id2<NDIRAC; id2++)
      for(id3=0; id3<NDIRAC; id3++)
       for(id4=0; id4<NDIRAC; id4++)
	R2[id1][id2][id3][id4] += S[a][a][id2][id1]*S[b][b][id4][id3] - S[b][a][id4][id1]*S[a][b][id2][id3];
 };
   
 for(id1=0; id1<NDIRAC; id1++)
  for(id2=0; id2<NDIRAC; id2++)
   R1[id1][id2] = R1[id1][id2]/(t_real)VOL;
 
 for(id1=0; id1<NDIRAC; id1++)
  for(id2=0; id2<NDIRAC; id2++)
   for(id3=0; id3<NDIRAC; id3++)
    for(id4=0; id4<NDIRAC; id4++)
     R2[id1][id2][id3][id4] = R2[id1][id2][id3][id4]/(t_real)VOL;
}

void single_op_av(int n_op_idx, t_complex *op, double* av)
{
 int i, id1, id2;    
 
 for(i=0; i<n_op_idx; i++)
  av[i] = 0.0;
 
 for(i=0; i<n_op_idx; i++)
  for(id1=0; id1<NDIRAC; id1++)
   for(id2=0; id2<NDIRAC; id2++)
    av[i] += (double)CCALL(creal)(R1[id1][id2]*op[i*NDIRAC*NDIRAC + id1*NDIRAC + id2]);
}

void double_op_av(int n_op_idx, t_complex *op1, t_complex *op2, double* av)
{
 int i, id1, id2, id3, id4;
 
 for(i=0; i<n_op_idx; i++)
  av[i] = 0.0;
 
 for(i=0; i<n_op_idx; i++)
  for(id1=0; id1<NDIRAC; id1++)
   for(id2=0; id2<NDIRAC; id2++)
    for(id3=0; id3<NDIRAC; id3++)
     for(id4=0; id4<NDIRAC; id4++)
      av[i] += (double)CCALL(creal)(R2[id1][id2][id3][id4]*op1[i*NDIRAC*NDIRAC + id1*NDIRAC + id2]*op2[i*NDIRAC*NDIRAC + id3*NDIRAC + id4]);
}

int main (int argc, char * argv[])
{ 
 int ifile, imass, i, k;
 
 t_complex* revecs;
 t_complex* revals;
#ifdef MU 
 t_complex* levecs;
#endif 
 int incl_zm = 0; // TODO: input by user
 
 //Split inames
 int nfiles = argc - 1;
 if(!nfiles)
  panic("%s: no input files specified!!!\n", __func__);
  
 initialize_gamma(); 
 read_sim_params("params.in");
 
 // Singles - current, magnetization 
#define N_SINGLE_OPERATORS 3
 value*      single_op_vals[N_SINGLE_OPERATORS][NMASSES];
 int         n_single_indices[N_SINGLE_OPERATORS];
 t_complex*  single_operators[N_SINGLE_OPERATORS];

 // Initializing single operators
 for(imass=0; imass<NMASSES; imass++)
  single_op_vals[0][imass]   = init_value(4, "j", "GeV^3", pow(0.197326/spacing, 3));
 n_single_indices[0] = 4;
 single_operators[0] = &(dgamma[0][0][0]);

 for(imass=0; imass<NMASSES; imass++)
  single_op_vals[1][imass]   = init_value(6, "s", "GeV^3", pow(0.197326/spacing, 3));
 n_single_indices[1] = 6;
 single_operators[1] = &(sigma[0][0][0]);
 
 for(imass=0; imass<NMASSES; imass++)
  single_op_vals[2][imass]   = init_value(1, "sigma", "GeV^3", pow(0.197326/spacing, 3));
 n_single_indices[2] = 1;
 single_operators[2] = &(dunity[0][0]);

 // Doubles - chirality x chirality, current x current, magnetization x magnetization, chirality x magnetization
#define N_DOUBLE_OPERATORS 5
 value*      double_op_vals[N_DOUBLE_OPERATORS][NMASSES];
 int         n_double_indices[N_DOUBLE_OPERATORS];
 t_complex*  double_operators[N_DOUBLE_OPERATORS][2];

 // Initializing double operators
 for(imass=0; imass<NMASSES; imass++)
  double_op_vals[0][imass]      = init_value(1, "r52", "GeV^6", pow(0.197326/spacing, 6));
 n_double_indices[0]    = 1;
 double_operators[0][0] = &(gamma5[0][0]);
 double_operators[0][1] = &(gamma5[0][0]);

 for(imass=0; imass<NMASSES; imass++)
  double_op_vals[1][imass]      = init_value(4, "j2", "GeV^6", pow(0.197326/spacing, 6));
 n_double_indices[1]    = 4;
 double_operators[1][0] = &(dgamma[0][0][0]);
 double_operators[1][1] = &(dgamma[0][0][0]);

 for(imass=0; imass<NMASSES; imass++)
  double_op_vals[2][imass]      = init_value(6, "s2", "GeV^6", pow(0.197326/spacing, 6));
 n_double_indices[2]    = 6;
 double_operators[2][0] = &(sigma[0][0][0]);
 double_operators[2][1] = &(sigma[0][0][0]);

 t_complex* g5tmp       = (t_complex *)malloc(6*NDIRAC*NDIRAC*sizeof(t_complex));
 for(i=0; i<6; i++)
  memcpy(&(g5tmp[i*NDIRAC*NDIRAC]), &(gamma5[0][0]), NDIRAC*NDIRAC*sizeof(t_complex));
 for(imass=0; imass<NMASSES; imass++)
  double_op_vals[3][imass]      = init_value(6, "r5s", "GeV^6", pow(0.197326/spacing, 6));
 n_double_indices[3]    = 6;
 double_operators[3][0] = &(sigma[0][0][0]);
 double_operators[3][1] = g5tmp;

 for(imass=0; imass<NMASSES; imass++)
  double_op_vals[4][imass]      = init_value(4, "r5j", "GeV^6", pow(0.197326/spacing, 6));
 n_double_indices[4]    = 4;
 double_operators[4][0] = &(dgamma[0][0][0]);
 double_operators[4][1] = g5tmp;

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
  
  
   
#ifdef MU
  res = read_ev(argv[ifile + 1], &evd, LM_SECTION);
  if((res & E_READ_CANTREAD) || (res & E_READ_CANTOPEN) || (res & E_READ_SECNOTHERE))
   panic("%s: cannot read the LM section from the file %s!!!\n", __func__, argv[ifile + 1]);
  /* Sorting eigenvalues so that left/right agree */ 
  int k, ki = 0, exchg = 0;
  t_cds_vector* tmp;
  if(posix_memalign((void*)&tmp, 16, sizeof(t_cds_vector)))
   panic("%s: can't allocate memory for tmp (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
  for(i=0; i<evd.nev-1; i++)
  {
   exchg = 0;
   t_real mindist = 100000000;
   for(k=i; k<evd.nev; k++)
   {
    t_real dist = CCALL(cabs)(revals[k] - CCALL(conj)(evd.evals[i]));
    if(dist < mindist)
    {
     mindist = dist;
     ki = k;   
    };
   };
   k = ki;
   memcpy(                               tmp, &(evd.evecs[i*VOL*NCOLORS*NDIRAC]), VOL*NCOLORS*NDIRAC*sizeof(t_complex));
   memcpy(&(evd.evecs[i*VOL*NCOLORS*NDIRAC]), &(evd.evecs[k*VOL*NCOLORS*NDIRAC]), VOL*NCOLORS*NDIRAC*sizeof(t_complex));
   memcpy(&(evd.evecs[k*VOL*NCOLORS*NDIRAC]),                                tmp, VOL*NCOLORS*NDIRAC*sizeof(t_complex));
   t_complex tmpval = evd.evals[i];
   evd.evals[i] = evd.evals[k];
   evd.evals[k] = tmpval;
#ifdef NOISY_OUTPUT
   fprintf(stdout, "Exchanging %i and %i ...\n", i, k);
   fflush(stdout);
#endif
  };
  free(tmp);
  levecs = (t_complex *)malloc(evd.nev*VOL*NCOLORS*NDIRAC*sizeof(t_complex)); 
  memcpy(levecs, &(evd.evecs[0]), evd.nev*VOL*NCOLORS*NDIRAC*sizeof(t_complex));
#endif 

#ifdef MU
  reorthogonalize(levecs, revecs, evd.nev);
#endif

  project_ev(&(revals[0]), evd.nev);
#ifdef MU
  int nzero = sort_ev(&(revecs[0]), &(levecs[0]), &(revals[0]), evd.nev); //TODO: sort_ev for complex spectrum???
#else
  int nzero = sort_ev(&(revecs[0]), &(revals[0]), evd.nev); //TODO: sort_ev for complex spectrum???
#endif
  int first_mode = (incl_zm)? 0 : nzero;
  int n_use_evals = 2*max_nz_evals + ((first_mode)? 0 : nzero);
  if(n_use_evals>evd.nev)
   panic("%s: I do not have %i evs, only %i!!!\n", __func__, n_use_evals, evd.nev);

  for(imass=0; imass<NMASSES; imass++)
  {
#ifdef MU
   init_props(&(revecs[first_mode*VOL*NCOLORS*NDIRAC]), &(levecs[first_mode*VOL*NCOLORS*NDIRAC]), &(revals[first_mode]), n_use_evals, mass_latt + 0.5*mass_latt*(imass - 1));
#else
   init_props(&(revecs[first_mode*VOL*NCOLORS*NDIRAC]), &(revals[first_mode]), n_use_evals, mass_latt + 0.5*mass_latt*(imass - 1));
#endif
   //Calculating averages of singles
   for(i=0; i<N_SINGLE_OPERATORS; i++)
   {
    double *av_tmp = (double *)malloc(n_single_indices[i]*sizeof(double));
    single_op_av(n_single_indices[i], single_operators[i], av_tmp);
    add_point(single_op_vals[i][imass], av_tmp);
    free(av_tmp);
   };

   //Calculating averages of doubles
   for(i=0; i<N_DOUBLE_OPERATORS; i++)
   {
    double *av_tmp = (double *)malloc(n_double_indices[i]*sizeof(double));
    double_op_av(n_double_indices[i], double_operators[i][0], double_operators[i][1], av_tmp);
    add_point(double_op_vals[i][imass], av_tmp);
    free(av_tmp);
   };
  };
  
  //freeing up the memory
  free_evd(&evd);
  free(revecs);
  free(revals);
#ifdef MU
  free(levecs);
#endif
 };
 
 //Saving/displaying data 
 for(i=0; i<N_SINGLE_OPERATORS; i++)
 {
  double *mX = (double *)malloc(n_single_indices[i]*sizeof(double));
  double *dX = (double *)malloc(n_single_indices[i]*sizeof(double));
  double *aX = (double *)malloc(n_single_indices[i]*sizeof(double));
  double *eX = (double *)malloc(n_single_indices[i]*sizeof(double));
  
  set_dzero(aX, n_single_indices[i]);
  set_dzero(eX, n_single_indices[i]);
  for(imass=0; imass<NMASSES; imass++)
  {
   double w = (double)(4 - NMASSES*imass)/(double)NMASSES;
   get_value(single_op_vals[i][imass], mX, dX);
   for(k=0; k<n_single_indices[i]; k++)
   {
    aX[k] += w*mX[k];
    eX[k]  = dX[k];
   };
  };
  fprintf(stdout, "%s: ", single_op_vals[i][0]->name);
  for(k=0; k<n_single_indices[i]; k++)
   fprintf(stdout, "%2.4E +/- %2.4E, ", aX[k], eX[k]);
  fprintf(stdout, "%s\n\n", single_op_vals[i][0]->units);

  char outfname[300];
  sprintf(outfname, "%s_%s.dat", single_op_vals[i][0]->name, outfname_suffix);
  FILE *outfile = fopen(outfname, "a");
  fprintf(outfile, "%+02.4lf   \t", Hphys);
#ifdef MU
  fprintf(outfile, "%+04.4lf   \t", MuPhys);
#endif
  for(k=0; k<n_single_indices[i]; k++)
   fprintf(outfile, "%+02.4E  \t%+02.4E  \t", aX[k], eX[k]);
  fprintf(outfile, "\n");
  fclose(outfile);
  
  free(mX);
  free(dX);
  free(aX);
  free(eX);
 };

 for(i=0; i<N_DOUBLE_OPERATORS; i++)
 {
  double *mX = (double *)malloc(n_double_indices[i]*sizeof(double));
  double *dX = (double *)malloc(n_double_indices[i]*sizeof(double));
  double *aX = (double *)malloc(n_double_indices[i]*sizeof(double));
  double *eX = (double *)malloc(n_double_indices[i]*sizeof(double));
  set_dzero(aX, n_double_indices[i]);
  set_dzero(eX, n_double_indices[i]);
  for(imass=0; imass<NMASSES; imass++)
  {
   double w = (double)(4 - NMASSES*imass)/(double)NMASSES;
   get_value(double_op_vals[i][imass], mX, dX);
   for(k=0; k<n_double_indices[i]; k++)
   {
    aX[k] += w*mX[k];
//   TODO: remove after debugging
//     if(i==2 && (k==0||k==5))
//      fprintf(stdout, "mX[%i] at imass=%i: %6.8E\n", k, imass, (double)mX[k]);
    eX[k] += w*w*dX[k]*dX[k];
   };
  };
  for(k=0; k<n_double_indices[i]; k++)
   eX[k] = sqrt(eX[k]);
  fprintf(stdout, "%s: ", double_op_vals[i][0]->name);
  for(k=0; k<n_double_indices[i]; k++)
   fprintf(stdout, "%2.4E +/- %2.4E, ", aX[k], eX[k]);
  fprintf(stdout, "%s\n\n", double_op_vals[i][0]->units);

  char outfname[300];
  sprintf(outfname, "%s_%s.dat", double_op_vals[i][0]->name, outfname_suffix);
  FILE *outfile = fopen(outfname, "a");
  fprintf(outfile, "%+02.4lf   \t", Hphys);
#ifdef MU
  fprintf(outfile, "%+04.4lf   \t", MuPhys);
#endif
  for(k=0; k<n_double_indices[i]; k++)
   fprintf(outfile, "%+02.4E  \t%+02.4E  \t", aX[k], eX[k]);
  fprintf(outfile, "\n");
  fclose(outfile);
  
  free(mX);
  free(dX);
  free(aX);
  free(eX);
 };
 
 free(g5tmp);

 //freeing up the memory
 for(i=0; i<N_SINGLE_OPERATORS; i++)
  for(imass=0; imass<NMASSES; imass++)
   free_value(single_op_vals[i][imass]);
 for(i=0; i<N_DOUBLE_OPERATORS; i++)
  for(imass=0; imass<NMASSES; imass++)
   free_value(double_op_vals[i][imass]);
 
 return EXIT_SUCCESS;
}
