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
  A[i] = 0.0;
}

char      outfname_suffix[300];

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

int main(int argc, char *argv[]) // format of input: meson prefix suffix datafiles
{ 
 t_op_ev evd;
 init_evd(&evd);
 
 int i, x, a, id;
 //name of the input file
 char infilename[300];
 char outfilename[300];
 
 if(argc!=3)
 {
  fprintf(stderr, "ERROR: wrong argument!!!\n");
  fflush(stderr);
  return E[M ]:XIT_FAILURE;
 };
  
 strcpy(infilename,  argv[1]);
 strcpy(outfilename, argv[2]);
 fprintf(stdout, "Input file: %s\n", infilename);
 fprintf(stdout, "Output file: %s\n", outfilename);
 
 
 //initialize the structure which contains Dirac eigenvectors and eigenval//Read the data file: eigenvectors and eigenvalues are loaded in the arrays evd.evecs and evd.evals, respectively
 //evd.nev is the total number of eigenvalue + eigenvector pairs
 int res = read_ev(infilename, &evd, SM_SECTION);
 
 if((res & E_READ_CANTREAD) || (res & E_READ_CANTOPEN) || (res & E_READ_SECNOTHERE))
  panic("%s: cannot read the file %s!!!\n", __func__, infilename);
 
 //loop over all eigenvalues
 //CCALL is a macro for calling type-dependent complex number routines
 for(i=0; i<evd.nev; i++)
  fprintf(stdout, "%2.6lf %2.6lf\n", CCALL(creal)(evd.evals[i]), CCALL(cimag)(evd.evals[i]));
 
 int xc[DIM];
 
 xc[2] = (LS - 1)/2;
 xc[3] = (LS - 1)/2;
 
 FILE *outfile = fopen(outfilename, "w");
 
 for(xc[0] = 0; xc[0]<LS; xc[0]++)
 {
  for(xc[1] = 0; xc[1]<LS; xc[1]++)
  {
   t_complex rho = 0.0 + I*0.0;
   x = INDEX_X_4D(xc);
   for(a=0; a<NCOLORS; a++)
    for(id=0; id<NDIRAC; id++)
     rho += CCALL(conj)(evd.evecs[i*EV_SIZE + x*NCOLORS*NDIRAC + NDIRAC*a + id])*evd.evecs[i*EV_SIZE + x*NCOLORS*NDIRAC + NDIRAC*a + id];
    fprintf(outfile, "%4.6lf\t", (double)CCALL(creal)(rho));     
    fprintf(stdout,  "%4.6lf\t", (double)CCALL(creal)(rho));     
  };
  fprintf(outfile, "\n");
  fprintf(stdout,  "\n");
  fflush(stdout);
 };
  
 
 fclose(outfile);
  
 
 //evd.evecs[i*EV_SIZE + x*NCOLORS*NDIRAC + NDIRAC*a + id];
 //i -  no. of eigenvalues
 //x  - spatial coordinate, = 0 ... VOL - 1
 //a  - color index
 //id - Dirac index
 
 //t_complex *X; X = (t_complex* )malloc(EV_SIZE*sizeof(t_complex));
 
 //malloc();
 
 //free();  
 
  
 free_evd(&evd);
 
 return EXIT_SUCCESS;
}
