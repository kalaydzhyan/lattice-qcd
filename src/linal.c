/*{{{*/
/*!
 * \file linal.c
 *
 * \brief
 *
 *
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \date   ðÎÄ ïËÔ 25 14:14:11 MSD 2004
 */
/*}}}*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "panic.h"
#include "linal.h"
#include "lattice.h"

void axpby(t_complex *x, const t_real a, const t_complex* y, const t_real b)
{
#if (EXP_SSE && SINGLE_PREC)
    __asm__ __volatile__ (
            "movss  %3, %%xmm6              /* a -> xmm6 */\n\t"
            "shufps $0, %%xmm6, %%xmm6\n\t"
            "movss  %4, %%xmm7              /* b -> xmm7 */\n\t"
            "shufps $0, %%xmm7, %%xmm7\n\t"
            "movl   %0, %%esi               /* x -> esi */\n\t"
            "movl   %2, %%edi               /* y -> edi */\n\t"
            "mov    %5, %%ecx\n\t"
            ".LOOP_AXPBY:\n\t"
            "movaps (%%esi), %%xmm0\n\t"
            "movaps 16(%%esi), %%xmm1\n\t"
            "mulps  %%xmm6, %%xmm0\n\t"
            "mulps  %%xmm6, %%xmm1\n\t"
            "movaps (%%edi), %%xmm2\n\t"
            "movaps 16(%%edi), %%xmm3\n\t"
            "mulps  %%xmm7, %%xmm2\n\t"
            "mulps  %%xmm7, %%xmm3\n\t"
            "addps  %%xmm2, %%xmm0\n\t"
            "addps  %%xmm3, %%xmm1\n\t"
            "movaps %%xmm0, (%%esi)\n\t"
            "movaps %%xmm1, 16(%%esi)\n\t"
            "addl   $32, %%esi\n\t"
            "addl   $32, %%edi\n\t"
            "subl   $1, %%ecx\n\t"
            "jnz    .LOOP_AXPBY"
            : "=m" (x)
            : "m" (x), "m" (y), "m" (a), "m" (b), "i" (VOL*NDIRAC*NCOLORS/4)
            : "%ecx", "%esi", "%edi", "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm6", "%xmm7"
            );
#else
    int i;

    for (i = 0; i < VOL*NDIRAC*NCOLORS; i++)
        x[i] = a*x[i] + b*y[i];
#endif
}

void axpbyn(t_complex *x, const t_real a, const t_complex* y, const t_real b, int n)
{
    int i;
    for(i = 0; i < n; i++)
     x[i] = a*x[i] + b*y[i];
}

void axpbynd(t_double_complex *x, const t_double_real a, const t_double_complex* y, const t_double_real b, int n)
{
    int i;
    for(i = 0; i < n; i++)
     x[i] = a*x[i] + b*y[i];
}

inline void xcpy(t_complex* x, const t_complex* y)
{
    int i;

    for (i = 0; i < VOL*NDIRAC*NCOLORS; i++)
        x[i] = y[i];
}

void xpby(t_complex* x, const t_complex* y, const t_real b)
{
#if (EXP_SSE && SINGLE_PREC)
    __asm__ __volatile__ (
            "movss  %3, %%xmm7              /* b -> xmm7 */\n\t"
            "shufps $0, %%xmm7, %%xmm7\n\t"
            "movl   %0, %%esi               /* x -> esi */\n\t"
            "movl   %2, %%edi               /* y -> edi */\n\t"
            "mov    %4, %%ecx\n\t"
            ".LOOP_XPBY:\n\t"
            "movaps (%%esi), %%xmm0\n\t"
            "movaps 16(%%esi), %%xmm1\n\t"
            "movaps (%%edi), %%xmm2\n\t"
            "movaps 16(%%edi), %%xmm3\n\t"
            "mulps  %%xmm7, %%xmm2\n\t"
            "mulps  %%xmm7, %%xmm3\n\t"
            "addps  %%xmm2, %%xmm0\n\t"
            "addps  %%xmm3, %%xmm1\n\t"
            "movaps %%xmm0, (%%esi)\n\t"
            "movaps %%xmm1, 16(%%esi)\n\t"
            "addl   $32, %%esi\n\t"
            "addl   $32, %%edi\n\t"
            "subl   $1, %%ecx\n\t"
            "jnz    .LOOP_XPBY"
            : "=m" (x)
            : "m" (x), "m" (y), "m" (b), "i" (VOL*NDIRAC*NCOLORS/4)
            : "%ecx", "%esi", "%edi", "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm7"
            );
#else
    int i;

    for (i = 0; i < VOL*NDIRAC*NCOLORS; i++)
        x[i] += b*y[i];
#endif
}

inline void xpcby(t_complex* x, const t_complex* y, const t_complex b)
{
    int i;

    for (i = 0; i < VOL*NDIRAC*NCOLORS; i++)
        x[i] += b*y[i];
}


inline void ag5xpby(t_complex* x, const t_real a, const t_complex* y, const t_real b)
{
    int i;

    for (i = 0; i < VOL*NCOLORS*NDIRAC; i += NDIRAC) {
        x[i+0] =  a*x[i+0] + b*y[i+0];
        x[i+1] =  a*x[i+1] + b*y[i+1];
        x[i+2] = -a*x[i+2] + b*y[i+2];
        x[i+3] = -a*x[i+3] + b*y[i+3];
    }
}

inline void ag5xpcby(t_complex* x, const t_complex a, const t_complex* y, const t_complex b)
{
    int i;

    for (i = 0; i < VOL*NCOLORS*NDIRAC; i += NDIRAC) {
        x[i+0] =  a*x[i+0] + b*y[i+0];
        x[i+1] =  a*x[i+1] + b*y[i+1];
        x[i+2] = -a*x[i+2] + b*y[i+2];
        x[i+3] = -a*x[i+3] + b*y[i+3];
    }
}

inline void xbyg5(t_complex* x)
{
    int i;

    for (i = 0; i < VOL*NCOLORS*NDIRAC; i += NDIRAC) {
        x[i+0] =  x[i+0];
        x[i+1] =  x[i+1];
        x[i+2] =  x[i+2];
        x[i+3] =  x[i+3];
    }
}

inline void g5x(t_complex* x, t_complex* y)
{
    int i;

    for (i = 0; i < VOL*NCOLORS*NDIRAC; i += NDIRAC) {
        y[i+0] =  x[i+0];
        y[i+1] =  x[i+1];
        y[i+2] = -x[i+2];
        y[i+3] = -x[i+3];
    }
}

/*!
 * \brief Evaluate Hermitian inner product
 *
 * \f[ result = \sum_{a=0}^{V} x^*_a y_a,\ \ V = NDIRAC*NCOLORS*VOL. \f]
 *
 * \param  t_complex *x
 * \param  t_complex *y
 *
 * \return t_complex
 */
inline t_complex innprod(t_complex* x, t_complex* y)
{
    int i;
    double complex sum = 0.0 + I*0.0;

    for (i = 0; i < VOL*NCOLORS*NDIRAC; i++) {
        sum += CCALL(conj)(x[i])*y[i];
    }

    return (t_complex) sum;
}

inline t_real vnorm(t_complex* x)
{
 int i;
 double complex sum = 0.0 + I*0.0;
 for(i=0; i<VOL*NCOLORS*NDIRAC; i++)
 {
  sum += CCALL(conj)(x[i])*x[i];
 };
 return (t_real)sqrt(creal(sum));
}

inline t_real    adiffn(t_complex* x, t_complex* y, int n)
{
 int i;
 t_real sum = 0.0;
 for(i=0; i<n; i++)
  sum += CCALL(creal)(x[i] - y[i])*CCALL(creal)(x[i] - y[i]) + CCALL(cimag)(x[i] - y[i])*CCALL(cimag)(x[i] - y[i]);
 return sqrt(sum);
}

t_double_real    adiffnd(t_double_complex* x, t_double_complex* y, int n)
{
 int i;
 t_double_real sum = 0.0;
 for(i=0; i<n; i++)
#ifdef LONG_INTERNAL_ARITHMETIC
  sum += creall(x[i] - y[i])*creall(x[i] - y[i]) + cimagl(x[i] - y[i])*cimagl(x[i] - y[i]);
#else
  sum += creal(x[i] - y[i])*creal(x[i] - y[i]) + cimag(x[i] - y[i])*cimag(x[i] - y[i]);
#endif
 return sqrt(sum);
}

void  ax(t_complex *x, const t_real a)
{
 int i;
 for(i=0; i<VOL*NCOLORS*NDIRAC; i++)
  x[i] *= a;
}

void axc(t_complex *x, const t_complex a)
{
 int i;
 for(i=0; i<VOL*NCOLORS*NDIRAC; i++)
  x[i] *= a;
}

void set_zero(t_complex* x)
{
 int i;
 for(i=0; i<VOL*NCOLORS*NDIRAC; i++)
  x[i] = 0.0 + 0.0*I;
}

void lu_factorize(t_double_complex *A, t_double_complex *L, t_double_complex *U, int size) //Performs LU factorization of (size x size) matrix A
{
 int i, j, k, q;
 t_double_complex S1, S2;
 if(A[0]==0)
  panic("%s: can't perform LU decomposition", __func__);
 for(i=0; i<size; i++)
  L[size*i] = A[size*i];
 for(j=0; j<size; j++)
  U[j] = A[j]/L[0];
 for(q=1; q<size; q++)
 {
  //L elements
  for(i=q; i<size; i++)
  {
   S1 = 0.0 + 0.0*I;
   // The sum
   for(j=0; j<q; j++)
    S1 += L[size*i + j]*U[size*j + q];
   // L elements
   L[size*i + q] = A[size*i + q] - S1;
  };
  //U elements
  for(j=q; j<size; j++)
  {
   S2 = 0.0 + 0.0*I;
   //The sum
   for(k=0; k<q; k++)
    S2 += L[size*q + k]*U[size*k + j];
   // U elements
   if(L[size*q + q]==0)
    panic("%s: can't perform LU decomposition", __func__);
   U[size*q + j] = (A[size*q + j] - S2)/L[size*q + q];
  };
 };
}

void l_solve(t_double_complex *L, t_double_complex *F, t_double_complex* R, int size)
{
 int i, j;
 t_double_complex Q;
 R[0] = F[0]/L[0];
 for(i=0; i<size; i++)
 {
  Q = F[i];
  for(j=0; j<i; j++)
   Q -= L[size*i + j]*R[j];
  R[i] = Q/L[size*i + i];
 };
}

void u_solve(t_double_complex *U, t_double_complex *F, t_double_complex* R, int size)
{
 int i, j;
 t_double_complex Q;
 R[size-1] = F[size-1];
 for(i = size-2; i>=0; i--)
 {
  Q = F[i];
  for(j = i+1; j<size; j++)
   Q -= U[size*i + j]*R[j];
  R[i] = Q;
 };
}

void lu_solve(t_double_complex *L, t_double_complex *U, t_double_complex *F, t_double_complex* R, int size)
{
 t_double_complex *Y = (t_double_complex *)malloc(size*sizeof(t_double_complex));
 l_solve(L, F, Y, size);
 u_solve(U, Y, R, size);
 free(Y);
}

void inverse(t_double_complex *A, t_double_complex *R, int size) //Returns the inverse of A in R
{
 int i, j;
 t_double_complex *L = (t_double_complex *)malloc(size*size*sizeof(t_double_complex));
 t_double_complex *U = (t_double_complex *)malloc(size*size*sizeof(t_double_complex));
 t_double_complex *Q = (t_double_complex *)malloc(size*sizeof(t_double_complex));
 t_double_complex *D = (t_double_complex *)malloc(size*sizeof(t_double_complex));

 for(i=0; i<size; i++)
 {
  Q[i] = 0.0 + 0.0*I;
  D[i] = 0.0 + 0.0*I;
 };

 lu_factorize(A, L, U, size);

 for(i=0; i<size; i++)
 {
  if(i>0)
  {
   D[i-1] = 0.0 + 0.0*I;
   D[i]   = 1.0 + 0.0*I;
  }
  else
   D[i] = 1.0 + 0.0*I;

  lu_solve(L, U, D, Q, size);
  for(j=0; j<size; j++)
   R[size*j + i] = Q[j];
 };
 free(L);
 free(U);
 free(Q);
 free(D);
}

void reorthogonalize(t_complex *L, t_complex *R, int n_vecs)
{
    int i, j;
    t_double_complex *metrics  = (t_double_complex *)malloc(n_vecs*n_vecs*sizeof(t_double_complex));
    t_double_complex *imetrics = (t_double_complex *)malloc(n_vecs*n_vecs*sizeof(t_double_complex));

    for(i = 0; i < n_vecs; i++)
     for(j = 0; j < n_vecs; j++)
      metrics[n_vecs*i + j] = (t_double_complex)innprod(&(L[i*VOL*NCOLORS*NDIRAC]),&(R[j*VOL*NCOLORS*NDIRAC]));

    inverse(metrics, imetrics, n_vecs);

    t_complex* tmp = (t_complex *)malloc(n_vecs*VOL*NCOLORS*NDIRAC*sizeof(t_complex));
    for(i = 0; i < n_vecs; i++)
    {
     set_zero(&(tmp[i*VOL*NCOLORS*NDIRAC]));
     for(j  = 0; j < n_vecs; j++)
      xpcby(&(tmp[i*VOL*NCOLORS*NDIRAC]), &(R[j*VOL*NCOLORS*NDIRAC]), (t_complex)imetrics[j*n_vecs + i]);
    };
    memcpy(&(R[0]), &(tmp[0]), n_vecs*VOL*NCOLORS*NDIRAC*sizeof(t_complex));

    free(tmp);
    free(metrics);
    free(imetrics);
}

void reorthogonalize_gs(t_complex *vecs, int n_vecs)
{
 //TODO: gram-schmidt
}

#define EPS_ORTHO       10e-6   //!< allowed nonorthogonality
#define EPS_EVDIFF      10e-6   //!< allowed \f$|\lambda_i - \lambda_j| \f$ < EPS_EVDIFF
/*{{{*/
/*!
 * \brief
 *
 * \param nev
 * \param vol
 * \param evals
 * \param evecs
 *
 * \return void
 */
/*}}}*/
void check_ortho(int nev, int vol, t_complex* evals, t_complex* evecs) { /*{{{*/
    int i, j;
    t_complex ip;

    for (i = 0; i < nev; i++) {
        for (j = 0; j < nev; j++) {
            if (i == j)
                continue;
            ip = innprod(&(evecs[i*vol]), &(evecs[j*vol]));
            if ((fabs(CCALL(creal)(evals[i]) - CCALL(creal)(evals[j])) > EPS_EVDIFF &&
                    fabs(CCALL(cimag)(evals[i]) - CCALL(cimag)(evals[j])) > EPS_EVDIFF) &&
                    sqrt(CCALL(creal)(ip)*CCALL(creal)(ip)+CCALL(cimag)(ip)*CCALL(cimag)(ip)) > EPS_ORTHO)
                warn("%s: eigenvectors (lambda[%i]=%8.6f+I %8.6f) and (lambda[%i]=%8.6f+I %8.6f) not orthogonal (ip = %8.6f+I %8.6f",
                        __func__, i, CCALL(creal)(evals[i]), CCALL(cimag)(evals[i]),
                        j, CCALL(creal)(evals[j]), CCALL(cimag)(evals[j]),
                        CCALL(creal)(ip), CCALL(cimag)(ip));
        }
    }
}/*}}}*/

void init_CDSV(t_cds_vector* csdv, int* x, int c, int d) /*{{{*/
{
    int idx;

    for (idx = 0; idx < VOL*NDIRAC*NCOLORS; idx++)
        ((t_complex*)csdv)[idx] = 0.0;
    idx = INDEX_X(x);
    (*csdv)[idx][c][d] = 1.0;
}/*}}}*/

