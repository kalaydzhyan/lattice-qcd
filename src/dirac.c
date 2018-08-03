/*{{{1*/
/*!
 * \file dirac.c
 *
 * \brief
 *
 *
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \author Pavel Buividovich, email: gbuividovich@gmail.com (implemented background magnetic field, chemical potential, SU(3) gauge group in 2008 - 2009)
 */
/*}}}1*/
#define _XOPEN_SOURCE   600
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "defs.h"
#include "types.h"
#include "lattice.h"
#include "arpack.h"
#include "panic.h"
#include "minmax.h"
#include "linal.h"
#include "dirac.h"
#include "timing.h"
#include "arn_log.h"

#ifdef SU2
#include "su2math.h"
#endif

#ifdef SU3
#include "su3math.h"
#endif

/*{{{*/
/*!
 * \brief Hopping parameter r
 *
 * For massless Wilson Dirac Fermions \f$r = \frac{1}{8\kappa}\f$
 *//*}}}*/
t_real hopping_r = 1.0;

/*{{{*/
/*!
 * \brief \f$ \kappa = \frac{1}{-2\rho + 8 r}\f$
 *
 * \param rho
 *
 * \return t_real \f$\kappa\f$
 *//*}}}*/
t_real rho2kappa(const t_real rho)/*{{{*/
{
    return (1.0/(-2.0*rho + 8.0*hopping_r));
}/*}}}*/

/*{{{*/
/*!
 * \brief \f$ \rho = 4 r - \frac{1}{2\kappa} \f$
 *
 * \param kappa
 *
 * \return t_real \f$ \rho \f$
 *//*}}}*/
t_real kappa2rho(const t_real kappa)/*{{{*/
{
    return (4.0*hopping_r - 1.0/(2.0*kappa));
}/*}}}*/

#if (DIM == 4)
/*{{{1*/
/*! \brief \f$ \gamma_\mu\f$ in chiral representation
 *
 * \f[
 * \gamma_0 = \left(
 * \begin{tabular}{cc} $0$ & $-i\sigma_1$ \\ $i\sigma_10$ & $0$ \\ \end{tabular}
 * \right)
 * \f]
 * \f[
 * \gamma_1 = \left(
 * \begin{tabular}{cc} $0$ & $-i\sigma_2$ \\ $i\sigma_2$ & $0$ \\ \end{tabular}
 * \right)
 * \f]
 * \f[
 * \gamma_2 = \left(
 * \begin{tabular}{cc} $0$ & $-i\sigma_3$ \\ $i\sigma_3$ & $0$ \\ \end{tabular}
 * \right)
 * \f]
 * \f[
 * \gamma_3 = \left(\begin{tabular}{cc} $0$ & $1$ \\ $1$ & $0$ \\ \end{tabular}\right)
 * \f],
 * where \f$\sigma_i,\ \ \ i=1,2,3\f$ are Pauli matrices
 */
/*}}}1*/
t_complex d_gamma[DIM][NDIRAC][NDIRAC] = /*{{{*/
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
};/*}}}*/
/*{{{1*/
/*! \brief \f$ \gamma_5 \f$ in chiral representation
 *
 * \f[
 * \gamma_5 = \left(\begin{tabular}{cc} $1$ & $0$ \\ $0$ & $-1$ \\ \end{tabular}\right)
 * \f]
 */
/*}}}1*/
t_complex d_gamma5[NDIRAC][NDIRAC] = /*{{{*/
{{1.0, 0.0, 0.0, 0.0},
 {0.0, 1.0, 0.0, 0.0},
 {0.0, 0.0,-1.0, 0.0},
 {0.0, 0.0, 0.0,-1.0}
};/*}}}*/
#endif

/*{{{1*/
/*! \brief  Multiply ColorDirac vector by \f$(\pm \gamma_\mu+r)\f$
 *
 * \param mu    \f$ \gamma_\mu, \ \ \ \mu\f$, index
 * \param fb    fb = 1(forward) \f$ -\gamma_\mu \f$,
 * fb = 0(backward) \f$ \gamma_\mu \f$
 * \param out   out += \f$ (\pm \gamma_\mu + r) \f$ in
 * \param in    ColorDirac vector
 *
 * \return void
 */
/*}}}1*/

static inline void gpr_mul_CDSV(int mu, int fb, t_complex * out, t_complex * in) /*{{{*/
{
    int a;          /* color index */
    t_real gsign;   /* +/- \gamma_\mu */

    gsign = (fb == BWD ? 1.0 : -1.0);
    /*  + r \delta_{\alpha\beta} \Psi_\beta  */
#if (EXP_SSE && SINGLE_PREC)
    /*{{{*/
    __asm__ __volatile__ (
            "movl   %0, %%ecx\n\t"
            "movl   %1, %%esi\n\t"
            "movl   %2, %%edi\n\t"
            ".LOOP_M_HOPPING_R:\n\t"
            "movups (%%esi), %%xmm0\n\t"
            "movups (%%edi), %%xmm1\n\t"
            "movss  %3, %%xmm2\n\t /* hopping_r -> xmm2 */"
            "shufps $0, %%xmm2, %%xmm2\n\t"
            "mulps  %%xmm2, %%xmm1\n\t"
            "addps  %%xmm1, %%xmm0\n\t"
            "movups %%xmm0, (%%esi)\n\t"
            "addl   $16, %%esi\n\t"
            "addl   $16, %%edi\n\t"
            "subl   $1, %%ecx\n\t"
            "jnz    .LOOP_M_HOPPING_R"
            :
            : "i" (NDIRAC*NCOLORS/2), "g" (out), "g" (in), "g" (hopping_r)
            : "%ecx", "%esi", "%edi", "%xmm0", "%xmm1", "%xmm2"
        );
    /*}}}*/
#else
    for (a = 0; a < NCOLORS*NDIRAC; a++)
        out[a] += hopping_r * (in[a]);
#endif
    /* + (\gamma_\mu)_{\alpha\beta} \Psi_\beta */
    for (a = 0; a < NCOLORS; a++) { /*{{{*/
        switch (mu) /*{{{*/
        {
            case 0:
                out[a*NDIRAC]   += gsign*in[a*NDIRAC + 3]*(-I);
                out[a*NDIRAC+1] += gsign*in[a*NDIRAC + 2]*(-I);
                out[a*NDIRAC+2] += gsign*in[a*NDIRAC + 1]*(+I);
                out[a*NDIRAC+3] += gsign*in[a*NDIRAC + 0]*(+I);
                break;
            case 1:
                out[a*NDIRAC]   += gsign*in[a*NDIRAC + 3]*(-1.0);
                out[a*NDIRAC+1] += gsign*in[a*NDIRAC + 2];
                out[a*NDIRAC+2] += gsign*in[a*NDIRAC + 1];
                out[a*NDIRAC+3] += gsign*in[a*NDIRAC + 0]*(-1.0);
                break;
            case 2:
                out[a*NDIRAC]   += gsign*in[a*NDIRAC + 2]*(-I);
                out[a*NDIRAC+1] += gsign*in[a*NDIRAC + 3]*(+I);
                out[a*NDIRAC+2] += gsign*in[a*NDIRAC + 0]*(+I);
                out[a*NDIRAC+3] += gsign*in[a*NDIRAC + 1]*(-I);
                break;
            case 3:
                out[a*NDIRAC]   += gsign*in[a*NDIRAC + 2];
                out[a*NDIRAC+1] += gsign*in[a*NDIRAC + 3];
                out[a*NDIRAC+2] += gsign*in[a*NDIRAC + 0];
                out[a*NDIRAC+3] += gsign*in[a*NDIRAC + 1];
                break;
        }/*}}}*/
    }/*}}}*/
}/*}}}*/

#if (EXP_SSE && SINGLE_PREC)

#define gpr_mul_CDSV_fwd(_mu, _out, _in)/*{{{*/                     \
{                                                                   \
    int (_a);           /* color index */                           \
    /*  + r \delta_{\alpha\beta} \Psi_\beta  */                     \
    __asm__ __volatile__ (                                          \
        "movl   %0, %%ecx\n\t"                                      \
        "movl   %1, %%esi\n\t"                                      \
        "movl   %2, %%edi\n\t"                                      \
        ".LOOP_M_HOPPING_RF:\n\t"                                   \
        "movaps (%%esi), %%xmm0\n\t"                                \
        "movaps (%%edi), %%xmm1\n\t"                                \
        "movss  %3, %%xmm2\n\t /* hopping_r -> xmm2 */"             \
        "shufps $0, %%xmm2, %%xmm2\n\t"                             \
        "mulps  %%xmm2, %%xmm1\n\t"                                 \
        "addps  %%xmm1, %%xmm0\n\t"                                 \
        "movaps %%xmm0, (%%esi)\n\t"                                \
        "addl   $16, %%esi\n\t"                                     \
        "addl   $16, %%edi\n\t"                                     \
        "subl   $1, %%ecx\n\t"                                      \
        "jnz    .LOOP_M_HOPPING_RF"                                 \
        :                                                           \
        : "i" (NDIRAC*NCOLORS/2), "g" ((_out)), "g" ((_in)),        \
         "g" (hopping_r)                                            \
        : "%ecx", "%esi", "%edi", "%xmm0", "%xmm1", "%xmm2"         \
    );                                                              \
    /* - (\gamma_\mu)_{\alpha\beta} \Psi_\beta */                   \
    for ((_a) = 0; (_a) < NCOLORS; (_a)++) {                        \
        switch ((_mu)) {                                            \
            case 0:                                                 \
                (_out)[_a*NDIRAC]   -= (_in)[_a*NDIRAC + 3]*(-I);   \
                (_out)[_a*NDIRAC+1] -= (_in)[_a*NDIRAC + 2]*(-I);   \
                (_out)[_a*NDIRAC+2] -= (_in)[_a*NDIRAC + 1]*(+I);   \
                (_out)[_a*NDIRAC+3] -= (_in)[_a*NDIRAC + 0]*(+I);   \
                break;                                              \
            case 1:                                                 \
                (_out)[_a*NDIRAC]   -= (_in)[_a*NDIRAC + 3]*(-1.0); \
                (_out)[_a*NDIRAC+1] -= (_in)[_a*NDIRAC + 2];        \
                (_out)[_a*NDIRAC+2] -= (_in)[_a*NDIRAC + 1];        \
                (_out)[_a*NDIRAC+3] -= (_in)[_a*NDIRAC + 0]*(-1.0); \
                break;                                              \
            case 2:                                                 \
                (_out)[_a*NDIRAC]   -= (_in)[_a*NDIRAC + 2]*(-I);   \
                (_out)[_a*NDIRAC+1] -= (_in)[_a*NDIRAC + 3]*(+I);   \
                (_out)[_a*NDIRAC+2] -= (_in)[_a*NDIRAC + 0]*(+I);   \
                (_out)[_a*NDIRAC+3] -= (_in)[_a*NDIRAC + 1]*(-I);   \
                break;                                              \
            case 3:                                                 \
                (_out)[_a*NDIRAC]   -= (_in)[_a*NDIRAC + 2];        \
                (_out)[_a*NDIRAC+1] -= (_in)[_a*NDIRAC + 3];        \
                (_out)[_a*NDIRAC+2] -= (_in)[_a*NDIRAC + 0];        \
                (_out)[_a*NDIRAC+3] -= (_in)[_a*NDIRAC + 1];        \
                break;                                              \
        }                                                           \
    }                                                               \
}/*}}}*/

#define gpr_mul_CDSV_bwd(_mu, _out, _in)/*{{{*/                     \
{                                                                   \
    int (_a);           /* color index */                           \
    /*  + r \delta_{\alpha\beta} \Psi_\beta  */                     \
    __asm__ __volatile__ (                                          \
        "movl   %0, %%ecx\n\t"                                      \
        "movl   %1, %%esi\n\t"                                      \
        "movl   %2, %%edi\n\t"                                      \
        ".LOOP_M_HOPPING_RB:\n\t"                                   \
        "movaps (%%esi), %%xmm0\n\t"                                \
        "movaps (%%edi), %%xmm1\n\t"                                \
        "movss  %3, %%xmm2\n\t /* hopping_r -> xmm2 */"             \
        "shufps $0, %%xmm2, %%xmm2\n\t"                             \
        "mulps  %%xmm2, %%xmm1\n\t"                                 \
        "addps  %%xmm1, %%xmm0\n\t"                                 \
        "movaps %%xmm0, (%%esi)\n\t"                                \
        "addl   $16, %%esi\n\t"                                     \
        "addl   $16, %%edi\n\t"                                     \
        "subl   $1, %%ecx\n\t"                                      \
        "jnz    .LOOP_M_HOPPING_RB"                                 \
        :                                                           \
        : "i" (NDIRAC*NCOLORS/2), "g" ((_out)), "g" ((_in)),        \
         "g" (hopping_r)                                            \
        : "%ecx", "%esi", "%edi", "%xmm0", "%xmm1", "%xmm2"         \
    );                                                              \
    /* - (\gamma_\mu)_{\alpha\beta} \Psi_\beta */                   \
    for ((_a) = 0; (_a) < NCOLORS; (_a)++) {                        \
        switch ((_mu)) {                                            \
            case 0:                                                 \
                (_out)[_a*NDIRAC]   += (_in)[_a*NDIRAC + 3]*(-I);   \
                (_out)[_a*NDIRAC+1] += (_in)[_a*NDIRAC + 2]*(-I);   \
                (_out)[_a*NDIRAC+2] += (_in)[_a*NDIRAC + 1]*(+I);   \
                (_out)[_a*NDIRAC+3] += (_in)[_a*NDIRAC + 0]*(+I);   \
                break;                                              \
            case 1:                                                 \
                (_out)[_a*NDIRAC]   += (_in)[_a*NDIRAC + 3]*(-1.0); \
                (_out)[_a*NDIRAC+1] += (_in)[_a*NDIRAC + 2];        \
                (_out)[_a*NDIRAC+2] += (_in)[_a*NDIRAC + 1];        \
                (_out)[_a*NDIRAC+3] += (_in)[_a*NDIRAC + 0]*(-1.0); \
                break;                                              \
            case 2:                                                 \
                (_out)[_a*NDIRAC]   += (_in)[_a*NDIRAC + 2]*(-I);   \
                (_out)[_a*NDIRAC+1] += (_in)[_a*NDIRAC + 3]*(+I);   \
                (_out)[_a*NDIRAC+2] += (_in)[_a*NDIRAC + 0]*(+I);   \
                (_out)[_a*NDIRAC+3] += (_in)[_a*NDIRAC + 1]*(-I);   \
                break;                                              \
            case 3:                                                 \
                (_out)[_a*NDIRAC]   += (_in)[_a*NDIRAC + 2];        \
                (_out)[_a*NDIRAC+1] += (_in)[_a*NDIRAC + 3];        \
                (_out)[_a*NDIRAC+2] += (_in)[_a*NDIRAC + 0];        \
                (_out)[_a*NDIRAC+3] += (_in)[_a*NDIRAC + 1];        \
                break;                                              \
        }                                                           \
    }                                                               \
}/*}}}*/

#endif

/*{{{1*/
/*!
 * \brief Evaluate \f$\Psi_{out} = Q \Psi_{in} \f$
 *
 * We present here some theoretical background about Wilson Dirac operator.
 * Free massive Wilson Dirac operator is defined in a following way:
 * \f[
 * (D_w)_{x,y} =
 * \frac{1}{2}\sum_\mu\gamma_\mu
 * \left(\partial^f_{\mu,x,y} + \partial^b_{\mu,x,y}\right) +
 * m \delta_{x,y} -
 * \frac{ra}{2}\sum_{\mu,a}\partial^f_{\mu,x,a}\partial^b_{\mu,a,y},
 * \f]
 * where
 * \f[
 * \partial^f_{\mu,x,y} = \frac{1}{a}(\delta_{x+\hat{\mu},y} - \delta_{x,y}),\;
 * \partial^b_{\mu,x,y} = \frac{1}{a}(\delta_{x,y} - \delta_{x-\hat{\mu},y}).
 * \f]
 * About parameter \f$ r \f$ you can read in Montvay, Munster
 * "Quantum Fields on a Lattice". It's important to say that
 * \f[
 * 0 < r \le 1.
 * \f]
 * After trivial algebra one obtains
 * \f[
 * (D_w)_{x,y} = \frac{am + 4 r}{a} \delta_{x,y} - \frac{1}{2a}\sum_\mu\left[
 * (r - \gamma_\mu) \delta_{x+\hat{\mu},y} +
 * (r + \gamma_\mu) \delta_{x-\hat{\mu},y}
 * \right].
 * \f]
 * For non-free Wilson Dirac operator we should write
 * \f[
 * (D_w)_{x,y} = \frac{am + 4 r}{a} \delta_{x,y} - \frac{1}{2a}\sum_\mu\left[
 * (r - \gamma_\mu) \delta_{x+\hat{\mu},y}U_{x,\mu} +
 * (r + \gamma_\mu) \delta_{x-\hat{\mu},y}U^\dagger_{x-\hat{\mu}, \mu}
 * \right].
 * \f]
 * This expression can also be written as
 * \f[
 * {D_w}_{x,y} = \frac{1}{2a\kappa} Q_{x,y},
 * \f]
 * where
 * \f[
 * \kappa = \frac{1}{2am+8r},
 * \f]
 * \f[
 * Q_{x,y} = I_{x,y} - \kappa M_{x,y},
 * \f]
 * \f[
 * I_{x,y} = \delta_{x,y},
 * \f]
 * \f[
 * M_{x,y} = \sum_{\mu=1}^4 \left[
 * (r-\gamma_\mu)\delta_{x+\hat{\mu},y}U_{x,\mu} +
 * (r+\gamma_\mu)\delta_{x-\hat{\mu},y}U^\dagger_{x-\hat{\mu}, \mu}
 * \right].
 * \f]
 * We will name matrix \f$ M \f$ hopping matrix.
 *
 * \param data      Wilson's \f$ \kappa \f$
 * \param gv        gauge field configuration
 * \param out       \f$ \Psi_{out} \f$
 * \param in        \f$ \Psi_{in} \f$
 *
 * \return void
 */
/*}}}1*/

void wdirac(void* data, t_gauge_vector * gv, t_cds_vector * out, t_cds_vector * in)/*{{{*/
{
    int idx;
    t_real kappa __attribute__ ((aligned(16)));

#ifdef MU
    t_real emuplus  = exp(Mu);
    t_real emuminus = exp(-Mu);
#endif

#ifdef TIMING /*{{{*/
    start_time(TIMING_WDIRAC);
#endif /*}}}*/
#if (EXP_SSE && SINGLE_PREC) /*{{{*/
    kappa = (-1.0)*(*((t_real*)data)); /*}}}*/
#else /*{{{*/
    kappa = *((t_real*)data);
#endif /*}}}*/
    /* I * \Psi_{in} */
#if (EXP_SSE && SINGLE_PREC) /*{{{*/
    __asm__ __volatile__ (
            "movl   %2, %%ecx\n\t"
            "movl   %0, %%esi\n\t"
            "movl   %1, %%edi\n\t"
            ".LOOP_MEM_CPY_WD:\n\t"
            "movaps (%%edi), %%xmm0\n\t"
            "movaps %%xmm0, (%%esi)\n\t"
            "addl   $16, %%edi\n\t"
            "addl   $16, %%esi\n\t"
            "subl   $1, %%ecx\n\t"
            "jnz    .LOOP_MEM_CPY_WD"
            :
            : "g" (out), "g" (in), "i" (NCOLORS*NDIRAC*VOL/2)
            : "%ecx", "%esi", "%edi", "%xmm0"
            );/*}}}*/
#else /*{{{*/
    memcpy(&((*out)[0][0][0]), &((*in)[0][0][0]), sizeof(t_complex)*NCOLORS*NDIRAC*VOL);
#endif /*}}}*/
    /* hopping matrix M */
// Main Piece of Code
 t_complex col_tmp[NDIRAC*NCOLORS] __attribute__ ((aligned(16)));
 t_complex tmp[NDIRAC*NCOLORS] __attribute__ ((aligned(16)));
 int idx_moved;
 int mu;
 int i;
#ifdef _OPENMP
#pragma omp parallel for private(col_tmp, tmp, idx_moved, mu, i)
#endif
    for (idx = 0; idx < VOL; idx++) /*{{{*/
    {
      for(i = 0; i < NDIRAC*NCOLORS; i++)
       tmp[i] = 0.0;
      for(mu = 0; mu < DIM; mu++)
      {
       /*{{{*/
       idx_moved = (*lat_mov)[idx][mu][FWD];
#ifdef TIMING /*{{{*/
       start_time(TIMING_A_MUL_CDSV);
#endif/*}}}*/
       MATCALL(A_mul_CDSV)(col_tmp, (*gv)[idx][mu], &((*in)[idx_moved][0][0]));
#ifdef EFIELD
           for(i = 0; i<NDIRAC*NCOLORS; i++)
            col_tmp[i] = col_tmp[i]*ev[idx][mu];
#endif
       if(mu==0)
        for(i = 0; i<NDIRAC*NCOLORS; i++)
#ifdef MU
         col_tmp[i] *= emuplus*anisotropy_factor;
#else
         col_tmp[i] *= anisotropy_factor;
#endif
#ifdef TIMING /*{{{*/
       end_time(TIMING_A_MUL_CDSV);
#endif /*}}}*/
#if (EXP_SSE && SINGLE_PREC)
       gpr_mul_CDSV_fwd(mu, tmp, col_tmp);
#else
       gpr_mul_CDSV(mu, FWD, tmp, col_tmp);
#endif
       idx_moved = (*lat_mov)[idx][mu][BWD];
#ifdef TIMING /*{{{*/
       start_time(TIMING_AD_MUL_CDSV);
#endif /*}}}*/
       MATCALL(Ad_mul_CDSV)(col_tmp, (*gv)[idx_moved][mu], &((*in)[idx_moved][0][0]));
#ifdef EFIELD
           for(i = 0; i<NDIRAC*NCOLORS; i++)
            col_tmp[i] = col_tmp[i]*conj(ev[idx_moved][mu]);
#endif
       if(mu==0)
        for(i = 0; i<NDIRAC*NCOLORS; i++)
#ifdef MU
         col_tmp[i] *= emuminus*anisotropy_factor;
#else
         col_tmp[i] *= anisotropy_factor;
#endif
#ifdef TIMING /*{{{*/
       end_time(TIMING_AD_MUL_CDSV);
#endif /*}}}*/
#if (EXP_SSE && SINGLE_PREC)
       gpr_mul_CDSV_bwd(mu, tmp, col_tmp);
#else
       gpr_mul_CDSV(mu, BWD, tmp, col_tmp);
#endif
        }/*}}}*/
#if (EXP_SSE && SINGLE_PREC) /*{{{*/
        __asm__ __volatile__ (
                "movss  %0, %%xmm2              /* -kappa -> xmm2 */\n\t"
                "shufps $0, %%xmm2, %%xmm2\n\t"
                "mov    %1, %%ecx\n\t"
                "mov    %2, %%esi               /* &tmp -> esi */\n\t"
                "mov    %3, %%edi               /* out -> edi */\n\t"
                ".LOOP_M_KAPPA:\n\t"
                "movaps (%%edi), %%xmm0\n\t"
                "movaps (%%esi), %%xmm1\n\t"
                "mulps  %%xmm2, %%xmm1          /* tmp * (-kappa) */\n\t"
                "addps  %%xmm1, %%xmm0          /* out -= kappa*tmp */\n\t"
                "movaps %%xmm0, (%%edi)\n\t"
                "addl   $16, %%esi\n\t"
                "addl   $16, %%edi\n\t"
                "subl   $1, %%ecx\n\t"
                "jnz    .LOOP_M_KAPPA"
                :
                : "g" (kappa), "i" (NDIRAC*NCOLORS/2), "g" (&tmp),
                "g" ((t_complex *)out+idx*NDIRAC*NCOLORS)
                : "%ecx", "%esi", "%edi"
            );
        /*}}}*/
#else   /*{{{*/
        for (i = 0; i < NDIRAC*NCOLORS; i++)
        {
            /* I - \kappa M */
            *((t_complex *)out + idx*NDIRAC*NCOLORS + i)  -= kappa*tmp[i];
        };
#endif /*}}}*/
    };/*}}}*/
#ifdef TIMING /*{{{*/
    end_time(TIMING_WDIRAC);
#endif /*}}}*/
;
}/*}}}*/

/*{{{*/
/*!
 * \brief Multiply fermionic vector by Hermitian Wilson Dirac Operator
 *
 * \f[
 * H_{x,y} = \gamma_5 Q_{x,y}
 * \f]
 *
 * \param data Wilson's \f$\kappa\f$
 * \param gv
 * \param out
 * \param in
 *
 * \return void
 */ /*}}}*/
void h_wdirac(void*, t_gauge_vector*, t_cds_vector*, t_cds_vector*) __attribute__ ((noinline)); /*{{{*/
void h_wdirac(void* data, t_gauge_vector * gv, t_cds_vector * out, t_cds_vector * in)
{
#if (EXP_SSE && SINGLE_PREC)
    float tmp;
#else
    int idx, a;
#endif

#ifdef TIMING /*{{{*/
    start_time(TIMING_H_WDIRAC);
#endif /*}}}*/
    wdirac(data, gv, out, in);
#ifdef TIMING /*{{{*/
    start_time(TIMING_H_WDIRAC_LN);
#endif/*}}}*/
#if (EXP_SSE && SINGLE_PREC) /*{{{*/
    tmp = -1.0;
    __asm__ __volatile__ (
            "movss  %3, %%xmm7\n\t"
            "shufps $0, %%xmm7, %%xmm7\n\t"
            "movl   %0, %%esi\n\t"
            "movl   %2, %%ecx\n\t"
            ".LOOP_HW_G5:\n\t"
            "movaps 16(%%esi), %%xmm0\n\t"
            "movaps 48(%%esi), %%xmm1\n\t"
            "movaps 80(%%esi), %%xmm2\n\t"
            "movaps 112(%%esi), %%xmm3\n\t"
            "mulps  %%xmm7, %%xmm0\n\t"
            "mulps  %%xmm7, %%xmm1\n\t"
            "mulps  %%xmm7, %%xmm2\n\t"
            "mulps  %%xmm7, %%xmm3\n\t"
            "movaps %%xmm0, 16(%%esi)\n\t"
            "movaps %%xmm1, 48(%%esi)\n\t"
            "movaps %%xmm2, 80(%%esi)\n\t"
            "movaps %%xmm3, 112(%%esi)\n\t"
            "addl   $128, %%esi\n\t"
            "subl   $1, %%ecx\n\t"
            "jnz    .LOOP_HW_G5"
            : "=m" (out)
            : "m" (out), "i" (VOL*NCOLORS/4), "m" (tmp)
            : "%esi", "%ecx", "%xmm7", "%xmm0", "%xmm1", "%xmm2", "%xmm3"
            ); /*}}}*/
#else /*{{{*/
    for (idx = 0; idx < VOL; idx++)
        for (a = 0; a < NCOLORS; a++) {
            (*out)[idx][a][2] *= -1.0;
            (*out)[idx][a][3] *= -1.0;
        }
#endif /*}}}*/
#ifdef TIMING /*{{{*/
    end_time(TIMING_H_WDIRAC_LN);
#endif /*}}}*/
#ifdef TIMING /*{{{*/
    end_time(TIMING_H_WDIRAC);
#endif/*}}}*/
}/*}}}*/

/*{{{*/
/*!
 * \brief
 *
 * \param data Wilson's \f$\kappa\f$
 * \param gv
 * \param out
 * \param in
 * \return void
 *//*}}}*/

 /*{{{*/

void h_wdirac_sq(void*, t_gauge_vector*, t_cds_vector*, t_cds_vector*) __attribute__ ((noinline));
void h_wdirac_sq(void* data, t_gauge_vector* gv, t_cds_vector* out, /*{{{*/
        t_cds_vector* in)
{
    t_cds_vector* tmp;

#ifdef TIMING /*{{{*/
    start_time(TIMING_H_WDIRAC_SQ);
#endif /*}}}*/
    if (posix_memalign((void*)&tmp, 16, sizeof(t_cds_vector)))
        panic("%s: can't allocate memory for tmp (%i Kb)", __func__,
                sizeof(t_cds_vector)/1024);
    h_wdirac(data, gv, tmp, in);
    h_wdirac(data, gv, out, tmp);
    free(tmp);
#ifdef TIMING /*{{{*/
    end_time(TIMING_H_WDIRAC_SQ);
#endif /*}}}*/
}/*}}}*/


/*{{{*/
/*!
 * \brief  Evaluate \f$ \Psi_{out} = \tilde{M}_{ov} \Psi_{in} \f$
 *
 * We will present here some theoretical background about Neuberger's
 * Dirac operator (overlap Dirac operator) here. This is an explicit
 * construction. Introduce
 * \f[
 * A = a D_w^{m=0} - \rho = a D_w^{m=-\rho/a}
 * \f]
 * where \f$D_w^{m=-\rho/a}\f$ is a Wilson Dirac operator given
 * in the description of wdirac() but with a negative mass term
 * \f$-\frac{\rho}{a}\f$
 * \f[
 * 0 < \rho \le 2.
 * \f]
 * Let us rewrite \f$A \f$ in a following way
 * \f[
 * A_{x,y} = \frac{1}{2\kappa}\, Q_{x,y},
 * \f]
 * \f[
 * 0 < \frac{1}{8r} \le \left(\kappa = \frac{1}{-2 \rho+8r}\right)
 * < \frac{1}{8r - 4}.
 * \f]
 * Then we define
 * \f[
 * D_{ov} = \frac{\rho}{a}\left(1+\frac{A}{\sqrt{AA^\dagger}}\right) =
 * \frac{\rho}{a}\left(
 * 1+\frac{Q}{\sqrt{QQ^\dagger}}
 * \right) =
 * \frac{\rho}{a}\left(1+ \gamma_5\, sign(H)\right),\, H = \gamma_5 Q.
 * \f]
 *
 * The advantage of this definition is that the normalization of the wave
 * functio and hence renormalization constant is
 * \f[
 * Z = 1 + O(g_0^2).
 * \f]
 * This can be seen from looking at the free field case in the limit
 * \f$ p \rightarrow 0\f$. Then we have
 * \f[
 * D_w^{m=0} = i \gamma_\mu p_\mu + O(p^2),\,
 * A = i\gamma_\mu a p_\mu - \rho + O(a^2p^2),
 * \f]
 * \f[
 * D_{ov} = \frac{\rho}{a}\left(1 +
 * \frac{i \gamma_\mu a p_\mu - \rho + O(a^2p^2)}{
 * \sqrt{\rho^2 + O(a^2p^2)}}\right) = i \gamma_\mu p_\mu + O(ap^2).
 * \f]
 * This operator satisfies a slightly modified GW relation
 * \f[
 * \gamma_5 D_{ov} + D_{ov} \gamma_5 = \frac{a}{\rho}D_{ov}\gamma_5 D_{ov},
 * \f]
 * which is important for us to analyse spectrum of \f$D_{ov}\f$.
 *
 * And now let's construct massive overlap operator
 * \f[
 *  M_{ov} = \left(1 - \frac{am_q}{2 \rho}\right)D_{ov}(\rho) + m_q =
 *  \frac{1}{a}(\rho - \frac{am_q}{2})(1+\gamma_5 sign(H)) + m_q,
 * \f]
 * \f[
 * M_{ov} = \frac{\rho}{a}\left[
 * \left(1+\frac{am_q}{2\rho}\right)+(1-\frac{am_q}{2\rho})\gamma_5 sign(H)
 * \right] = \frac{\rho}{a} \tilde{M}_{ov}.
 * \f]
 * To determine normalization of the wave functions let us write
 * fermion action
 * \f[
 * S_F = a^4 \sum \bar{\Psi}M_{ov}\Psi =
 * a^3 \sum \bar{\Psi} \rho \tilde{M}_{ov}\Psi =
 * \sum \bar{\Psi}^{dimless} \tilde{M}_{ov} \Psi^{dimless},
 * \f]
 * hence \f$ \Psi^{dimless} = a^{3/2} \sqrt{\rho} \Psi\f$.
 *
 *
 * And now some technical details go:
 *
 * The key ingredient of the overlap operator is the \f$sign(H)\f$ function
 * \f[
 * sign(H) = \frac{H}{\sqrt{H^\dagger H}},
 *  \f]
 * and \f$ H \f$ is a hermitian operator, \f$ H = H^\dagger \f$. It's
 * eigenvalues are all real, i.e. \f$ spec(H) \in {\cal R}\f$.
 *
 * Out first task is to construct approximation for \f$ sign(H) \f$. What does
 * function of matrix mean? Function of matrix \f$f(A)\f$ is any polynomial
 * \f$g(A)\f$ which has exactly the same values on the spectrum of
 * matrix \f$A\f$, i.e. if \f$\Lambda_A = spec(A)\f$, then function of \f$ A \f$
 * is defined in such a way that \f$ f(\Lambda_A) = g(\Lambda_A) \f$.
 *
 * So, out approximation shoul be valid on the whole spectrum of
 * matrix \f$ H \f$ which is
 * \f[
 * spec(H) \in [\lambda_{min};\lambda_{max}] \in {\cal R}
 * \f]
 * Now we use the relation
 * \f[
 * sign(H) \equiv sign\left(\frac{H}{\|H\|}\right) = sign(W),\ \ \
 * \|H\| = \sqrt{\lambda_{max}^2} = \lambda_{max},
 * \f]
 * \f[
 * spec(W) \in \left[\frac{\lambda_{min}}{\lambda_{max}};1\right] \in {\cal R}
 * \f]
 * We are going to use minmax polynomial approximation (minmax.c) for
 * \f$\frac{1}{\sqrt{H^2}} \f$ in the range \f$ \sqrt{\epsilon} \le H \le 1 \f$. This means
 * that for the approximation which is valid on the whole spectrum of \f$ H \f$
 * we should choose
 * \f$ \epsilon = \frac{\lambda^2_{min}}{\lambda^2_{max}} \f$.
 * In this range we have
 * \f[
 * sign(H) = sign(W) \approx W P_n(W^2) =
 * \frac{H}{\|H\|}P_n\left(\frac{H}{\|H\|^2}\right) =
 * P_n\left(\frac{H}{\|H\|}\right) \frac{H}{\|H\|} = V
 * \f]
 * So we need a procedure to multiply vector \f$ \Psi_{in} \f$ by an
 * operator \f$ V \f$. Steps:
 *
 * - 1. \f$\Psi_0 = \frac{H}{\|H\|}\Psi_{in}\f$
 * - 2. \f$\Psi_{out} = P_n\left(\frac{H^2}{\|H\|^2}\right) \Psi_0\f$
 * \f[
 * \Psi_{out} =
 * \sum_{i=0}^n c_i T_i\left(z\left(\frac{H^2}{\|H\|^2}\right)\right) \Psi_0,
 * \f]
 * \f[
 * z(x) = \frac{2x-1-\epsilon}{1-\epsilon},
 * \f]
 * \f[
 * z\left(\frac{H^2}{\|H\|^2}\right) =
 * \frac{2 H^2 - (\lambda_{max}^2 + \lambda_{min}^2)}{\lambda_{max}^2 -
 * \lambda_{min}^2} = a H^2 + b
 * \f]
 * \f[
 * a =\frac{2}{\lambda_{max}^2-\lambda_{min}^2},\ \ \
 * b = -\ \frac{\lambda_{max}^2+\lambda_{min}^2}{\lambda_{max}^2-\lambda_{min}^2}
 * \f]
 * - 3. \f$\Psi_1 = T_1\left(z\left(\frac{H^2}{\|H\|^2}\right)\right)\Psi_0 =
 * (a H^2 +b)\Psi_0\f$
 * - 4. \f$\Psi_{out} = c_1 \Psi_1 + c_0 \Psi_0\f$
 * - 5. \f$\Psi_{out}\  +=
 * c_i T_i\left(z\left(\frac{H^2}{\|H\|^2}\right)\right)\Psi_0,\ \ \ i\ge 2.\f$
 * We use recurrence formula
 * \f$ T_i(x) = 2xT_{i-1}(x) - T_{i-2}(x) \f$.
 * It'll be better to use Clenshaw's Recurrence Formula.
 * \f[
 * T_i\left(z\left(\frac{H^2}{\|H\|^2}\right)\right) \Psi_0 =
 * T_i(aH^2+b)\Psi_0 = 2(aH^2+b)T_{i-1}(aH^2+b)\Psi_0 - T_{i-2}(aH^2+b)\Psi_0=
 * 2(aH^2+b)\Psi_{i-1} - \Psi_{i-2}
 * \f]
 *
 * Also we are using low-mode projectors. For details refer to hep-lat/0212012.
 *
 * \param data
 * \param gv    gauge configuration
 * \param out   \f$ \Psi_{out} \f$
 * \param in    \f$ \Psi_{in} \f$
 *
 * \return void
 *//*}}}*/
void ov_dirac_mm(void * data, t_gauge_vector* gv, t_cds_vector* out, /*{{{*/
        t_cds_vector* in)
{
    int i;
    t_ov_data_mm* ov_data;
    ov_data = (t_ov_data_mm*) data;
#ifdef MU
    fprintf(stdout, "\nSIGN > ");
    fflush(stdout);
#ifdef NOISY_OUTPUT
    fprintf(stdout, "Allocating %i Mb of memory for Lanczos vectors ... \n", (ov_data->deg + 5)*VOL*NDIRAC*NCOLORS*sizeof(t_complex)/(1024*1024) );
    fflush(stdout);
#endif

    t_complex* v = malloc(ov_data->deg*VOL*NDIRAC*NCOLORS*sizeof(t_complex));
    if(v==NULL)
     panic("%s: can't allocate memory for v (%i Kb)", __func__, ov_data->deg*VOL*NDIRAC*NCOLORS*sizeof(t_complex)/1024);
    t_complex* w = malloc(3*VOL*NDIRAC*NCOLORS*sizeof(t_complex));
    if(w==NULL)
     panic("%s: can't allocate memory for w (%i Kb)", __func__, 3*VOL*NDIRAC*NCOLORS*sizeof(t_complex)/1024);
    t_complex* rin = malloc(ov_data->n_proj*sizeof(t_complex));
    if(rin==NULL)
     panic("%s: can't allocate memory for rin (%i Kb)", __func__, ov_data->n_proj*sizeof(t_complex)/1024);
    t_complex* lin = malloc(ov_data->n_proj*sizeof(t_complex));
    if(lin==NULL)
     panic("%s: can't allocate memory for lin (%i Kb)", __func__, ov_data->n_proj*sizeof(t_complex)/1024);

#ifdef NOISY_OUTPUT
    fprintf(stdout, "DONE: Allocating %i Mb of memory for Lanczos vectors ... \n", (ov_data->deg + 5)*VOL*NDIRAC*NCOLORS*sizeof(t_complex)/(1024*1024) );
    fflush(stdout);
#endif

    /* initializing the left and right Lanczos vectors */
#ifdef NOISY_OUTPUT
    fprintf(stdout, "Initializing the Lanczos procedure ... \n");
    fflush(stdout);
#endif

    memcpy(&(v[0]), in, VOL*NDIRAC*NCOLORS*sizeof(t_complex));
    memcpy(&(w[1*VOL*NDIRAC*NCOLORS]), in, VOL*NDIRAC*NCOLORS*sizeof(t_complex));
    for(i=0; i<ov_data->n_proj; i++)
    {
     lin[i] = innprod(&(ov_data->levecs[i*VOL*NDIRAC*NCOLORS]), &(v[0]));
     rin[i] = innprod(&(ov_data->revecs[i*VOL*NDIRAC*NCOLORS]), &(w[1*VOL*NDIRAC*NCOLORS]));
     xpcby(&(v[0]), &(ov_data->revecs[i*VOL*NDIRAC*NCOLORS]), -1.0*lin[i]);
     xpcby(&(w[1*VOL*NDIRAC*NCOLORS]), &(ov_data->levecs[i*VOL*NDIRAC*NCOLORS]), -1.0*rin[i]);
    };
    t_real vvnorm = vnorm(&(v[0]));
    t_real wwnorm = vnorm(&(w[1*VOL*NDIRAC*NCOLORS]));
    ax(&(v[0]),1.0/vvnorm);
    ax(&(w[1*VOL*NDIRAC*NCOLORS]),1.0/wwnorm);

#ifdef NOISY_OUTPUT
    fprintf(stdout, "DONE: Initializing the Lanczos procedure.\n");
    fflush(stdout);
#endif

    /* Variables for the Lanczos algorithm */
#ifdef NOISY_OUTPUT
    fprintf(stdout, "Allocating auxiliary variables for the Lanczos procedure ... \n");
    fflush(stdout);
#endif

    t_complex* alpha = malloc(ov_data->deg*sizeof(t_complex));
    if(alpha==NULL)
     panic("%s: can't allocate memory for alpha (%i Kb)", __func__, ov_data->deg*sizeof(t_complex)/1024);
    t_complex* beta = malloc((ov_data->deg - 1)*sizeof(t_complex));
    if(beta==NULL)
     panic("%s: can't allocate memory for beta (%i Kb)", __func__, ov_data->deg*sizeof(t_complex)/1024);
    t_complex* gamma = malloc((ov_data->deg - 1)*sizeof(t_complex));
    if(gamma==NULL)
     panic("%s: can't allocate memory for gamma (%i Kb)", __func__, ov_data->deg*sizeof(t_complex)/1024);
    t_complex* avtmp = malloc(VOL*NDIRAC*NCOLORS*sizeof(t_complex));
    if(avtmp==NULL)
     panic("%s: can't allocate memory for avtmp (%i Kb)", __func__, VOL*NDIRAC*NCOLORS*sizeof(t_complex)/1024);
    t_complex* awtmp = malloc(VOL*NDIRAC*NCOLORS*sizeof(t_complex));
    if(awtmp==NULL)
     panic("%s: can't allocate memory for awtmp (%i Kb)", __func__, VOL*NDIRAC*NCOLORS*sizeof(t_complex)/1024);
#ifdef NOISY_OUTPUT
    fprintf(stdout, "DONE: Allocating auxiliary variables for the Lanczos procedure ... \n");
    fflush(stdout);
#endif

    /* Starting the two-sided Lanczos algorithm for asymmetric matrices */
#ifdef NOISY_OUTPUT
    fprintf(stdout, "Starting the Lanczos iterations ... \n");
    fflush(stdout);
#endif

    for(i = 0; i<ov_data->deg-1; i++)
    {
     // TODO: add reorthogonalization at each step
     t_cds_vector* arg1 = (t_cds_vector *)&(v[(i+1)*VOL*NDIRAC*NCOLORS]);
     t_cds_vector* arg2 = (t_cds_vector *)&(v[(i)*VOL*NDIRAC*NCOLORS]);
     h_wdirac(&(ov_data->kappa), gv, arg1, arg2);
     alpha[i] = innprod(&(w[1*VOL*NDIRAC*NCOLORS]), &(v[(i+1)*VOL*NDIRAC*NCOLORS]));
     Mu = -1.0*Mu; // conjugate operator
     arg1 = (t_cds_vector *)&(w[0*VOL*NDIRAC*NCOLORS]);
     arg2 = (t_cds_vector *)&(w[1*VOL*NDIRAC*NCOLORS]);
     h_wdirac(&(ov_data->kappa), gv, arg1, arg2);
     Mu = -1.0*Mu; // normal operator again
     xpcby(&(v[(i+1)*VOL*NDIRAC*NCOLORS]), &(v[i*VOL*NDIRAC*NCOLORS]), -1.0*alpha[i]);
     if(i>0)
      xpcby(&(v[(i+1)*VOL*NDIRAC*NCOLORS]), &(v[(i-1)*VOL*NDIRAC*NCOLORS]), -1.0*gamma[i-1]);
     xpcby(&(w[0*VOL*NDIRAC*NCOLORS]), &(w[1*VOL*NDIRAC*NCOLORS]), -1.0*CCALL(conj)(alpha[i]));
     if(i>0)
      xpcby(&(w[0*VOL*NDIRAC*NCOLORS]), &(w[2*VOL*NDIRAC*NCOLORS]), -1.0*CCALL(conj)(beta[i-1]));
     beta[i] = vnorm(&(v[(i+1)*VOL*NDIRAC*NCOLORS]));
     gamma[i] = innprod(&(w[0*VOL*NDIRAC*NCOLORS]), &(v[(i+1)*VOL*NDIRAC*NCOLORS]));
     gamma[i] = gamma[i]/beta[i];
     axc(&(v[(i+1)*VOL*NDIRAC*NCOLORS]), 1.0/beta[i]);
     axc(&(w[0*VOL*NDIRAC*NCOLORS]),     1.0/CCALL(conj)(gamma[i]));
     memcpy(&(w[2*VOL*NDIRAC*NCOLORS]), &(w[1*VOL*NDIRAC*NCOLORS]), VOL*NDIRAC*NCOLORS*sizeof(t_complex));
     memcpy(&(w[1*VOL*NDIRAC*NCOLORS]), &(w[0*VOL*NDIRAC*NCOLORS]), VOL*NDIRAC*NCOLORS*sizeof(t_complex));
     fprintf(stdout, ".");
     fflush(stdout);
    };

    t_cds_vector* arg1 = (t_cds_vector *)&(w[0*VOL*NDIRAC*NCOLORS]);
    t_cds_vector* arg2 = (t_cds_vector *)&(v[(ov_data->deg - 1)*VOL*NDIRAC*NCOLORS]);
    h_wdirac(&(ov_data->kappa), gv, arg1, arg2);
    alpha[(ov_data->deg - 1)] = innprod(&(w[1*VOL*NDIRAC*NCOLORS]), &(w[0*VOL*NDIRAC*NCOLORS]));

#ifdef NOISY_OUTPUT
    fprintf(stdout, "\n Lanczos iterations DONE. \n");
    fflush(stdout);
#endif
    /* Now we have the coefficients of the tridiagonal Lanczos matrix and the Lanczos vectors */
    /* Now the sign function of the Lanczos matrix should be calculated */
    /* The sign function is the matrix sgn */
#ifdef NOISY_OUTPUT
    fprintf(stdout, "Calculating the sign function of a tridiagonal matrix ... \n");
    fflush(stdout);
#endif
    t_double_complex* sgn = malloc(ov_data->deg*ov_data->deg*sizeof(t_double_complex));
    if(sgn==NULL)
     panic("%s: can't allocate memory for sgn (%i Kb)", __func__, ov_data->deg*ov_data->deg*sizeof(t_double_complex)/1024);
    // Initializing the tridiagonal matrix used in Robert's method
    for(i=0; i<ov_data->deg*ov_data->deg; i++)
     sgn[i] = 0.0 + I*0.0;
    sgn[0] = (t_double_complex)alpha[0];
    sgn[1] = (t_double_complex)gamma[0];
    for(i=1; i<ov_data->deg-1; i++)
    {
     sgn[i*ov_data->deg + i] = (t_double_complex)alpha[i];
     sgn[i*ov_data->deg + i - 1] = (t_double_complex)beta[i-1];
     sgn[i*ov_data->deg + i + 1] = (t_double_complex)gamma[i];
    };
    sgn[ov_data->deg*ov_data->deg - 1] = (t_double_complex)alpha[ov_data->deg - 1];
    sgn[ov_data->deg*ov_data->deg - 2] = (t_double_complex)beta[ov_data->deg - 2];

    t_double_complex* nsgn = malloc(ov_data->deg*ov_data->deg*sizeof(t_double_complex));
    if(nsgn==NULL)
     panic("%s: can't allocate memory for nsgn (%i Kb)", __func__, ov_data->deg*ov_data->deg*sizeof(t_double_complex)/1024);
    // Starting the Robert's iterations
#ifdef NOISY_OUTPUT
    fprintf(stdout, "Starting the Robert's iterations S_n+1 = 0.5(S_n + S_n^-1)... \n");
    fflush(stdout);
#else
    fprintf(stdout, "|>");
#endif
    t_double_real err = 0;
    int count = 0;
    while((err > ov_data->sign_tol) || (count == 0))
    {
     inverse(sgn, nsgn, ov_data->deg);
     axpbynd(sgn, 0.5, nsgn, 0.5, ov_data->deg*ov_data->deg);
     err = adiffnd(sgn, nsgn, ov_data->deg*ov_data->deg);
     count++;
#ifdef NOISY_OUTPUT
     fprintf(stdout,"Err: %8.8Lf, sign_tol: %8.8Lf, ov_data->deg: %i\n", err, ov_data->sign_tol, ov_data->deg);
     fflush(stdout);
#else
     fprintf(stdout,".");
     fflush(stdout);
#endif
    };

#ifdef NOISY_OUTPUT
    fprintf(stdout, "\n\n DONE: Calculating the sign function of a tridiagonal matrix.\n");
    fflush(stdout);
#endif
    /* Finally calculating the sign */

    /* Critical eigenvalues */
#ifdef NOISY_OUTPUT
    fprintf(stdout, "Deflating the critical eigenvalues ... \n");
    fflush(stdout);
#endif

    set_zero((t_complex *)&(out[0]));
    for(i = 0; i<ov_data->n_proj; i++)
    {
     t_real f = CCALL(creal)(ov_data->evals[i]) >= 0 ? 1.0 : -1.0;
     xpcby((t_complex *)&(out[0]), &(ov_data->revecs[i*VOL*NDIRAC*NCOLORS]), f*lin[i]);
    };

#ifdef NOISY_OUTPUT
    fprintf(stdout, "DONE: Deflating the critical eigenvalues. \n");
    fflush(stdout);
#endif
    /* Krylov subspace */
#ifdef NOISY_OUTPUT
    fprintf(stdout, "Deflating the Krylov subspace ... \n");
    fflush(stdout);
#endif

    set_zero(&(w[0*VOL*NDIRAC*NCOLORS])); // This will temporary store the Krylov part
    for(i = 0; i<ov_data->deg; i++)
    {
     xpcby(&(w[0*VOL*NDIRAC*NCOLORS]), &(v[i*VOL*NDIRAC*NCOLORS]), (t_complex)sgn[ov_data->deg*i]);
    };
    /* re-projecting the Lanczos vectors */
#ifdef NOISY_OUTPUT
    fprintf(stdout, "Re-projecting the Lanczos vectors ... \n");
    fflush(stdout);
#endif

    for(i = 0; i<ov_data->n_proj; i++)
    {
     t_complex vl = innprod(&(ov_data->levecs[i*VOL*NDIRAC*NCOLORS]),&(w[0*VOL*NDIRAC*NCOLORS]));
     xpcby(&(w[0*VOL*NDIRAC*NCOLORS]), &(ov_data->revecs[i*VOL*NDIRAC*NCOLORS]), -1.0*vl);
    };
    xpby((t_complex *)&(out[0]), &(w[0*VOL*NDIRAC*NCOLORS]), vvnorm);

#ifdef NOISY_OUTPUT
    fprintf(stdout, "DONE: Deflating the Krylov subspace. \n");
    fflush(stdout);
#endif
    /* Freeing up the memory */
    free(sgn);
    free(nsgn);
    free(awtmp);
    free(avtmp);
    free(alpha);
    free(beta);
    free(gamma);
    free(v);
    free(w);
    free(lin);
    free(rin);
#else
    t_real a, b;
    t_cds_vector* psi[3];
    t_cds_vector* psi_tmp;
    t_real norm;
#ifdef OV_DIRAC_CONV_LOG /*{{{*/
    t_complex tnorm;
    FILE* tfile;
#endif /*}}}*/
#ifdef OV_DIRAC_ORTHO_PROJ /*{{{*/
    t_complex* dotpr;
#endif /*}}}*/
#ifdef OV_DIRAC_PROG_BAR /*{{{*/
    int jpb;
    double oneperc;
#endif/*}}}*/

#ifdef TIMING /*{{{*/
    start_time(TIMING_OVDIRAC);
#endif /*}}}*/
    for (i = 0; i < 3; i++)
        if (posix_memalign((void*)&(psi[i]), 16, sizeof(t_cds_vector)))
            panic("%s: can't allocate memory for psi[%i] (%i Kb)", i, __func__, sizeof(t_cds_vector)/1024);
    a = 2.0/(ov_data->l_sq_max - ov_data->l_sq_min);
    b = (ov_data->l_sq_max + ov_data->l_sq_min)/(ov_data->l_sq_min - ov_data->l_sq_max);
    /* \Psi_0 */
    h_wdirac(&(ov_data->kappa), gv, psi[0], in);
#ifdef OV_DIRAC_ORTHO_PROJ /*{{{*/
    if ((dotpr = malloc(sizeof(t_complex)*(ov_data->n_proj))) == NULL)
        panic("%s: can't allocate memory for dotpr(n_proj=%i) (%i Kb)", __func__,
                ov_data->n_proj,
                sizeof(t_complex)*(ov_data->n_proj)/1024);
    proj_ortho_comp((t_complex*)(psi[0]), ov_data->evecs, ov_data->n_proj, dotpr);
    for (i = 0; i < ov_data->n_proj; i++)
         dotpr[i] = innprod((ov_data->evecs + i*VOL*NDIRAC*NCOLORS), (t_complex*)(in));
#endif /*}}}*/
    /* normalize */
    norm = 1.0/sqrt(fabs(ov_data->l_sq_max));
    for (i = 0; i < VOL*NDIRAC*NCOLORS; i++)
        *((t_complex*)(psi[0]) + i) *= norm;
    /* \Psi_1 */
    h_wdirac_sq(&(ov_data->kappa), gv, psi[1], psi[0]);
    axpby((t_complex*)(psi[1]), a, (t_complex*)(psi[0]), b);
    /* \Psi_{out} = c_1 \Psi_1 + c_0 \Psi_0 */
    xcpy((t_complex*)out, (t_complex*)(psi[0]));
    axpby((t_complex*)out, (t_real)(ov_data->coef[0]), (t_complex*)(psi[1]), (t_real)(ov_data->coef[1]));
#ifdef OV_DIRAC_CONV_LOG /*{{{*/
    if ((tfile = fopen(OV_DIRAC_CONV_FNAME, "wb")) == NULL)
        panic("%s: can't open file %s for writing", __func__, OV_DIRAC_CONV_FNAME);
    fprintf(tfile, "# itnum norm(out) norm(psi_i) coef_i\n");
    tnorm = innprod((t_complex*)(out), (t_complex*)(out));
    fprintf(tfile, "%5i\t%12.10f\t", 1, creal(tnorm));
    tnorm = innprod((t_complex*)(psi[1]), (t_complex*)(psi[1]));
    fprintf(tfile, "%12.10f\t%12.10f\n", creal(tnorm), ov_data->coef[1]);
#endif /*}}}*/
    /* \Psi_{out} += ... */
#ifdef OV_DIRAC_PROG_BAR /*{{{*/
    oneperc = 100.0/ov_data->deg;
#else /*}}}*/
 fprintf(stdout, "OVERLAP ...");
 fflush(stdout);
#endif
    for (i = 2; i <= ov_data->deg; i++) {
#ifdef OV_DIRAC_PROG_BAR /*{{{*/
        fprintf(stdout, "\rOV: ");
        for (jpb = 0; jpb < 70*oneperc*i/100; jpb++)
            fprintf(stdout, ".");
        fprintf(stdout, "%3.1f%%", oneperc*i);
        fflush(stdout);
#endif /*}}}*/
        h_wdirac_sq(&(ov_data->kappa), gv, psi[2], psi[1]);
#ifdef TIMING /*{{{*/
        start_time(TIMING_OV_DIRAC_LN);
#endif /*}}}*/
        axpby((t_complex*)(psi[2]),  a, (t_complex*)(psi[1]), b);
        axpby((t_complex*)(psi[2]), 2.0, (t_complex*)(psi[0]), -1.0);
        xpby((t_complex*)out, (t_complex*)(psi[2]), (t_real)(ov_data->coef[i]));
#ifdef OV_DIRAC_CONV_LOG /*{{{*/
        tnorm = innprod((t_complex*)(out), (t_complex*)(out));
        fprintf(tfile, "%5i\t%12.10f\t", i, creal(tnorm));
        tnorm = innprod((t_complex*)(psi[2]),(t_complex*)(psi[2]));
        fprintf(tfile, "%12.10f\t%12.10f\n", creal(tnorm), ov_data->coef[i]);
#endif /*}}}*/
#ifdef TIMING /*{{{*/
        end_time(TIMING_OV_DIRAC_LN);
#endif /*}}}*/
        /* psi[2] -> psi[1] -> psi[0] -> psi[2](new) */
        psi_tmp = psi[0];
        psi[0] = psi[1];
        psi[1] = psi[2];
        psi[2] = psi_tmp;
    }
#ifdef OV_DIRAC_ORTHO_PROJ /*{{{*/
    /* if you are calculating free or degenerate spectrum -- comment line below*/
    /* clear explanation is still missing */
    proj_ortho_comp((t_complex*)out, ov_data->evecs, ov_data->n_proj, NULL);
    for (i = 0; i < ov_data->n_proj; i++) {
        t_real mult;
        mult = CCALL(creal)(ov_data->evals[i]) > 0.0 ? 1.0 : -1.0;
        xpcby((t_complex*)out, (t_complex*)(ov_data->evecs+i*VOL*NDIRAC*NCOLORS), mult*dotpr[i]);
//      mult = 1.0/fabs(CCALL(creal)(ov_data->evals[i]));
//      xpcby((t_complex*)out, (t_complex*)(ov_data->evecs+i*VOL*NDIRAC*NCOLORS), mult*dotpr[i]);
    }
    free(dotpr);
#endif /*}}}*/
    for(i = 0; i < 3; i++)
     free(psi[i]);
#ifdef OV_DIRAC_CONV_LOG /*{{{*/
    fclose(tfile);
#endif /*}}}*/
#ifdef OV_DIRAC_PROG_BAR /*{{{*/
        fprintf(stdout, "\n");
#else /*}}}*/
        fprintf(stdout, " FINISHED\n");
        fflush(stdout);
#endif
#ifdef TIMING /*{{{*/
    end_time(TIMING_OVDIRAC);
#endif/*}}}*/
#endif
    ag5xpby((t_complex*)out, (1.0-ov_data->mass), (t_complex*)in, (1.0+ov_data->mass));
}/*}}}*/

void ov_dirac_unitary(void *data, t_gauge_vector* gv, t_cds_vector* out, t_cds_vector* in)
{
 t_ov_data_mm* ov_data;
 ov_data = (t_ov_data_mm*) data;
 ov_dirac_mm(data, gv, out, in);
 axpby((t_complex*)out, 1.0/(1.0 - ov_data->mass), (t_complex*)in, -(1.0 + ov_data->mass)/(1.0 - ov_data->mass));
}

/*{{{*/
/*!
 * \brief Evaluete \f$ \Psi_{out} = \tilde{M_{ov}} \Psi_{in} \f$.
 *
 * The same as ov_dirac_mm() but Clenshaw's recurrence formula is used.
 * For more details refer to Numerical Recipes.
 *
 * \param data
 * \param gv
 * \param out
 * \param in
 *
 * \return void
 *//*}}}*/
void ov_dirac_mm_cc(void * data, t_gauge_vector* gv, t_cds_vector* out, /*{{{*/
        t_cds_vector* in)
{
#ifdef MU

#else
    int i;
    t_ov_data_mm* ov_data;
    t_real a, b;
    t_cds_vector* psi[3];
    t_cds_vector* psi_tmp;
    t_real norm;
#ifdef OV_DIRAC_ORTHO_PROJ /*{{{*/
    t_complex* dotpr;
#endif /*}}}*/

#ifdef TIMING
    start_time(TIMING_OVDIRAC);
#endif
    for (i = 0; i < 3; i++)
        if (posix_memalign((void*)&(psi[i]), 16, sizeof(t_cds_vector)))
            panic("%s: can't allocate memory for psi[%i] (% Kb)", i, __func__,
                    sizeof(t_cds_vector)/1024);
    ov_data = (t_ov_data_mm*) data;
    a = 2.0/(ov_data->l_sq_max - ov_data->l_sq_min);
    b = (ov_data->l_sq_max + ov_data->l_sq_min)/(ov_data->l_sq_min - ov_data->l_sq_max);
    /* \Psi_0 */
    h_wdirac(&(ov_data->kappa), gv, out, in);
#ifdef OV_DIRAC_ORTHO_PROJ /*{{{*/
    if ((dotpr = malloc(sizeof(t_complex)*(ov_data->n_proj))) == NULL)
        panic("%s: can't allocate memory for dot products dotpr (% Kb)", __func__,
                sizeof(t_complex)*(ov_data->n_proj)/1024);
    proj_ortho_comp((t_complex*)(out), ov_data->evecs, ov_data->n_proj, dotpr);
#endif /*}}}*/
    /* normalize */
    norm = 1.0/sqrt(ov_data->l_sq_max);
    for (i = 0; i < VOL*NDIRAC*NCOLORS; i++) {
        *((t_complex*)out + i) *= norm;
        *((t_complex*)psi[1] + i) = 0.0;
        *((t_complex*)psi[2] + i) = 0.0;
    }
    /* \Psi_{out} += ... */
    for (i = ov_data->deg; i >= 1; i--)
    {
        /* z(H^2/|H|^2) */
        h_wdirac_sq(&(ov_data->kappa), gv, psi[0], psi[1]);
        axpby((t_complex*)psi[0],  a, (t_complex*)psi[1], b);
        /* */
        axpby((t_complex*)psi[0], 2.0, (t_complex*)psi[2], -1.0);
        xpby((t_complex*)(psi[0]), (t_complex*)out, (t_real) ov_data->coef[i]);
        psi_tmp = psi[2];
        psi[2] = psi[1];
        psi[1] = psi[0];
        psi[0] = psi_tmp;
    }
    /* d0 */
    h_wdirac_sq(&(ov_data->kappa), gv, psi[0], psi[1]);
    axpby((t_complex*)psi[0],  a, (t_complex*)psi[1], b);
    xpby((t_complex*)psi[0], (t_complex*)psi[2], (t_real)(-1.0));
    axpby((t_complex*)out, (t_real)(ov_data->coef[0]), (t_complex*)psi[0], (t_real)(1.0));
    /* */
#ifdef OV_DIRAC_ORTHO_PROJ /*{{{*/
    proj_ortho_comp((t_complex*)out, ov_data->evecs, ov_data->n_proj, NULL);
    for (i = 0; i < ov_data->n_proj; i++) {
        t_real mult;
        mult = 1.0/fabs(CCALL(creal)(ov_data->evals[i]));
        xpcby((t_complex*)out, (t_complex*)(ov_data->evecs+i*VOL*NDIRAC*NCOLORS),
                mult*dotpr[i]);
    }
    free(dotpr);
#endif /*}}}*/
    ag5xpby((t_complex*)out, (1.0-ov_data->mass), (t_complex*)in, (1.0+ov_data->mass));
    for (i = 0; i < 3; i++)
        free(psi[i]);
#ifdef TIMING /*{{{*/
    end_time(TIMING_OVDIRAC);
#endif/*}}}*/
#endif
}/*}}}*/

/*{{{1*/
/*!
 * \brief Calculate eigenvalues of Dirac Operator specified in evd structure
 *
 * TODO Eigenvector should also be stored
 *
 * \param gv        gauge configuration
 * \param evd       t_ov_ev
 *
 * \return void
 */
/*}}}1*/
void arnoldi(t_gauge_vector * gv, t_op_ev* evd, int ncv_factor_nom, int ncv_factor_denom) /*{{{*/
{
    int i;
    int ido;                    /* reverse communication interface */
    char bmat = 'I';            /* standart eigenvalue problem */
    int n = VOL*NDIRAC*NCOLORS; /* Dimension of the eigenproblem */
    char which[2];
    /* number of eigenvales to be computed. 0 < nev < n-1*/
    int nev = -1;
    t_real tol;
    int ncv;                    /* 2 <= ncv-nev and ncv <= n */
    t_complex * resid = NULL;
    t_complex * v = NULL;
    t_complex * workd = NULL;
    t_complex * workl = NULL;
    t_real * rwork = NULL;
    int ldv = VOL*NDIRAC*NCOLORS;
    int iparam[11];
    int ipntr[14];
    int lworkl;
    int info = 0;
    int cont, count;
    /* neupd variables */
    int rvec; /* 0 - eigenvalues, 1 - eigenvectors+eigenvalues */
    int *select, ldz;
    char howmny;
    t_complex *d, *z, sigma, *workev;

    if (evd == NULL)
        panic("%s: evd == NULL", __func__);
    if (evd->which != NULL)
        strncpy(which, evd->which, 2*sizeof(char));
    else
        panic("%s: evd->which == NULL", __func__);
    if (evd->nev <= 0 || evd->nev >= n-1)
        panic("%s: nev(=%i) should be 0 < nev(=%i) < n(=%i)-1", __func__,
                evd->nev, evd->nev, n);
    else
        nev = evd->nev;

    ncv = ncv_factor_nom*nev/ncv_factor_denom;

    if (ncv-nev < 2)
        ncv = 2+nev;
    if (ncv > n)
        ncv = n;
    tol = (t_real) evd->tol;
    /* naupd arguments */
    if ((resid = malloc(sizeof(t_complex)*VOL*NDIRAC*NCOLORS)) == NULL)
        panic("%s: can't allocate memory for resid (%i Kb)", __func__,
                sizeof(t_complex)*VOL*NDIRAC*NCOLORS/1024);
    if (posix_memalign((void*)&v, 16, sizeof(t_complex)*VOL*NDIRAC*NCOLORS*ncv))
        panic("%s: can't allocate memory for v (%i Kb)", __func__,
                sizeof(t_complex)*VOL*NDIRAC*NCOLORS*ncv/1024);
    if (posix_memalign((void*)&workd, 16, sizeof(t_complex)*VOL*NDIRAC*NCOLORS*3))
        panic("%s: can't allocate memory for workd (%i Kb)", __func__,
                sizeof(t_complex)*VOL*NDIRAC*NCOLORS*3/1024);
    lworkl = 3*ncv*ncv+5*ncv;
    if ((workl = malloc(sizeof(t_complex)*lworkl))== NULL)
        panic("%s: can't allocate memory for workl (%i Kb)", __func__,
                sizeof(t_complex)*lworkl/1024);
    if ((rwork = malloc(sizeof(double)*ncv)) == NULL)
        panic("%s: can't allocate memory for rwork (%i Kb)", __func__,
                sizeof(double)*ncv/1024);
    for (i = 0; i < 11; i++)
        iparam[i] = 0;
    iparam[0] = 1;
    iparam[2] = MAX_ARNOLDI_ITER;
    iparam[3] = 1;
    iparam[6] = 1;
    /* ARPACK arnoldi routine */
    cont = 1;
    count = 0;
    ido = 0; /* first call to reverse communication interface */
#ifdef ARNOLDI_LOG
    UNDERSC2(start_log_arnoldi)();
#endif
    while(cont == 1)
    {
        count++;
#if (TIMING && DEBUG)
        fprintf(stderr, "count = %i\n", count);
#endif
        naupd(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv,
                iparam, ipntr, workd, workl, &lworkl, rwork, &info, 1, 2);
        if (abs(ido) == 1)
            (*(evd->op))(evd->data, gv, (t_cds_vector*)(workd+ipntr[1]-1),
                    (t_cds_vector*)(workd + ipntr[0]-1));
        else
            cont = 0;
    }
    /* get eigenvalues */
    if (info != 0)
        warn("%s: some errors occure during naupd (info = %i)", __func__, info);
    /* neupd arguments */
    if (evd->val_vec == 0)
        rvec = 0;
    else if (evd->val_vec == 1)
        rvec = 1;
    else
        panic("%s: val_vec(=%i) should be 0 or 1", __func__, evd->val_vec);
    howmny = 'A';
    if ((select = malloc(sizeof(int)*ncv)) == NULL)
        panic("%s: can't allocate memory for select (%i)", __func__,
                sizeof(int)*ncv);
    if (evd->evals != NULL)
        panic("%s: evd->evals != NULL", __func__);
    else
        if ((d = malloc(sizeof(t_complex)*(nev+1))) == NULL)
            panic("%s: can't allocate memory for d (%i Kb)", __func__,
                    sizeof(t_complex)*(nev+1));
    /* we can use already allocated array v to store eigenvectors if
     * we are not interested in Schur basis. Read neupd.f for details */
    if (evd->evecs != NULL)
        panic("%s: evd->evecs != NULL", __func__);
    else
        z = v;
    ldz = n;
    sigma = 0.0 + I*0.0;
    if ((workev = malloc(sizeof(t_complex)*2*ncv)) == NULL)
        panic("%s: can't allocate memory for workev (%i)",__func__,
                sizeof(t_complex)*2*ncv);
    neupd(&rvec, &howmny, select, d, z, &ldz, &sigma, workev, &bmat, &n, which,
            &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl,
            &lworkl, rwork, &info, 1, 2, 1);
    if (info != 0)
        warn("%s: info(%i) != 0", __func__, info);
    evd->evals = d;
    /* eigenvectors */
    evd->evecs = z;
    /* memory cleanup */
    free(resid);
    free(workd);
    free(workl);
    free(rwork);
    free(select);
    free(workev);
#ifdef ARNOLDI_LOG
    UNDERSC2(end_log_arnoldi)();
#endif
}/*}}}*/

#define IDSTR_LEN   6
#define SEC_LEN     7
#define SEC_IH      0
#define SEC_VAL     1
#define SEC_VEC     2
#define SEC_NEV     3
#define SEC_TS      4
#define SEC_OVAL    5
#define SEC_OVEC    6

/*{{{*/
/*!
 * \brief Read eigenvalues/eigenvectors from file
 *
 * Format:
 *
 * [HEADER]
 *
 * - <IDSTR:"EVDATA">
 * - <CFGNAME:char[256]>
 * - <DIM:int>
 * - <LS:int>
 * - <LT:in>
 * - <NCOLORS:int>
 * - <NDIRAC:int>
 *
 * [SECTION_SM]
 *
 * - <ISHERE:int>
 * - <ISVAL:int>
 * - <ISVEC:int>
 * - <NEV:int>
 * - <TYPESIZE:int>
 * - <OFFSET_EVALS:int>
 * - <OFFSET_EVECS:int>
 *
 * [SECTION_LM] = [SECTION_SM]
 *
 * [DATA]
 *
 * - [PSI_0]...[PSI_j]...[PSI_N]
 *
 * j = (x_0 + LS*x_1 + LS^2*x_2 + LS^3*x_3)*NDIRAC*NCOLORS + a,
 * a \in [0, NDIRAC*NCOLORS]
 * [PSI_j] = Re(PSI_j),Im(PSI_j),
 * sizeof(Re(PSI_j))=sizeof(Im(PSI_j)) = TYPESIZE
 * N = LS^3*LT*NDIRAC*NCOLORS*2
 *
 * TODO Better error recognition stuff (mainly for read_ev program)
 *
 * \param fname read from file fname
 * \param evd   put eigenvaules/eigenvectors here. evd->evals and evd->evecs
 * should be NULLs.
 * \param sm_lm read <SECTION_SM> or <SECTION_LM> of file
 */
/*}}}*/
int read_ev(char * fname, t_op_ev * evd, int sm_lm) /*{{{*/
{
    char IDSTR[IDSTR_LEN+1] = "EVDATA\0";
    FILE* dat;
    char tmp[256] = "";
    int slen, itmp, i, j;
    int my_sec_data[SEC_LEN] = {0, 0, 0, 0, 0, 0, 0};
    int sec_off;
    int dim;
    int fullvol = 1;
    int result = E_READ_NOERROR;

    if ((dat = fopen(fname, "rb")) == NULL)
        return E_READ_CANTOPEN;
    slen = fread(tmp, sizeof(char), IDSTR_LEN, dat);
    if (slen == 0)
        return E_READ_CANTREAD;
    if (strcmp(tmp, IDSTR))
        return E_READ_UNSUPFORM;
    slen = fread(tmp, sizeof(char), 256, dat);
    if (slen < 256)
        return E_READ_CANTREAD;
    /* CHECKSUM TEST SHOULD BE HERE */
    if ((slen = fread(&dim, sizeof(int), 1, dat)) != 1) {
        warn("%s: can't read <DIM> from %s", __func__, fname);
        return E_READ_CANTREAD;
    }
    if (dim != DIM) {
        warn("%s: <DIM>(=%i) != DIM(=%i)", __func__, itmp, DIM);
        result |= E_READ_BADDIM;
    }
    evd->dim = dim;
    /* LS */
    if ((slen = fread(&itmp, sizeof(int), 1, dat)) != 1) {
        warn("%s: can't read <LS> from %s", __func__, fname);
        return E_READ_CANTREAD;
    }
    if (itmp != LS) {
        warn("%s: <LS>(=%i) != LS(=%i)", __func__, itmp, LS);
        result |= E_READ_BADLS;
    }
    evd->ls = itmp;
    for (i = 0; i < dim-1; i++)
        fullvol *= itmp;
    /* LT */
    if ((slen = fread(&itmp, sizeof(int), 1, dat)) != 1) {
        warn("%s: can't read <LT> from %s", __func__, fname);
        return E_READ_CANTREAD;
    }
    if (itmp != LT) {
        warn("%s: <LT>(=%i) != LT(=%i)", __func__, itmp, LT);
        result |= E_READ_BADLT;
    }
    evd->lt = itmp;
    fullvol *= itmp;
    /* NCOLORS */
    if ((slen = fread(&itmp, sizeof(int), 1, dat)) != 1) {
        warn("%s: can't read <NCOLORS> from %s", __func__, fname);
        return E_READ_CANTREAD;
    }
    if (itmp != NCOLORS) {
        warn("%s: <NCOLORS>(=%i) != NCOLORS(=%i)", __func__, itmp, NCOLORS);
        result |= E_READ_BADNCOLORS;
    }
    evd->ncolors = itmp;
    fullvol *= itmp;
    /* NDIRAC */
    if ((slen = fread(&itmp, sizeof(int), 1, dat)) != 1) {
        warn("%s: can't read <NDIRAC> from %s", __func__, fname);
        return E_READ_CANTREAD;
    }
    if (itmp != NDIRAC) {
        warn("%s: <NDIRAC>(=%i) != NDIRAC(=%i)", __func__, itmp, NDIRAC);
        result |= E_READ_BADNDIRAC;
    }
    evd->ndirac = itmp;
    fullvol *= itmp;
    if (sm_lm == 0) {
        sec_off = (IDSTR_LEN+256)*sizeof(char)+5*sizeof(int);
    } else if (sm_lm == 1) {
        sec_off = (IDSTR_LEN+256)*sizeof(char)+5*sizeof(int)+7*sizeof(int);
    } else
        panic("%s: unsupported section(=%i)", __func__, sm_lm);
    fseek(dat, sec_off, SEEK_SET);
    /* READ SM/LM SECTION */
    if ((slen = fread(&itmp, sizeof(int), 1, dat)) != 1)
        return E_READ_CANTREAD;
    if (itmp != 1)
        return E_READ_SECNOTHERE;
    if ((slen = fread(&(my_sec_data[1]), sizeof(int), SEC_LEN-1, dat)) != SEC_LEN-1)
            panic("%s: can't read %s section", __func__, (sm_lm == 1 ? "LM" : "SM"));
    evd->nev = my_sec_data[SEC_NEV];
    evd->data_type = my_sec_data[SEC_TS];
    if (my_sec_data[SEC_VAL]) { /*{{{*/
        void * rtmp;
        if (evd->evals != NULL)
            panic("%s: evd->evals should be NULL", __func__);
        if ((evd->evals = (t_complex*)malloc(sizeof(t_complex)*evd->nev)) == NULL)
            panic("%s: can't allocate memory for evd->evals", __func__);
        fseek(dat, my_sec_data[SEC_OVAL], SEEK_SET);
        rtmp = (void*)malloc(sizeof(char)*my_sec_data[SEC_TS]*2);
        for (i = 0; i < my_sec_data[SEC_NEV]; i++) {
            fread(rtmp, my_sec_data[SEC_TS], 2, dat);
            if (my_sec_data[SEC_TS] == 4) {
                ((t_complex*)(evd->evals))[i] = ((float*)rtmp)[0] + I*((float*)rtmp)[1];
            } else if (my_sec_data[SEC_TS] == 8) {
                ((t_complex*)(evd->evals))[i] = ((double*)rtmp)[0] + I*((double*)rtmp)[1];
            }
        }
    }/*}}}*/
    if (my_sec_data[SEC_VEC]) { /*{{{*/
        void * rtmp;
        if (evd->evecs != NULL)
            panic("%s: evd->evecs should be NULL", __func__);
        if ((evd->evecs = (t_complex*)malloc(sizeof(t_complex)*evd->nev*fullvol)) == NULL)
            panic("%s:can't allocate memory for evd->evecs", __func__);
        fseek(dat, my_sec_data[SEC_OVEC], SEEK_SET);
        rtmp = (void*) malloc(sizeof(char)*my_sec_data[SEC_TS]*2);
        for (i = 0; i < my_sec_data[SEC_NEV]; i++) {
#ifdef CHECK_NORM
            t_complex norm = 0;
#endif
            for (j = 0; j < fullvol; j++) {
                slen = fread(rtmp, my_sec_data[SEC_TS], 2, dat);
                if (my_sec_data[SEC_TS] == 4) {
                    *(evd->evecs+i*fullvol+j) = ((float*)rtmp)[0] + I*((float*)rtmp)[1];
                } else if (my_sec_data[SEC_TS] == 8) {
                    *(evd->evecs+i*fullvol+j) = ((double*)rtmp)[0] + I*((double*)rtmp)[1];
                }
            }
#ifdef CHECK_NORM
            norm = innprod(evd->evecs+i*fullvol, evd->evecs+i*fullvol);
            fprintf(stderr, "norm = %f + I %f\n", creal(norm), cimag(norm));
#endif
        }
    }/*}}}*/
    fclose(dat);

    return result;
}/*}}}*/

/*{{{*/
/*!
 * \brief  Write eigenvalues/eigenvectorss to file
 *
 * \param fname save to this file
 * \param cname name of configuration file
 * \param evd   eigenvalues/eigenvectors
 * \param sm_lm write (sm_lm == 0 ? <SECTION_SM> : <SECTION_LM>
 */
/*}}}*/
void write_ev(char * fname, char * cname, t_op_ev * evd, int sm_lm) /*{{{*/
{
    FILE* dat;
    char IDSTR[IDSTR_LEN+1] = "EVDATA\0";
    char tmp[256] = "";
    int i, j, slen, itmp;
    int head_len;
    int sec_len;
    int data_off;
    int my_sec_data[SEC_LEN] = {0, 0, 0, 0, 0, 0, 0};
    int other_sec_data[SEC_LEN] = {0, 0, 0, 0, 0, 0, 0};

    head_len = (IDSTR_LEN+256)*sizeof(char)+5*sizeof(int);
    sec_len = (7*sizeof(int));
    data_off = head_len+2*sec_len;
    if ((dat = fopen(fname, "r+b")) == NULL) {
        if ((dat = fopen(fname, "w+b")) == NULL)
            panic("%s: can't open file %s for reading and writing", __func__, fname);
    }
    rewind(dat);
    slen = fread(tmp, sizeof(char), IDSTR_LEN, dat);
    if (slen == 0 || strcmp(tmp, IDSTR))
    {
        /* write header */
        rewind(dat);
        /* ID string */
        fwrite(IDSTR, sizeof(char)*IDSTR_LEN, 1, dat);
        /* <CFNAME> */
        slen = strlen(cname);
        fwrite(cname, sizeof(char), slen, dat);
        itmp = 0;
        for (i = 0; i < 256-slen; i++)
            fwrite(&itmp, sizeof(char), 1, dat);
        /* <DIM> */
        itmp = DIM;
        fwrite(&itmp, sizeof(int), 1, dat);
        /* <LS> */
        itmp = LS;
        fwrite(&itmp, sizeof(int), 1, dat);
        /* <LT> */
        itmp = LT;
        fwrite(&itmp, sizeof(int), 1, dat);
        /* <NCOLORS> */
        itmp = NCOLORS;
        fwrite(&itmp, sizeof(int), 1, dat);
        /* <NDIRAC> */
        itmp = NDIRAC;
        fwrite(&itmp, sizeof(int), 1, dat);
    }
    /* TODO check correct header (ls, lt, dim, ncolors, ndirac) */

    /* SM and LM sections */
    fseek(dat, head_len, SEEK_SET);
    if (sm_lm == LM_SECTION)
    {
     slen = fread(&itmp, sizeof(int), 1, dat);
     if(slen != 1 || itmp != 1)
     {
      /* create empty SM section */
      fseek(dat, head_len, SEEK_SET);
      itmp = 0;
      for(i = 0; i < sec_len; i++)
       fwrite(&itmp, sizeof(char), 1, dat);
     }
     else
     {
      /* SM section is already here */
      other_sec_data[SEC_IH] = itmp;
      fread(&(other_sec_data[1]), sizeof(int), SEC_LEN-1, dat);
     };
     fseek(dat, head_len+sec_len, SEEK_SET);
    }
    else
     if(sm_lm == SM_SECTION)
     {
      i = ftell(dat);
      fseek(dat, head_len+sec_len, SEEK_SET);
      slen = fread(&itmp, sizeof(int), 1, dat);
      if(slen == 1 && itmp == 1)
       fread(&(other_sec_data[1]), sizeof(int), SEC_LEN-1, dat);
      fseek(dat, i, SEEK_SET);
     }
     else
      panic("%s: unsupported section type(=%i)", __func__, sm_lm);
    /* <ISHERE> */
    i = ftell(dat);
    slen = fread(&itmp, sizeof(int), 1, dat);
    if(slen != 0)
    {
     if(itmp == 1)
     {
      warn("%s: <ISHERE> == 1", __func__);
      my_sec_data[SEC_IH] = 1;
      slen = fread(&(my_sec_data[1]), sizeof(int), SEC_LEN-1, dat);
      if(my_sec_data[SEC_IH] && my_sec_data[SEC_NEV] != evd->nev)
       panic("%s: envd->nev(=%i) != <NEV>(=%i)",__func__, evd->nev, my_sec_data[SEC_NEV]);
     }
     fseek(dat, i, SEEK_SET);
    }
    my_sec_data[SEC_IH] = itmp = 1;
    fwrite(&itmp, sizeof(int), 1, dat);
    /* <ISVAL> */
    my_sec_data[SEC_VAL] = itmp = 1;
    fwrite(&itmp, sizeof(int), 1, dat);
    /* <ISVEC> */
    if(evd->val_vec == 0)
     my_sec_data[SEC_VEC] = itmp = 0;
    else
     my_sec_data[SEC_VEC] = itmp = 1;
    fwrite(&itmp, sizeof(int), 1, dat);
    /* NEV */
    my_sec_data[SEC_NEV] = itmp = evd->nev;
    fwrite(&itmp, sizeof(int), 1, dat);
    /* <TYPESIZE> */
    my_sec_data[SEC_TS] = itmp = sizeof(t_real);
    fwrite(&itmp, sizeof(int), 1, dat);
    /* <OFFSET_VAL> */
    /* if my data is alreade here */
    if (my_sec_data[SEC_VAL] && my_sec_data[SEC_OVAL])
     itmp = my_sec_data[SEC_OVAL];
    else
    {
     itmp = data_off;
     if(itmp <= my_sec_data[SEC_OVEC])
      itmp = my_sec_data[SEC_OVEC]+sizeof(t_complex)*VOL*NDIRAC*NCOLORS*my_sec_data[SEC_NEV];
     if(itmp <= other_sec_data[SEC_OVAL])
      itmp = other_sec_data[SEC_OVAL]+sizeof(t_complex)*other_sec_data[SEC_NEV];
     if(itmp <= other_sec_data[SEC_OVEC])
      itmp = other_sec_data[SEC_OVEC]+sizeof(t_complex)*VOL*NDIRAC*NCOLORS*other_sec_data[SEC_NEV];
    }
    fwrite(&itmp, sizeof(int), 1, dat);
    my_sec_data[SEC_OVAL] = itmp;
    /* <OFFSET_VEC> */
    if (evd->val_vec == 0)
     itmp = 0;
    else
    {
     /* if my data is alreade here */
     if (my_sec_data[SEC_VEC] && my_sec_data[SEC_OVEC])
      itmp = my_sec_data[SEC_OVEC];
     else
     {
      itmp = data_off;
      if(itmp <= my_sec_data[SEC_OVAL])
       itmp = my_sec_data[SEC_OVAL]+sizeof(t_complex)*my_sec_data[SEC_NEV];
      if(itmp <= other_sec_data[SEC_OVAL])
       itmp = other_sec_data[SEC_OVAL]+sizeof(t_complex)*other_sec_data[SEC_NEV];
      if(itmp <= other_sec_data[SEC_OVEC])
       itmp = other_sec_data[SEC_OVEC]+sizeof(t_complex)*VOL*NDIRAC*NCOLORS*other_sec_data[SEC_NEV];
     }
    }
    fwrite(&itmp, sizeof(int), 1, dat);
    my_sec_data[SEC_OVEC] = itmp;
    /* <EVAL> */
    if (evd->evals == NULL)
     panic("%s: evd->evals == NULL", __func__);
    fseek(dat, my_sec_data[SEC_OVAL], SEEK_SET);
    for(i = 0; i < evd->nev; i++)
    {
     t_real tre, tim;
     tre = CCALL(creal)(evd->evals[i]);
     tim = CCALL(cimag)(evd->evals[i]);
     fwrite(&tre, sizeof(t_real), 1, dat);
     fwrite(&tim, sizeof(t_real), 1, dat);
    }
    /* <EVEC> */
    fseek(dat, my_sec_data[SEC_OVEC], SEEK_SET);
    if (evd->val_vec == 1)
    {
     t_real tre, tim;
     for(i = 0; i < evd->nev; i++)
     {
#ifdef CHECK_NORM
      t_complex norm;
#endif
      for (j = 0; j < VOL*NDIRAC*NCOLORS; j++)
      {
       tre = CCALL(creal)(*(evd->evecs + i*VOL*NDIRAC*NCOLORS + j));
       tim = CCALL(cimag)(*(evd->evecs + i*VOL*NDIRAC*NCOLORS + j));
       fwrite(&tre, sizeof(t_real), 1, dat);
       fwrite(&tim, sizeof(t_real), 1, dat);
      }
#ifdef CHECK_NORM
      norm = innprod(evd->evecs+i*VOL*NDIRAC*NCOLORS, evd->evecs+i*VOL*NDIRAC*NCOLORS);
      fprintf(stderr, "write norm = %f + I %f\n", creal(norm), cimag(norm));
#endif
     }
    }
    fclose(dat);
} /*}}}*/

#undef  IDSTR_LEN
#undef  SEC_LEN
#undef  SEC_IH
#undef  SEC_VAL
#undef  SEC_VEC
#undef  SEC_NEV
#undef  SEC_TS
#undef  SEC_OVAL
#undef  SEC_OVEC

/*{{{*/
/*!
 * \brief Compare complex numbers by their magnitude
 *
 * \param a
 * \param b
 *
 * \return int (|a| < |b| ? -1 : (|a| ==|b| ? 0 : 1))
 */
/*}}}*/
int cmplx_cmp_mod(const void * a, const void * b) /*{{{*/
{
    t_real tmpa;
    t_real tmpb;

    tmpa = (*((t_complex*)a))*(CCALL(conj)(*((t_complex*)a)));
    tmpb = (*((t_complex*)b))*(CCALL(conj)(*((t_complex*)b)));

    return (tmpa < tmpb ? -1 : (tmpa == tmpb ? 0 : 1));
}/*}}}*/

/*{{{*/
//assumes o_vec are normalized
/*}}}*/
void proj_ortho_comp(t_complex* in, t_complex* evecs, int n_proj, t_complex* dp)/*{{{*/
{
    int i;
    t_complex dotpr;
    t_complex* tmp;

    for ( i = 0; i < n_proj; i++) {
        tmp = (t_complex*)(evecs+i*VOL*NDIRAC*NCOLORS);
        dotpr = innprod(tmp, in);
        if (dp != NULL)
            dp[i] = dotpr;
        xpcby(in, tmp, (t_complex)(-1.0*dotpr));
    }
}/*}}}*/



void scalar(void* data, t_gauge_vector * gv, t_cds_vector * out, t_cds_vector * in)/*{{{*/
{
    int idx;
    int idx_moved;
    int mu;
    int i;
    t_complex col_tmp[NDIRAC*NCOLORS] __attribute__ ((aligned(16)));

    for (i = 0; i < VOL*NDIRAC*NCOLORS; i++) {
        ((t_complex *)out)[i]  = +8.0*((t_complex*)in)[i];
    }
    for (idx = 0; idx < VOL; idx++) /*{{{*/
    {
        for (mu = 0; mu < DIM; mu++) {/*{{{*/
            idx_moved = (*lat_mov)[idx][mu][FWD];
            MATCALL(A_mul_CDSV)(col_tmp, (*gv)[idx][mu], &((*in)[idx_moved][0][0]));
#ifdef EFIELD
                for(i = 0; i<NDIRAC*NCOLORS; i++)
                 col_tmp[i] = col_tmp[i]*ev[idx][mu];
#endif

            for (i = 0; i < NDIRAC*NCOLORS; i++)
            {
                *((t_complex *)out + idx*NDIRAC*NCOLORS + i)  -= col_tmp[i];
            };
            idx_moved = (*lat_mov)[idx][mu][BWD];
            MATCALL(Ad_mul_CDSV)(col_tmp, (*gv)[idx_moved][mu], &((*in)[idx_moved][0][0]));
#ifdef EFIELD
                for(i = 0; i<NDIRAC*NCOLORS; i++)
                 col_tmp[i] = col_tmp[i]*conj(ev[idx_moved][mu]);
#endif
            for (i = 0; i < NDIRAC*NCOLORS; i++) {
                *((t_complex *)out + idx*NDIRAC*NCOLORS + i)  -= col_tmp[i];
            }
        }/*}}}*/
    }/*}}}*/
}/*}}}*/

int shumr(void *data, t_gauge_vector *gv, t_cds_vector *source, t_cds_vector *solution, t_real tol, int imax) //x is the solution
{
 t_ov_data_mm* ov_data;
 ov_data = (t_ov_data_mm*) data;
 t_real c1 = 1 + ov_data->mass;
 t_real c2 = 1 - ov_data->mass;

 int ret = 0;

 t_complex c_k, s_k, t11, mu, nu, mu_k, omega, theta, gamma, L11, iqv, eps;
 t_real L21, rnorm_p;

 t_cds_vector *q_tilde, *Dp, *v, *w, *s, *p, *q, *r, *rp, *xp, *q_old, *v_old;
 t_cds_vector *w_old, *s_old, *p1, *p2, *Dp1, *Dp2;

 if(posix_memalign((void*)&q_tilde, 16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for q_tilde (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&Dp,    16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for Dp      (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&v,     16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for v       (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&w,     16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for w       (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&s,     16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for s       (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&p,     16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for p       (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&q,     16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for q       (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&r,     16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for r       (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&rp,    16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for rp      (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&xp,    16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for xp      (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&q_old, 16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for q_old   (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&v_old, 16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for v_old   (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&w_old, 16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for w_old   (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&s_old, 16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for s_old   (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&p1,    16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for p1      (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&p2,    16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for p2      (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&Dp1,   16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for Dp1     (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&Dp2,   16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for Dp2     (%i Kb)", __func__, sizeof(t_cds_vector)/1024);

 set_zero((t_complex *)solution);

 memcpy(r, source, sizeof(t_cds_vector));

 t_real    rho   = vnorm((t_complex *)r);
 t_real    rnorm = rho;
 t_complex alpha = rho;

 memcpy(q, source, sizeof(t_cds_vector));
 ax((t_complex *)q, 1.0/rho);

 t_complex u12       = 0.0 + I*0.0;
 t_complex beta      = 1.0 + I*0.0;
 t_complex L11_tilde = 1.0 + I*0.0;
 t_complex c_km1 = 1.0 + I*0.0;
 t_complex s_km1 = 0.0 + I*0.0;
 t_complex c_km2 = 0.0 + I*0.0;
 t_complex s_km2 = 0.0 + I*0.0;

 set_zero((t_complex *)q_old);
 set_zero((t_complex *)v_old);
 set_zero((t_complex *)w_old);
 set_zero((t_complex *)s_old);
 set_zero((t_complex *)xp);
 set_zero((t_complex *)p1);
 set_zero((t_complex *)p2);
 set_zero((t_complex *)Dp1);
 set_zero((t_complex *)Dp2);

 int counter = 0;

 while((rnorm > tol) && (counter<=imax))
 {
  //V_op(q, v);
  ov_dirac_unitary(data, gv, v, q);
  if(counter > 0)
  {
   iqv = innprod((t_complex *)q_old, (t_complex *)v_old);
   //if(CCALL(cabs)(iqv)>1E-10) //TODO:  compare with tol?
   u12 = - 1.0*innprod((t_complex *)q_old, (t_complex *)v)/iqv;
   /*else
   {
    printf("iqv = 0!!!\n");
    return;
   };*/
  };
  gamma = -1.0*c1*u12;
  L11 = innprod((t_complex *)q, (t_complex *)v) + u12*innprod((t_complex *)q, (t_complex *)v_old);
  memcpy(q_tilde, v, sizeof(t_cds_vector));
  xpcby((t_complex *)q_tilde,     (t_complex *)q, -1.0*L11);
  xpcby((t_complex *)q_tilde, (t_complex *)v_old, u12);
//L21 is real
  L21 = vnorm((t_complex *)q_tilde);
  if(L21 <= tol)
  {
   printf("SHUMR stopped at iteration %i: L21 = %4.6E <= tol = %4.6E\n", counter, L21, tol);
   ret = 1;
   break;
  };

//w = q + qold*u12 + wold*gamma/L11tilde;
  memcpy(w, q, sizeof(t_cds_vector));
  xpcby((t_complex *)w, (t_complex *)q_old, u12);
  xpcby((t_complex *)w, (t_complex *)w_old, gamma/L11_tilde);

//s=c1*(q+q_old*u12)+c2*(v+v_old*u12)+s_old*gamma/L11_tilde;
  memcpy(s, q, sizeof(t_cds_vector));
  axc((t_complex *)s, c1);
  xpcby((t_complex *)s, (t_complex *)q_old, c1*u12);
  xpcby((t_complex *)s,     (t_complex *)v,     c2);
  xpcby((t_complex *)s, (t_complex *)v_old, c2*u12);
  xpcby((t_complex *)s, (t_complex *)s_old, gamma/L11_tilde);

  L11_tilde = c1 + c2*L11 - beta*gamma/L11_tilde;
  alpha     = alpha*beta/L11_tilde;

//x=x+w*alpha;
//r=r-s*alpha;
  xpcby((t_complex *)solution, (t_complex *)w,      alpha);
  xpcby(       (t_complex *)r, (t_complex *)s, -1.0*alpha);
//q_old=q; v_old=v; w_old=w; s_old=s;
  memcpy(q_old, q, sizeof(t_cds_vector));
  memcpy(v_old, v, sizeof(t_cds_vector));
  memcpy(w_old, w, sizeof(t_cds_vector));
  memcpy(s_old, s, sizeof(t_cds_vector));

//q=q_tilde/L21;
  memcpy(q, q_tilde, sizeof(t_cds_vector));
  axc((t_complex *)q, 1.0/L21);

  beta = -1.0*c2*L21;
  rnorm = vnorm((t_complex *)r);

  t11 = c1 + c2*L11;
  mu  = t11*c_km1 + gamma*CCALL(conj)(s_km1)*c_km2;
  nu  = c2*L21;
  if(CCALL(cabs)(mu) > 1E-10) //TODO: tolerance???
  {
    c_k = CCALL(cabs)(mu)/sqrt(CCALL(cabs)(mu)*CCALL(cabs)(mu) + CCALL(cabs)(nu)*CCALL(cabs)(nu));
    s_k = CCALL(conj)(c_k*nu/mu);
  }
  else
  {
    printf("SHUMR: Exceptional |mu| = %6.6E encountered!!!\n", CCALL(cabs)(mu));
    c_k = 0.0 + I*0.0;
    s_k = 1.0 + I*0.0;
  };
  omega = nu*alpha*s_k;
  mu_k  = c_k*mu + s_k*nu;
  eps   = t11*s_km1 - gamma*c_km1*c_km2;
  theta = -1.0*gamma*s_km2;
  //p = q/mu_k + q_old*u12/mu_k - p1*eps/mu_k - p2*theta/mu_k;
  memcpy(p, q, sizeof(t_cds_vector));
  axc((t_complex *)p, 1.0/mu_k);
  xpcby((t_complex *)p, (t_complex *)q_old, u12/mu_k);
  xpcby((t_complex *)p,    (t_complex *)p1, -1.0*eps/mu_k);
  xpcby((t_complex *)p,    (t_complex *)p2, -1.0*theta/mu_k);
  //Dp = c1/mu_k*q + c1*u12/mu_k*q_old + c2/mu_k*v + c2*u12/mu_k*v_old - Dp1/mu_k*eps - Dp2*theta/mu_k;
  memcpy(Dp, q, sizeof(t_cds_vector));
  axc((t_complex *)Dp, c1/mu_k);
  xpcby((t_complex *)Dp, (t_complex *)q_old,     c1*u12/mu_k);
  xpcby((t_complex *)Dp,     (t_complex *)v,         c2/mu_k);
  xpcby((t_complex *)Dp, (t_complex *)v_old,     c2*u12/mu_k);
  xpcby((t_complex *)Dp,   (t_complex *)Dp1,   -1.0*eps/mu_k);
  xpcby((t_complex *)Dp,   (t_complex *)Dp2, -1.0*theta/mu_k);

  //rnorm_p = vnorm(r + omega*Dp);
  memcpy(rp, r, sizeof(t_cds_vector));
  xpcby((t_complex *)rp, (t_complex *)Dp, omega);
  rnorm_p = vnorm((t_complex *)rp);

  //xp=x-omega*p;
  memcpy(xp, solution, sizeof(t_cds_vector));
  xpcby((t_complex *)xp, (t_complex *)p, -1.0*omega);

  c_km2 = c_km1;
  s_km2 = s_km1;
  c_km1 = c_k;
  s_km1 = s_k;
  memcpy( p2,  p1, sizeof(t_cds_vector));
  memcpy(Dp2, Dp1, sizeof(t_cds_vector));
  memcpy( p1,   p, sizeof(t_cds_vector));
  memcpy(Dp1,  Dp, sizeof(t_cds_vector));

  counter++;
  printf("SHUMR: iteration %i, rnorm = %6.6E\n", counter, rnorm);
 };
 if(counter>=imax)
  ret = 2;
  
 free(q_tilde);
 free(Dp); 
 free(v);
 free(w);
 free(s);
 free(p);
 free(q);
 free(r);
 free(rp); 
 free(xp);
 free(q_old); 
 free(v_old);
 free(w_old);
 free(s_old);
 free(p1);
 free(p2);
 free(Dp1); 
 free(Dp2);
  
 return ret;
}

