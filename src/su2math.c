/*{{{*/
/*!
 * \file su2math.c
 *
 * \brief SU2 global constants and routines
 *
 *
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \date   Чтв Сен 23 12:11:22 MSD 2004
 *//*}}}*/
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "types.h"
#include "su2math.h"
#include "mt19937ar.h"

/*! I_SIGMA1_SU2 =
 * \f$ i \sigma_1 =
 * \left(\begin{tabular}{cc} 0 & i \\ i & 0 \\ \end{tabular}\right)\f$
 */
t_gmat I_SIGMA1_SU2 = {0.0 + I*0.0, 0.0 + I*1.0};
/*! I_SIGMA1_SU2 =
 * \f$ i \sigma_2 =
 * \left(\begin{tabular}{cc} 0 & 1 \\ -1 & 0 \\ \end{tabular}\right)\f$
 */
t_gmat I_SIGMA2_SU2 = {0.0 + I*0.0, 1.0 + I*0.0};
/*! I_SIGMA1_SU2 =
 * \f$ i \sigma_3 =
 * \left(\begin{tabular}{cc} i & 0 \\ 0 & -i \\ \end{tabular}\right)\f$
 */
t_gmat I_SIGMA3_SU2 = {0.0 + I*1.0, 0.0 + I*0.0};

#if (EXP2_SSE & SINGLE_PREC)
static float patd[4] __attribute__ ((aligned(16))) = {1.0,-1.0,-1.0,1.0};
#endif

/*{{{*/
/*!
 * \brief	Multiply ColorDiracSpace vector by \f$SU2^{\dagger}\f$ matrix
 *
 * \param _A t_gmat
 * \param _in (t_complex *)
 * \param _out (t_complex *)
 * 
 */
/*}}}*/
inline void Ad_mul_CDSV_SU2(t_complex* out, t_gmat A, t_complex* in) /*{{{*/
{
	int i;
#if (EXP2_SSE & SINGLE_PREC)
	for (i = 0; i < NDIRAC; i++) {
	__asm__ __volatile__ (
			"movlps	(%2,%3,8), %%xmm1			/* p0r, p0i -> lodword(xmm1) */\n\t"
			"addl	%4, %3						/* %3 + NDIRAC */\n\t"
			"movhps	(%2,%3,8), %%xmm1			/* p1r, p1i -> hidword(xmm1) */\n\t"
			"/* xmm1 = p0r, p0i, p1r, p1i */\n\t"
			"movss	(%1), %%xmm0				/* ar -> xmm0 */\n\t"
			"shufps	$0, %%xmm0, %%xmm0\n\t"
			"movaps	%%xmm1, %%xmm6\n\t"
			"mulps	%%xmm0, %%xmm6				/* xmm6 = ar * Psi */\n\t"
			"movss	4(%1), %%xmm0				/* xmm0 = ai */\n\t"
			"shufps	$0, %%xmm0, %%xmm0\n\t"
			"movaps (%5), %%xmm2				/* xmm2 = {1,-1,-1,1} */\n\t"
			"mulps	%%xmm2, %%xmm0				/* xmm0 = {ai,-ai,-ai,ai} */\n\t"
			"shufps	$177, %%xmm1, %%xmm1		/* xmm1 = {p0i, p0r, p1i, p1r} */\n\t"
			"movaps	%%xmm1, %%xmm7\n\t"
			"mulps	%%xmm0, %%xmm7				/* xmm7 = 'ai'*xmm1 */\n\t"
			"addps	%%xmm6, %%xmm7\n\t"
			"movss	8(%1), %%xmm0				/* xmm0 = br */\n\t"
			"shufps	$0, %%xmm0, %%xmm0\n\t"
			"shufps	$5, %%xmm2, %%xmm2			/* xmm2 = {-1,-1,1,1} */\n\t"
			"mulps	%%xmm2, %%xmm0				/* xmm0 = {-br,-br,br,br} */\n\t"
			"shufps	$27, %%xmm1, %%xmm1			/* xmm1 = {p1r, p1i, p0r, p0i}\n\t"
			"movaps	%%xmm1, %%xmm6\n\t"
			"mulps	%%xmm0,	%%xmm6\n\t"
			"addps	%%xmm6,	%%xmm7\n\t"
			"movaps	12(%1), %%xmm0				/* xmm0 = bi */\n\t"
			"shufps	$0, %%xmm0, %%xmm0\n\t"
			"shufps	$34, %%xmm2, %%xmm2			/* xmm2 = {1,-1,1,-1} */\n\t"
			"mulps	%%xmm2, %%xmm0				/* xmm0 = {bi,-bi,bi,-bi} */\n\t"
			"shufps	$177, %%xmm1, %%xmm1		/* xmm1 = {p1i, p1r, p0i, p0r} */\n\t"
			"mulps	%%xmm0, %%xmm1\n\t"
			"addps	%%xmm1, %%xmm7\n\t"
			"movl	%0, %%edi\n\t"
			"movhps	%%xmm7, (%%edi,%3,8)\n\t"
			"subl	%4, %3\n\t"
			"movlps	%%xmm7, (%%edi,%3,8)"
			: "=m" (out)
			: "r" (A), "r" (in), "r" (i), "i" (NDIRAC), "r" (patd)
			: "%xmm0", "%xmm1", "%xmm2", "%xmm6", "%xmm7", "%edi"
			);
	}
#else
	for (i = 0; i < NDIRAC; i++) {									
		out[i] = CCALL(conj)(A[0])*(in[i]) - A[1]*(in[i+NDIRAC]);
		out[i+NDIRAC] = CCALL(conj)(A[1])*(in[i]) +	A[0]*(in[i+NDIRAC]);									
	}																	
#endif
}/*}}}*/

#if (EXP2_SSE & SINGLE_PREC)
static float pat[4] __attribute__ ((aligned(16))) = {-1.0, 1.0, 1.0,-1.0};
#endif

/*{{{*/
/*!
 * \brief	Multiply ColorDiracSpace vector by \f$SU2\f$ matrix
 *
 * \param A t_gmat
 * \param in (t_complex *)
 * \param out (t_complex *)
 * 
 */
/*}}}*/
inline void A_mul_CDSV_SU2(t_complex* out, t_gmat A, t_complex* in) /*{{{*/
{
	int i;
#if (EXP2_SSE && SINGLE_PREC)
	for (i = 0; i < NDIRAC; i++) {
	__asm__ __volatile__ (
			"movlps	(%2,%3,8), %%xmm1			/* p0r, p0i -> lodword(xmm1) */\n\t"
			"addl	%4, %3						/* %3 + NDIRAC */\n\t"
			"movhps	(%2,%3,8), %%xmm1			/* p1r, p1i -> hidword(xmm1) */\n\t"
			"/* xmm1 = p0r, p0i, p1r, p1i */\n\t"
			"movss	(%1), %%xmm0				/* ar -> xmm0 */\n\t"
			"shufps	$0, %%xmm0, %%xmm0\n\t"
			"movaps	%%xmm1, %%xmm6\n\t"
			"mulps	%%xmm0, %%xmm6				/* xmm6 = ar * Psi */\n\t"
			"movss	4(%1), %%xmm0				/* xmm0 = ai */\n\t"
			"shufps	$0, %%xmm0, %%xmm0\n\t"
			"movaps (%5), %%xmm2				/* xmm2 = {-1, 1, 1,-1} */\n\t"
			"mulps	%%xmm2, %%xmm0				/* xmm0 = {-ai, ai, ai,-ai} */\n\t"
			"shufps	$177, %%xmm1, %%xmm1		/* xmm1 = {p0i, p0r, p1i, p1r} */\n\t"
			"movaps	%%xmm1, %%xmm7\n\t"
			"mulps	%%xmm0, %%xmm7				/* xmm7 = 'ai'*xmm1 */\n\t"
			"addps	%%xmm6, %%xmm7\n\t"
			"movss	8(%1), %%xmm0				/* xmm0 = br */\n\t"
			"shufps	$0, %%xmm0, %%xmm0\n\t"
			"shufps	$5, %%xmm2, %%xmm2			/* xmm2 = {1,1,-1,-1} */\n\t"
			"mulps	%%xmm2, %%xmm0				/* xmm0 = {br,br,-br,br} */\n\t"
			"shufps	$27, %%xmm1, %%xmm1			/* xmm1 = {p1r, p1i, p0r, p0i}\n\t"
			"movaps	%%xmm1, %%xmm6\n\t"
			"mulps	%%xmm0,	%%xmm6\n\t"
			"addps	%%xmm6,	%%xmm7\n\t"
			"movaps	12(%1), %%xmm0				/* xmm0 = bi */\n\t"
			"shufps	$0, %%xmm0, %%xmm0\n\t"
			"shufps	$34, %%xmm2, %%xmm2			/* xmm2 = {-1,1,-1,1} */\n\t"
			"mulps	%%xmm2, %%xmm0				/* xmm0 = {-bi,bi,-bi,bi} */\n\t"
			"shufps	$177, %%xmm1, %%xmm1		/* xmm1 = {p1i, p1r, p0i, p0r} */\n\t"
			"mulps	%%xmm0, %%xmm1\n\t"
			"addps	%%xmm1, %%xmm7\n\t"
			"movl	%0, %%edi\n\t"
			"movhps	%%xmm7, (%%edi,%3,8)\n\t"
			"subl	%4, %3\n\t"
			"movlps	%%xmm7, (%%edi,%3,8)"
			: "=m" (out)
			: "r" (A), "r" (in), "r" (i), "i" (NDIRAC), "r" (pat)
			: "%xmm0", "%xmm1", "%xmm2", "%xmm6", "%xmm7", "%edi"
			);
	}
#else
	for (i = 0; i < NDIRAC; i++) {
		out[i] = A[0]*(in[i]) + A[1]*(in[i+NDIRAC]);
		out[i+NDIRAC] = -CCALL(conj)(A[1])*(in[i]) + CCALL(conj)(A[0])*(in[i+NDIRAC]);
	}
#endif
}/*}}}*/

/*{{{*/
/*!
 * \brief Print \f$SU(2)\f$ matrix to stdout 
 *
 * \param m
 *
 * \return void 
 *//*}}}*/
void print_mat_SU2(t_complex * m)/*{{{*/
{
	fprintf(stdout, "\n");
	fprintf(stdout, "(%s%1.6f %s i %1.6f, %s%1.6f %s i %1.6f)\n",
		(CCALL(creal)(m[0]) >= 0.0 ? "+" : "-"),
		CCALL(fabs)(CCALL(creal)(m[0])),
		(CCALL(cimag)(m[0]) >= 0.0 ? "+" : "-"),
		CCALL(fabs)(CCALL(cimag)(m[0])),
		(CCALL(creal)(m[1]) >= 0.0 ? "+" : "-"),
		CCALL(fabs)(CCALL(creal)(m[1])),
		(CCALL(cimag)(m[1]) >= 0.0 ? "+" : "-"),
		CCALL(fabs)(CCALL(cimag)(m[1])));
	fprintf(stdout, "(%s%1.6f %s i %1.6f, %s%1.6f %s i %1.6f)\n",
		(CCALL(creal)(-m[1]) >= 0.0 ? "+" : "-"),
		CCALL(fabs)(CCALL(creal)(m[1])),
		(CCALL(cimag)(m[1]) >= 0.0 ? "+" : "-"),
		CCALL(fabs)(CCALL(cimag)(m[1])),
		(CCALL(creal)(m[0]) >= 0.0 ? "+" : "-"),
		CCALL(fabs)(CCALL(creal)(m[0])),
		(CCALL(cimag)(-m[0]) >= 0.0 ? "+" : "-"),
		CCALL(fabs)(CCALL(cimag)(m[0])));
}/*}}}*/

void print_part_cdsv(t_complex *v) /*{{{*/
{
	int i;
	
	fprintf(stderr, "(\n");
	for (i = 0; i < NCOLORS*DIM;i++)
		fprintf(stderr, "\t[%f + I*%f]\n", CCALL(creal)(*(v+i)), CCALL(cimag)(*(v+i)));
	fprintf(stderr, ")\n");
}/*}}}*/

void mulcdv_test() /*{{{*/
{
	t_complex * cdsv1;
	t_complex * cdsv2;
	t_gmat t;
	int i;
	
	cdsv1 = malloc(sizeof(t_complex)*NCOLORS*DIM);
	cdsv2 = malloc(sizeof(t_complex)*NCOLORS*DIM);
	memset(cdsv1, 0, sizeof(t_complex)*NCOLORS*DIM);
	for (i = 0; i < NCOLORS*DIM; i++)
		*(cdsv1+i) = 1.0;
	memset(cdsv2, 0, sizeof(t_complex)*NCOLORS*DIM);
	print_part_cdsv(cdsv1);
	MATCALL(rand_A)(t);
//	MATCALL(print_mat)(t);
//	MATCALL(Ad_mul_CDSV)(cdsv2, t, cdsv1);
//	gampr_mul_CDSV(3, 1, cdsv2, cdsv1);
	print_part_cdsv(cdsv2);
	free(cdsv1);
	free(cdsv2);
}/*}}}*/
