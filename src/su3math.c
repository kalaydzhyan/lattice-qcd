/*{{{*/
/*!
 * \file su3math.c
 *
 * \brief SU3 global constants and routines
 *
 *
 * $Id$
 *
 * \author Pavel Buividovich, email: gbuividovich@gmail.com (implemented background magnetic field, chemical potential, SU(3) gauge group in 2008 - 2009)
 *//*}}}*/
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "types.h"
#include "su3math.h"
#include "mt19937ar.h"

/*{{{*/
/*!
 * \brief	Multiply ColorDiracSpace vector by \f$SU3^{\dagger}\f$ matrix
 *
 * \param _A t_gmat
 * \param _in (t_complex *)
 * \param _out (t_complex *)
 * 
 */
/*}}}*/

inline void Ad_mul_CDSV_SU3(t_complex* out, t_gmat A, t_complex* in) /*{{{*/
{
 int i, ic, ic1;
 for(i = 0; i<NDIRAC; i++)
  for(ic = 0; ic < NCOLORS; ic++)
  {									
   out[i + NDIRAC*ic] = 0;
   for(ic1 = 0; ic1 < NCOLORS; ic1++)
    out[i + NDIRAC*ic] += CCALL(conj)(A[NCOLORS*ic1 + ic])*in[i + NDIRAC*ic1];    
  }; 
}/*}}}*/

/*{{{*/
/*!
 * \brief	Multiply ColorDiracSpace vector by \f$SU3\f$ matrix
 *
 * \param A t_gmat
 * \param in (t_complex *)
 * \param out (t_complex *)
 * 
 */
/*}}}*/
inline void A_mul_CDSV_SU3(t_complex* out, t_gmat A, t_complex* in) /*{{{*/
{
 int i, ic, ic1;
 for(i = 0; i<NDIRAC; i++)
  for(ic = 0; ic < NCOLORS; ic++)
  {									
   out[i + NDIRAC*ic] = 0;
   for(ic1 = 0; ic1 < NCOLORS; ic1++)
    out[i + NDIRAC*ic] += A[NCOLORS*ic + ic1]*in[i + NDIRAC*ic1];
  }; 
}/*}}}*/

/*{{{*/
/*!
 * \brief Print \f$SU(3)\f$ matrix to stdout 
 *
 * \param m
 *
 * \return void 
 *//*}}}*/
void print_mat_SU3(t_complex * m)/*{{{*/
{
    int i, j; 
	fprintf(stdout, "\n");
	for(i=0; i<NCOLORS; i++)
	{
     for(j=0; j<NCOLORS; j++)
      fprintf(stdout, "(%2.4lf + I*%2.4lf)   ", (double)CCALL(creal)(m[i*NCOLORS + j]), (double)CCALL(cimag)(m[i*NCOLORS + j]) );
     fprintf(stdout, "\n");                    
    };
	fprintf(stdout, "\n");
	fflush(stdout);
}/*}}}*/

void print_part_cdsv(t_complex *v) /*{{{*/
{
	int i;
	
	fprintf(stderr, "(\n");
	for (i = 0; i < NCOLORS*DIM;i++)
		fprintf(stderr, "\t[%f + I*%f]\n", CCALL(creal)(*(v+i)), CCALL(cimag)(*(v+i)));
	fprintf(stderr, ")\n");
}/*}}}*/

