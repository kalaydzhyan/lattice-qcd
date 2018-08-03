/*{{{*/
/*!
 * \file polacc.c
 *
 * \brief Polynomial acceleration.
 *
 * References: 
 * - Y.Saad. Numerical methods for large eigenvalue problems
 * - H.Neff et.al. On the low fermionic eigenmode... [hep-lat/0106016]
 * 
 * The idea for any Hermitian operator \f$ H \f$ is very simple. The 
 * spectrum of operator \f$ H \f$ is real. Suppose that
 * \f[
 * spec(H) = [-l, -s] \cup [s, l].
 * \f]
 * We are interested in small part of the spectrum in the vicinity
 * of \f$ -s \f$
 * \f[
 * spec(H_{int}) = [-m, -s] \cup [s, m].
 * \f]
 * Then we construct a mappinh which maps an uninteresting region
 * \f$ spec(H_{unint}) = [-l, -m] \cup [m, l] \f$ to the 
 * interval \f$ [-1;1] \f$ but interesting region is mapped outside
 * of this interval.
 * \f[
 * p(x) =\frac{1}{2}(A x^2 + B) = \frac{1}{2}\left(\frac{4}{l^2-m^2} x ^2
 * - \frac{2(l^2+m^2)}{l^2-m^2}\right).
 * \f]
 * \image html polacc_1.jpg
 * \image latex polacc_1.eps
 * Then we construct a polynomial which is small in the interval 
 * \f$ [-1,1] \f$ and large outside of this interval. The best choice is
 * Chebyshev polynomial of some degree.
 * \image html polacc_2.jpg
 * \image latex polacc_2.eps
 * 
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \date   ‘Œ ÔÀ‘ 29 09:08:31 MSD 2004
 *//*}}}*/
#define _XOPEN_SOURCE	600
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "panic.h"
#include "types.h"
#include "linal.h"
#include "polacc.h"

/*{{{*/
/*!
 * \brief Define coefficients of mapping \f$ p(x) \f$
 *
 * \param m
 * \param l
 * \param data
 *
 * \return void 
 *//*}}}*/
void build_mapping(double m, double l, t_polacc_data* data) /*{{{*/
{
  data->cf[0] = -2.0*((l*l + m*m)/(l*l - m*m));
  data->cf[1] = 4.0/(l*l - m*m);
}/*}}}*/

/*{{{*/
/*!
 * \brief Define coefficients of ellips
 *
 * \param lambda1 Estimate of smallest eigenvalue
 * \param xb	  See picture in polacc_ov()
 * \param data
 *
 * \return void 
 */
/*}}}*/
void build_mapping_ov(double lambda1, double xb, t_polacc_data* data) /*{{{*/
{
	double const R = 2;
	double const xc = 1;
	
	data->c = xb + sqrt(R*R - 1 + (xb-xc)*(xb-xc));
	data->e = data->c/R;
	data->sigma1 = data->e/(lambda1 - data->c);
}/*}}}*/

/*{{{*/
/*!
 * \brief Do polynomial acceleration for Overlap operator
 *
 * \param data
 * \param gv
 * \param out
 * \param in
 *
 * \return void 
 */
/*}}}*/
void pol_acc_ov(void* data, t_gauge_vector* gv, t_cds_vector* out, t_cds_vector* in) /*{{{*/
{
	int i, j;
	t_cds_vector* psi[3];
	t_cds_vector* psi_tmp;
	t_polacc_data* pad;
	double sigmak, sigmakp1;

	pad = (t_polacc_data*)data;
	for (i = 0; i < 3; i++)
		if (posix_memalign((void*)&(psi[i]), 16, sizeof(t_cds_vector)))
			panic("%s: can't allocate memory for psi[%i] (% Kb)", __func__, i,
					sizeof(t_cds_vector)/1024);
	xcpy((t_complex*)(psi[0]), (t_complex*)(in));
	/* psi[1] = sigma_1/e A*z0 (as in Y.Saad) */
	(*(pad->op))((void*)(pad->data), gv, psi[1], in);
	for (j = 0; j < VOL*NDIRAC*NCOLORS; j++)
		((t_complex*)(psi[1]))[j] *= pad->sigma1/pad->e;
	/* z1 = psi[1] = psi[1] - c/e z0 */
	xpby((t_complex*)(psi[1]), (t_complex*)(psi[0]), -pad->sigma1*pad->c/pad->e);
	sigmak = pad->sigma1;
	for (i = 1; i < pad->deg; i++) {
		sigmakp1 = 1.0/((2.0/pad->sigma1) - sigmak);
		(*(pad->op))((void*)(pad->data), gv, psi[2], psi[1]);
		for (j = 0; j < VOL*NDIRAC*NCOLORS; j++)
			((t_complex*)(psi[2]))[j] *= 2.0*sigmakp1/pad->e;
		xpby((t_complex*)(psi[2]), (t_complex*)(psi[1]), -2.0*sigmakp1*pad->c/pad->e);
		xpby((t_complex*)(psi[2]), (t_complex*)(psi[0]), -sigmakp1*sigmak);
		psi_tmp = psi[0];
		psi[0] = psi[1];	/* T_{n-1}(p(x)) */
		psi[1] = psi[2];	/* T_n(p(x)) */
		psi[2] = psi_tmp;	/* T_{n+1}(p(x)) */
		sigmak = sigmakp1;
	}
	xcpy((t_complex*)(out), (t_complex*)(psi[1]));
	for (i = 0; i < 3; i++)
		free(psi[i]);
}/*}}}*/

/*{{{*/
/*!
 * \brief Do Chebyshev iteration
 *
 * Evaluate \f$ \Psi_{out} = T_d(p(O)) \Psi_{in} \f$
 * 
 * \param gv	gauge field
 * \param out	\f$ \Psi_{out} \f$
 * \param in	\f$ \Psi_{in} \f$
 * \param  data \f$ O \f$ = data->op, \f$d\f$ = data-deg
 * 
 * \return void 
 *//*}}}*/
void pol_acc(void* data, t_gauge_vector* gv, t_cds_vector* out, t_cds_vector* in) /*{{{*/
{
	int i;
	t_cds_vector* psi[3];
	t_cds_vector* psi_tmp;
	t_polacc_data* pad;
	
	pad = (t_polacc_data*)data;
	for (i = 0; i < 3; i++)
		if (posix_memalign((void*)&(psi[i]), 16, sizeof(t_cds_vector)))
			panic("%s: can't allocate memory for psi[%i] (% Kb)", __func__, i,
					sizeof(t_cds_vector)/1024);
	/* calculate p(op) in */	
	xcpy((t_complex*)(psi[0]), (t_complex*)(in));
	(*(pad->op))((void*)(pad->data), gv, out, in);
	(*(pad->op))((void*)(pad->data), gv, psi[1], out);
	xcpy((t_complex*)(psi[2]), (t_complex*)(psi[1]));
	axpby((t_complex*)(psi[2]), (t_real)(0.5*pad->cf[1]),
			(t_complex*)(psi[0]), (t_real)(0.5*pad->cf[0]));
	/* Recursion for Chebyshev polynomial of degree pad->deg */
	xcpy((t_complex*)psi[1], (t_complex*)in);
	/*
	 * psi[0] = in
	 * psi[1] = in
	 * psi[2] = p(op) in
	 */
	for (i = 1; i < pad->deg; i++) {
		psi_tmp = psi[0];
		psi[0] = psi[1];	/* T_{n-1}(p(x)) */
		psi[1] = psi[2];	/* T_n(p(x)) */
		psi[2] = psi_tmp;	/* T_{n+1}(p(x)) */
		/* psi[2] = 2 p(x) psi[1] = 2 p(x) T_n(p(x)) */
		(*(pad->op))(pad->data, gv, out, psi[1]);
		(*(pad->op))(pad->data, gv, psi[2], out);
		axpby((t_complex*)(psi[2]), (t_real)(pad->cf[1]),
				(t_complex*)(psi[1]), (t_real)(pad->cf[0]));
		xpby((t_complex*)(psi[2]), (t_complex*)(psi[0]), (t_real)(-1.0));
	}
	xcpy((t_complex*)(out), (t_complex*)(psi[2]));
	for (i = 0; i < 3; i++)
		free(psi[i]);
}/*}}}*/

/*{{{*/
/*!
 * \brief Calculate original eigenvalues
 *
 * Here we are using the fact that if \f$ |\Psi > \f$ is the
 * eigenfunction of some operator \f$ Q \f$, then it is also the eigenfunction
 * of some polynomial of \f$ Q \f$, i.e.
 * \f$ p(Q) |\Psi> = p(\lambda) |\Psi>\f$.
 *
 * To extract the eigenvalue of \f$ Q \f$ we use a following relation
 *
 * \f$ \lambda = \frac{<\Psi|Q|\Psi>}{<\Psi|\Psi>}. \f$
 * 
 * \param gv gauge configuration
 * \param oev
 *
 * \return void 
 *//*}}}*/
void orig_eval_hw(t_gauge_vector* gv, t_op_ev* oev) /*{{{*/
{
	int i;
	t_cds_vector* tmp;
	t_polacc_data* pad;
 	t_cds_vector* cdstmp;
	
	pad = (t_polacc_data*)oev->data;
	if (posix_memalign((void*)&(tmp), 16, sizeof(t_cds_vector)))
		panic("%s: can't allocate memory for tmp (% Kb)", __func__,
				sizeof(t_cds_vector)/1024);
	for (i = 0; i < oev->nev; i++) {
		cdstmp =(t_cds_vector*)(oev->evecs+i*VOL*NDIRAC*NCOLORS);
		(*(pad->op))(pad->data, gv, tmp, cdstmp);
    	oev->evals[i] = innprod((t_complex*)cdstmp, (t_complex*)tmp);
		/* TODO: residual calculation */
	}
	free(tmp);
}/*}}}*/
