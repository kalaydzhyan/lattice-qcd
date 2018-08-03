/*{{{*/
/*!
 * \file minmax.c
 *
 * \brief MinMax polynomial approximation of \f$sign(X)\f$ function
 * 
 * We want to find the best approximation of the \f$ sign(X) \f$ function 
 * in the range \f$ \sqrt{\epsilon} \le |X| \le 1 \f$.
 * One can write
 * \f[
 * sign(X) = \frac{X}{\sqrt{X^2}} = X \frac{1}{\sqrt{y}},
 * \f]
 * where \f$ y = X^2, \ \epsilon \le y \le 1 \f$. Now we will try
 * to find the best
 * approximation of the function \f$ g(y) = \frac{1}{\sqrt{y}} \f$.
 * 
 * We are looking for a polynomial
 * \f$ P_n(y),\ y \in [\epsilon;1],\ \f$ of degree \f$n\f$
 * that minimizes the maximal relative error
 * \f[
 * \delta = max |h(y)|,
 * \f]
 * where
 * \f[
 * h(y) = \frac{\frac{1}{\sqrt{y}} - P_n(y)}{\frac{1}{\sqrt{y} } } = 
 * 1 - \sqrt{y}P_n(y)
 * \f]
 * For reasons of numerical stability, it is advantageous to represent 
 * the polynomial in the form of a series
 * \f[
 * P_n(y) = \sum_{k=0}^n c_k T_k(z),\ z = \frac{2y-1-\epsilon}{1-\epsilon},
 * \f]
 * of Chebyshev polynomials, which are defined in the range \f$ [-1;1] \f$.
 * 
 * Minmax polynomial can be shown to exist and to be uniquely determined
 * provided the function that is to be approximated does not vanish.
 * To construct the minmax polynomial we should know (at least) \f$ n+2 \f$
 * extrema of error function of equal height and alternating sign. Evidently
 * it is not possible to know in advance where these extrema are located.
 * 
 * To find \f$P_n(y)\f$ we are using Maehly "First Direct Method".
 *
 * Stage 1.
 * For a given set of critical points
 * \f[
 * \{y_l\} = {\epsilon \le y_0 < y_1 < ... < y_{n+1} \le 1},
 * \f]
 * find \f$P_n(y)\f$ such that
 * \f[
 * h(y_l) = (-1)^l u,\ l = 0,...,n+1,
 * \f]
 * by solving just \f$ n+2 \f$ linear equations in the unknowns
 * \f$c_0,..,c_n,u\f$.
 * 
 * Stage 2.
 *
 * Find the extrema of \f$ h(y) \f$, computed in Stage 1. These extrema 
 * are now a new set of critical points for the subsequent Stage 1.
 *
 * As a start set of critical points we are using extrema of Chebyshev
 * polynomial of degree \f$ n+1 \f$, which has exactly \f$n+2\f$ extrema.
 * 
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \date   Втр Сен 28 16:47:57 MSD 2004
 *//*}}}*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#include "defs.h"
#include "minmax.h"
#include "panic.h"

/* wrapper of LAPACK routines */
/* computes a QR factorization of a general rectangular matrix. */
void UNDERSC(dgeqrf)(int*, int*, double*, int*, double*, double*, int*, int*);
/*
 * Multiplies a general matrix by the orthogonal matrix
 * from a QR factorization determined by DGEQRF
 */
void UNDERSC(dormqr)(char*, char*, int*, int*,  int*, double*, int*, double*,
		double*, int*, double*, int*, int*, int, int);
/*
 * Solves a triangular system of linear equations AX=B,
 * A**T X=B or A**H X=B
 */ 
void UNDERSC(dtrtrs)(char*, char*, char*, int*, int*, double*, int*, double*,
		int*, int*, int, int, int);

/*{{{*/
/*!
 * \brief Calculate zeroes of Chebyshev polynomial
 *
 * The polynomial \f$ T_n(x) \f$ has \f$n\f$ zeros in the interval
 * \f$[-1:1]\ \ \f$, and they are located at the points
 * \f[
 *  x = \cos\left(
 *  \frac{\pi(k+\frac{1}{2})}{n}
 *  \right), k = 0,1,...,n-1
 * \f] 
 * 
 * \param n degree of Chebyshev polynomial \f$ T_n(x) \f$
 * \param x zeros
 *
 * \return void 
 */
/*}}}*/
static void cheb_pol_zeros(int n, double * x) /*{{{*/
{
	int k;
	
	for (k = 0; k < n; k++)
		x[n-k-1] = cos(((double)k+.5)*M_PI/((double) n));
}/*}}}*/

/*{{{*/
/*!
 * \brief Calculate extrema of Chebyshev polynomial
 *
 * The polynomial \f$ T_n(x) \f$ has \f$n+1\f$ extrema (minima and maxima)
 * in the interval
 * \f$[-1:1]\ \ \f$, and they are located at the points
 * \f[
 *  x_k = \cos\left(\frac{\pi(k)}{n}\right), k = 0,1,...,n
 * \f]
 * \verbatim
 * |---------0---------|
 * ^        <-         ^
 * |                   |
 * x_n=-1              x_0=1
 * \endverbatim
 * We do not calculate \f$x_n\f$ and \f$x_0\f$
 *
 * \param n degree of Chebyshev polynomial \f$ T_n(x) \f$
 * \param x extrema. double x[n-1].
 *
 * \return void 
 */
/*}}}*/
static void cheb_pol_extr(int n, double* x)/*{{{*/
{
	int k;
	
	for (k = 1; k < n; k++)
 		x[n-k-1]=cos(((double)k)*M_PI/((double) n));
}/*}}}*/

/*{{{*/
/*!
 * \brief Evaluate function approximated by polynomial of degree \f$n\f$
 *
 * \f[
 * P_n(y) = \sum_{k=0}^n c_k T_k(z),\ \ \
 * z = (2y-1-\epsilon)/(1-\epsilon) 
 * \f]
 * We are using folowing recurrence relation
 * \f[
 * T_{n+1}(z) = 2zT_n(z)-T_{n-1}(z), n \ge 1
 * \f]
 * 
 * TODO Use Clenshaw's recurrence formula for better precision
 * 
 * \param eps approximation interval \f$[\epsilon, 1]\f$
 * \param n	degree of polynomial
 * \param c	coefficients of approximation \f$c_k,\ \ k=0,1,...,n\f$
 * \param y 	
 * \param C	if (C != NULL) then C[i] = \f$ c_i T_i(z) \f$
 *
 * \return double 
 */
/*}}}*/
double eval_pol_cheb(double eps, int n, double* c, double y, double* C) /*{{{*/
{
	int i;
	double val;		/* T_{k+1}(z) */
	double val_p;	/* T_k(z) */
	double val_pp;	/* T_{k-1}(z) */
	double res;
	double x;

	x = (2.0*y - 1.0 - eps)/(1.0 - eps);
	val_p = 1.0;	/* T_0(z) = 1.0 */
	if (C != NULL)
    	C[0] = c[0];
	val = x;		/* T_1(z) = z */
	res = c[0] + (n > 0 ?  c[1]*val : 0.0);
	if (C != NULL && n > 0)
		C[1]=c[1]*val;
	/* summation */
	for (i = 1; i < n; i++)
	{
		val_pp = val_p;
		val_p = val;
		val = 2.0*x*val_p - val_pp;
		res += c[i+1]*val;
		if (C != NULL)
			C[i+1] = c[i+1]*val;
	}
  
	return res;
} /*}}}*/

double eval_pol_cheb_cc(double eps, int n, double* c, double y) /*{{{*/
{
	int i;
	double sv;
	double x, x2;
	double d = 0.0, dd = 0.0;

	x = (2.0*y - 1.0 - eps)/(1.0 - eps); /* x \in [-1:1] */
	x2 = 2.0*x;
	/* summation */
	for (i = n-1; i >= 1; i--) {
		sv = d;
		d = x2*d - dd + c[i];
		dd = sv;
	}
  
	return (x*d - dd + c[0]);
} /*}}}*/


/*{{{*/
/*!
 * \brief Chebyshev coefficients of the derivative of the function whose
 * coefficients are c
 * 
 * By explicit calculation one obtains
 *	\f[
 *	c^\prime_n = 0, \ \ 
 *	c^\prime_{n-1} = 2nc_n,\ \  c^\prime_{n-2} = 2(n-1)c_{n-1},\ \
 *	c^\prime_{k-1} = c^\prime_{k+1} + 2(k+1)c_{k+1},\ \ 
 *	c^\prime_0 = \frac{1}{2} c^\prime_{k-1}{\big|}_{k=1}
 *	\f]
 *
 * \param eps functions \f$P(y)\f$ and \f$P^\prime(y)\f$ are defined on the interval
 * \f$[\epsilon,1]\f$
 * \param n	 degree of polynomial 
 * \param c coefficients of the polynomial \f$P(y)=\sum_{k=0}^n c_k T_k(z)\f$
 * \param cpr coefficients of the polynomial
 * \f$P^\prime(y) = \sum_{k=0}^{n-1} c^\prime_{k} T_k{z}\f$
 *
 * \return void 
 */
/*}}}*/
static void cheb_pol_der(double eps, int n, double* c, double* cpr) /*{{{*/
{
	int i;

	for (i = 0; i < n; i++)
		cpr[i] = 0.0;
	if (n > 0)
		cpr[n-1] = 2.0*n*c[n];
	if (n > 1)
		cpr[n-2]=2.0*(n-1)*c[n-1];
	for (i = n-3; i >= 0; i--)
		cpr[i] = cpr[i+2] + 2.0*((double)i+1.0)*c[i+1];
	cpr[0] *= 0.5;
	for (i = 0; i < n; i++)
		cpr[i] *= 2.0/(1.0 - eps);
#ifdef DEBUG_CHEB_POL_DER
	for (i = 0; i < n; i++)
		fprintf(stderr, "c[%i] = %f\tcpr[%i] = %f\n", i, c[i], i, cpr[i]);
#endif
}/*}}}*/


/*{{{*/
/*!
 * \brief Solve Ax = y using the QR decomposition of A as implemented in LAPACK
 *
 * \param n	dim(A)
 * \param A
 * \param y
 * \param x
 *
 * \return void 
 */
/*}}}*/
static void solv_leq_qr(int n, double* A, double* y, double* x) /*{{{*/
{
	int i, j;
	double* R;
	double* work,* tau;
	int workl, info, int1;
	char charL, charT, charU, charN ;
	
	charL = 'L';
	charT = 'T';
	charU = 'U';
	charN = 'N';
	int1 = 1;
	workl = 200*n;
	if ((work = malloc(sizeof(double)*workl)) == NULL)
		panic("solv_leq_qr: can't allocate memory for work (%i Kb)",
				sizeof(double)*workl);
	if ((tau = malloc(sizeof(double)*n)) == NULL)
		panic("solv_leq_qr: can't allocate memory for work (%i Kb)",
				sizeof(double)*n);
	if ((R = malloc(sizeof(double)*n*n)) == NULL)
		panic("solv_leq_qr: can't allocate memory for work (%i Kb)",
				sizeof(double)*n*n);
	/* Fortran stores matrix in a row-major order format */
  	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			R[n*j+i]=A[n*i+j];
	/* calls of LAPACK routines */
  	UNDERSC(dgeqrf)(&n, &n, R, &n, tau, work, &workl, &info);
	if (info != 0)
		panic("%s: error in dgeqrf (info = %i)", __func__, info);
	for (i = 0; i < n; i++)
		x[i] = y[i];
	UNDERSC(dormqr)(&charL, &charT, &n, &int1, &n, R, &n, tau, x, &n, work, &workl, &info, 1, 1);
	if (info != 0)
		panic("%s: error in dormqr (info = %i)", __func__, info);
	UNDERSC(dtrtrs)(&charU, &charN, &charN, &n, &int1, R, &n, x, &n, &info, 1, 1, 1);
	if (info != 0)
		panic("%s: error in dtrtrs (info = %i)", __func__, info);
	/* free memory */
	free(R);
	free(tau);
	free(work);
}/*}}}*/


/*{{{*/
/*!
 * \brief Stage 1 of Maehly algorithm.
 *
 * Find polynomial which satisfies
 * \f$ h(y_i) = (-1)^i u,\ l = 0, ..., n+1 \f$
 * We need to solve a system of \f$ n+2 \f$ linear equation in the
 * unknowns \f$ c_0, ..., c_n, u \f$
 *
 * \f[
 * h(y_i) = 1 - \sqrt{y_i} P_n(y_i) =
 * 1 - \sqrt{y_i} \sum_{k=0}^n c_k T_k(z_i) = (-1)^i u,
 * \f]
 * \f[
 * \left\{
 * \begin{tabular}{ccc}
 * $\sum_{k=0}^n \sqrt{y_0} c_k T_k(z_0) + (-1)^0 u$ & = & 1 \\
 * ... & & \\
 * $\sum_{k=0}^n \sqrt{y_{n+1}} c_k T_k(z_{n+1}) + (-1)^{n+1} u $& = & 1 \\
 * \end{tabular}
 * \right.
 * \f]
 * Or in a matrix form \f$ A x = r\f$
 * \f[
 * \left(
 * \begin{tabular}{cccc}
 * $\sqrt{y_0}$ & ... & $\sqrt{y_0}T_n(y_0)$ & $(-1)^0$ \\
 * & & ... &  \\
 * $\sqrt{y_{n}}$ & ... & $\sqrt{y_{n}}T_n(y_{n})$ & $(-1)^{n}$ \\
 * $\sqrt{y_{n+1}}$ & ... & $\sqrt{y_{n+1}}T_n(y_{n+1})$ & $(-1)^{n+1}$ \\
 * \end{tabular}
 * \right)
 * \left(
 * \begin{tabular}{c}
 * $c_0$ \\
 * ... \\
 * $c_n$ \\
 * $u$ 
 * \end{tabular}
 * \right) =
 * \left(
 * \begin{tabular}{c}
 * $1$ \\
 * ... \\
 * $1$ \\
 * $1$ 
 * \end{tabular}
 * \right)
 * \f]
 * 
 * \param eps
 * \param deg
 * \param c
 * \param y
 *
 * \return void 
 */
/*}}}*/
void extr_pol(double eps, int deg, double* c, double* y)/*{{{*/
{
	int i, j;
	double* A, * r, * cf;

	if ((A = malloc(sizeof(double)*(deg+2)*(deg+2))) == NULL)
		panic("extr_pol: can't allocate memory for matrix A (%i Kb)",
				sizeof(double)*(deg+2)*(deg+2)/1024);
	if ((r = malloc(sizeof(double)*(deg+2))) == NULL)
		panic("extr_pol: can't allocate memory for y (%i Kb)",
				sizeof(double)*(deg+2)/1024);
	if ((cf = malloc(sizeof(double)*(deg+2))) == NULL)
		panic("extr_pol: can't allocate memory for cf (%i Kb)",
				sizeof(double)*(deg+2)/1024);
	/* construct vector r {{{ */
	for (i = 0; i < deg+2; i++)
		r[i] = 1.0;
	/*}}}*/
	/* construct matrix A {{{*/
	for (i = 0; i < deg+2; i++)
	{
		/*
		 * here we use the ability of eavl_pol_cheb() to
		 * evaluete c_i T_i(y), i = 0, ..., n,
		 * with c_i = \sqrt{y}
		 */ 
		for (j = 0; j < deg+1; j++)
			cf[j] = sqrt(y[i]);
		eval_pol_cheb(eps, deg, cf, y[i], A + i*(deg + 2));
		A[i*(deg + 2) + deg + 1] = ( i%2 == 0 ? 1.0 : -1.0);
	}/*}}}*/
	solv_leq_qr(deg+2, A, r, c);
	/* free memory */
	free(cf);
	free(r);
	free(A);
}/*}}}*/


/*{{{*/
/*!
 * \brief Stage 2 of Maehly algorithm
 *
 * Locate extrema (new critical points) of error function
 *
 * TODO Better documentation and even understanding of this routine
 * 
 * \param eps 
 * \param n	degree of polynomial
 * \param c coefficients \f$c_k,\ k = 0, ..., n \f$
 * \param x old set of critical points \f$x_i,\ i = 0, ..., n+1 \f$
 * \param e	new set of critical points \f$e_i,\ i = 0, ..., n+1 \f$,
 * i.e. new extrema
 *
 * \return double maximal relative error
 */
/*}}}*/
static double locate_extr(double eps, int n, double* c, double* x, double* e) /*{{{*/
{
	int i;
	double l, r, m, val_l, val_r, val_m, dval_m, * c_, max_err;
	double min;
	
	if ((c_ = malloc(sizeof(double)*n)) == NULL)
		panic("locate_extr: can't allocate memory for c_ (%i K)",
				sizeof(double)*n/1024);
	/*
	 * a = 0.5*(1.0 - eps);
	 * b = 0.5*(1.0 + eps);
	 */
	cheb_pol_der(eps, n, c, c_);
	min = 1.0;   /* min=-1/1: search for minimum/maximum */

	if (c[n+1] < 0.0)
		min *= -1.0;
	e[0] = eps;
	e[n+1] = 1.;
	/* error function h(y) = |1 - \sqrt(y) P(y)| */
	/* maxerr = h(\epsilon) */
	max_err = fabs(1.0 - sqrt(eps)*eval_pol_cheb(eps, n, c, eps, NULL));
	/* var_r = h(1.0) */
	val_r = 1.0 - sqrt(1.0)*eval_pol_cheb(eps, n, c, 1.0, NULL);
	if (fabs(val_r) > max_err)
		max_err = fabs(val_r);
	for (i = 1; i < n+1; i++)
	{
		min *= -1.0;
		l = (i == 0 ? eps : x[i-1]);
		r = (i == n+1 ? 1.0 : x[i+1]);
		val_l = 1.0 - sqrt(l)*eval_pol_cheb(eps, n, c, l, NULL);
		val_r = 1.0 - sqrt(r)*eval_pol_cheb(eps, n, c, r, NULL);
		while (fabs(r-l) > 1.e-12)
		{
			m = 0.5*(l+r);
			val_l = 1.0 - sqrt(l)*eval_pol_cheb(eps, n, c, l, NULL);
			val_r = 1.0 - sqrt(r)*eval_pol_cheb(eps, n, c, r, NULL);
			val_m = 1.0 - sqrt(m)*eval_pol_cheb(eps, n, c, m, NULL);
			dval_m = 1.0/(2.0*sqrt(m))*eval_pol_cheb(eps, n, c, m, NULL)+
				sqrt(m)*eval_pol_cheb(eps, n-1, c_, m, NULL);
			if (dval_m*min>0)
				r = m;
			else
				l = m;
		}
		m = 0.5*(l+r);
		val_m = 1.0 - sqrt(m)*eval_pol_cheb(eps, n, c, m, NULL);
		if (fabs(val_m) > max_err)
			max_err = fabs(val_m);
		e[i]=m;
  }

  return max_err;
}/*}}}*/


/*{{{*/
/*!
 * \brief Find the polynomial of degree \f$n\f$ in the range \f$[\epsilon;1]\f$
 *
 * \param eps
 * \param deg
 * \param c		coefficients of the polynomial 
 *
 * \return double 
 */
/*}}}*/
double minmax_pol_deg(double eps, int deg, double* c)/*{{{*/
{
	int i;
	double * extr, * extr2, * extr_tmp, max_err;
	int cont;
    
	if ((extr = malloc(sizeof(double)*(deg+2))) == NULL)
			panic("minmax_pol_deg: can't allocate memory for extr (%i K)",
				sizeof(double)*(deg+2)/1024);
	if ((extr2 = malloc(sizeof(double)*(deg+2))) == NULL)
			panic("minmax_pol_deg: can't allocate memory for extr2 (%i K)",
				sizeof(double)*(deg+2)/1024);
	/* 
	 * let deg = 2
	 * extr = {0.0, 0.0, 0.0, 0.0}
	 */
    extr[0] = eps;
    extr[deg+1] = 1.0;
	/*
	 * extr = {eps, 0.0, 0.0, 1.0}
	 *              ^^^  ^^^
	 *              will be modified
	 */
    cheb_pol_extr(deg+1, extr+1);
	/*
	 * extr ={exp. Textr_2, Textr_1, 1.0}
	 */ 
	/* map extrema Textr_{n-1}...Textr_1 from interval [-1:1] to interval [eps:1 */
    for (i = 1; i < deg+1; i++)
		extr[i] = 0.5*((1.0 - eps)*extr[i] + 1.0 + eps);
	/* should be removed -- testing only */
	max_err=locate_extr(eps, deg, c, extr, extr2);
	/* */
    do
    {
		/* Stage 1 of Maehly algorithm */
		extr_pol(eps, deg, c, extr);
		/* Stage 2 of Mayehly algorithm */
		max_err = locate_extr(eps, deg, c, extr, extr2);
		extr_tmp = extr;
		extr = extr2;
		extr2 = extr_tmp;
		cont = (max_err/c[deg+1]>1.0001);
    } while (cont);
	/* free memory */
	free(extr2);
	free(extr);

	return max_err;
}/*}}}*/

/*{{{*/
/*!
 * \brief Find minmax polynomial with specified precision \f$\delta\f$ in the range
 * \f$[\epsilon; 1]\f$
 * 
 * \param eps 	\f$ \epsilon \le y \le 1\f$
 * \param delta	wanted error \f$\delta = max |h(y)|\f$
 * \param deg	degree of minmax polynomial which is determined by asymptotic
 * formala (see hep-lat/0112012) 
 * \f[
 * \delta = A\exp^{-bn\sqrt{\epsilon}},\ \ A = 0.41,\ b = 0.21
 * \f]
 * \param c		coefficients of the approximation
 * should be NULL, overwise memory will be free and reallocated.
 *
 * \return void 
 */
/*}}}*/
void minmax_pol(double eps, double delta, int * deg, double** c) /*{{{*/
{
	double err;

	(*deg) = lround(-log(delta/0.41)/(2.1*sqrt(eps)));
#ifdef MINMAX_LOG
	fprintf(stderr, "# approximated degree of poynomial: %i\n", *deg);
#endif
	/* above formula works good for deg > 30 */
	if (*deg < 5)
	    *deg = 5;
	if ((*c) != NULL)
		free(*c);
	if ((*c = malloc(sizeof(double)*((*deg)+2))) == NULL)
		panic("%s: can't allocate memory for *c (%i Kb)", __func__,
				sizeof(double)*((*deg)+2)/1024);
  	err = minmax_pol_deg(eps, (*deg), *c);
#ifdef MINMAX_LOG
	fprintf(stderr, "# err of polynomial of degree %i: %12.10f; delta = %12.10f\n", *deg, err, delta);
#endif
	/* adjust degree of polynomial */
	if (err > delta)
		/* degree is too small */
		do
		{
			(*deg)++;
			if ((*c = realloc(*c, sizeof(double)*((*deg)+2))) == NULL)
				panic("minmax_pol: can't reallocate memory for *c (%i Kb)",
					sizeof(double)*((*deg)+2)/1024);
			err = minmax_pol_deg(eps, (*deg), *c);
#ifdef MINMAX_LOG
	fprintf(stderr, "# err of polynomial of degree %i: %12.10f; delta = %12.10f\t/\\\n", *deg, err, delta);
#endif
		} while (err > delta);
	else
	{
		/* degree is too big */
		do
		{
			(*deg)--;
			if ((*c = realloc(*c, sizeof(double)*((*deg)+2))) == NULL)
				panic("minmax_pol: can't reallocate memory for *c (%i Kb)",
					sizeof(double)*((*deg)+2)/1024);
			err = minmax_pol_deg(eps, (*deg), *c);
#ifdef MINMAX_LOG
	fprintf(stderr, "# err of polynomial of degree %i: %12.10f; delta = %12.10f\t\\/\n", *deg, err, delta);
#endif
		} while (err < delta);
		(*deg)++;
		if ((*c = realloc(*c, sizeof(double)*((*deg)+2))) == NULL)
			panic("minmax_pol: can't reallocate memory for *c (%i Kb)",
				sizeof(double)*((*deg)+2)/1024);
		err = minmax_pol_deg(eps, (*deg), *c);
	}
}/*}}}*/
