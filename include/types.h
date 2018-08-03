/*{{{*/
/*!
 * \file types.h
 *
 * \brief Types definition 
 *
 *
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \author Pavel Buividovich, email: gbuividovich@gmail.com (implemented background magnetic field, chemical potential, SU(3) gauge group in 2008 - 2009)
 *//*}}}*/
#ifndef _TYPES_H_
#define _TYPES_H_

#include <complex.h>

#include "defs.h"

#ifdef LONG_INTERNAL_ARITHMETIC
typedef long double         t_double_real;    // real type for high-precision internal arithmetics
typedef long double complex t_double_complex; // complex type for high-precision internal arithmetics
#else
typedef double              t_double_real;    // real type for high-precision internal arithmetics
typedef double complex      t_double_complex; // complex type for high-precision internal arithmetics
#endif

#ifdef SINGLE_PREC
typedef float t_real;	
typedef float complex t_complex;

/*! \brief Wrapper to complex.h function (float) */
#define CCALL(callname) (callname ## f)	

#else
/*! \brief Real number of double precision */
typedef double t_real;
/*! \brief Complex number of double precision  */
typedef double complex t_complex;

/*! \brief Wrapper to complex.h function (double) */
#define CCALL(callname) (callname)	

#endif

/*!
 * \brief Wrapper for all matricies routines
 *	Add gauge group suffix to the name of function
 * \param callname	Function name
 */
 
#ifdef SU2
#define NCOLORS 2
#define MATCALL(callname)	callname ## _SU2
#define N_GMAT_COMP	2	//!< Number of elements in gauge matrix to be stored
#endif

#ifdef SU3
#define NCOLORS 3
#define MATCALL(callname) 	callname ## _SU3
#define N_GMAT_COMP	9	//!< Number of elements in gauge matrix to be stored
#endif

/*! \brief SU2 Gauge MATrix
 *	
 *	Structure of the SU2 matrix \f$ U \f$:
 *  \f[ U = U_0 + i \sum_{i=1}^3 \sigma_i U_i =
 *	\left(\begin{tabular}{cc} 
 *  $U_{11}$ & $U_{12}$  \\
 *  $U_{21}$ & $U_{22}$  \\
 *  \end{tabular} \right)  = 
 *	\left(\begin{tabular}{cc} 
 *  $U_{11}$    & $U_{12}$		\\
 *  $-U_{12}^*$ & $U_{11}^*$	\\
 *  \end{tabular} \right) \f]
 * 	\f[ U_{11} = U_0 + i U_3 \f]
 * 	\f[ U_{12} = U_2 + i U_1 \f]
 *	
 *	We choose to store in memory only 2 elements \f$ U_{11} \f$ and \f$ U_{12}\f$
 *	because others can be easely reconstructed.
 */
typedef t_complex t_gmat[N_GMAT_COMP];

extern t_gmat ISIGMA1_SU2;
extern t_gmat ISIGMA2_SU2;
extern t_gmat ISIGMA3_SU2;

/*! 
 * \brief ColorDiracSpace vector
 */
typedef t_complex t_cds_vector[VOL][NCOLORS][NDIRAC];	
/*! 
 * \brief  Gauge field configuration
 */
typedef t_gmat t_gauge_vector[VOL][DIM];
typedef t_real t_extfield[VOL][DIM];
/*!
 * \brief data for minmax version of overlap operator
 */
typedef struct t_ov_data_mm {
	t_real kappa;		//!< wilson \f$\kappa = \frac{1}{2a\rho+8r}\f$
	t_real mass;		//!< Overlap mass. \f$ mass = \frac{am_q}{2} \f$
	int deg;			//!< degree of minmax polynomial for zero mu, or the size of Krylov subspace for nonzero mu
	int n_proj;			//!< number of projectors (the number of deflated critical eigenvalues)
    t_double_real sign_tol;	// tolerance in the calculation of the sign function of the inner matrix
	t_complex* evals;	//!< eigenvalues of Q
#ifdef MU
    t_complex* revecs;
    t_complex* levecs;
#else
	t_complex* evecs;	//!< eigenvectors of Q
	double* coef;		//!< coefficients of minmax polynomial
	t_real l_sq_min;	//!<\f$\lambda_{min}^2\f$
	t_real l_sq_max;	//!<\f$\lambda_{max}^2\f$
#endif
} t_ov_data_mm;

typedef void (*t_op)(void*, t_gauge_vector*, t_cds_vector*, t_cds_vector*);

typedef struct t_op_ev {
	int 		ls;
	int 		lt;
	int 		dim;
	int			ncolors;
	int			ndirac;
	t_op		op;
	int			nev;	//!< how many eigenvalues to compute. \f$ 0 < pnev < N-1 \f$
	/*! \brief which eigenvalues:
	 *
 	 * - 'LM' - want the nev eigenvalues of largest magnitude 
	 * - 'SM' - want the nev eigenvalues of smallest magnutide 
	 * - 'LR' - want the nev eigenvalues of largest real part
	 * - 'SR' - want the nev eigenvalues of smallest real part 
	 * - 'LI' - want the nev eigenvalues of largest imaginary part 
	 * - 'SI' - want the nev eigenvalues of smallest imaginary part 
	 */
	char		which[2];	
	int			val_vec;	//!<  = 0 - only values, = 1 - values & vectors
	/*! \brief 	Tolerance. 
	 * 
	 * \f$ |\lambda_c - \lambda_t| < ptol |\lambda_c|\f$, where 
	 * \f$\lambda_c\f$ - computed eigenvalue,
	 * \f$\lambda_t\f$ - eigenvalue of matrix nearest to \f$\lambda_c\f$
	 */
	t_real		tol;
	t_complex*	evals;		//!< on return: array of eigenvalues
	t_complex*	evecs;		//!< on return: should be NULL TODO
	/*! \brief t_ov_data_mm for Overlap or kappa for Wilson */
	void*		data;
	int			data_type;	//!< type of data (float or double) (used only in read)
} t_op_ev;

#define EVD_DATA_FLOAT	(4)
#define EVD_DATA_DOUBLE	(8)

#endif
