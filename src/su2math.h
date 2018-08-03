/*!
 * \file su2math.h
 *
 * \brief SU2 macroses
 *
 *
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \date   Чтв Сен 23 12:22:04 MSD 2004
 */
#ifndef _SU2MATH_H_
#define _SU2MATH_H_

#include <math.h>
#include "defs.h"
#include "types.h"

/*!
 * \brief Initialize SU2 matrix _m with zero matrix
 *
 * \param _m t_gmat
 *
 */
#define		LOADZERO_SU2(_m)	\
	(_m)[0] = 0.0 + I*0.0;		\
	(_m)[1] = 0.0 + I*0.0;

/*!
 * \brief Initialize SU2 matrix _m with UNITY matrix
 *
 * \param _m t_gmat
 *
 */
#define		LOADIDENTITY_SU2(_m)\
	(_m)[0] = 1.0 + I*0.0;		\
	(_m)[1] = 0.0 + I*0.0;

/*!
 * \brief Load vector _v to matrix _m
 *
 * \param _m	t_gmat
 * \param _v 	t_real[4] = {u0, u1, u2, u3}
 *
 */
#define		LOADVECTOR_SU2(_m, _v)	\
	(_m)[0] = (_v)[0] + I*(_v[3]);	\
	(_m)[1] = (_v)[2] + I*(_v[1]);										

/*!
 * \brief Load complex constant _c to matrix _m
 *
 * \param _m t_gmat
 * \param _c t_complex
 *
 */
#define		LOADCONST_SU2(_m, _c)	\
	(_m)[0] = (_c);					\
	(_m)[1] = (_c);

/*!
 * \brief M2 = M2 + M1
 *
 * \param _m2 t_gmat
 * \param _m1 t_gmat
 *
 */
#define 	A_add_B_SU2(_m2, _m1)	\
	(_m2)[0] += (_m1)[0];			\
	(_m2)[1] += (_m1)[1];												

/*!
 * \brief  M = a*M
 *
 * \param _a t_complex
 * \param _m t_gmat
 *
 */
#define 	a_mul_A_SU2(_a, _m)		\
	(_m)[0] *= (_a);				\
	(_m)[1] *= (_a);

/*!
 * \brief M3 = M1*M2
 *
 * \param _m1 t_gmat
 * \param _m2 t_gmat
 * \param _m3 t_gmat
 *
 */
#define 	C_eq_A_mul_B_SU2(_m3, _m1, _m2)							\
	(_m3)[0] = (_m1)[0]*(_m2)[0] - (_m1)[1]*CCALL(conj)((_m2)[1]);	\
	(_m3)[1] = (_m1)[0]*(_m2)[1] + (_m1)[1]*CCALL(conj)((_m2)[0]);

/*!
 * \brief M3 = M1+M2
 *
 * \param _m1 t_gmat
 * \param _m2 t_gmat
 * \param _m3 t_gmat
 *
 */
#define 	C_eq_A_add_B_SU2(_m3, _m1, _m2)	\
	(_m3)[0] = (_m1)[0]+(_m2)[0];			\
	(_m3)[1] = (_m1)[1]+(_m2)[1];			\

/*!
 * \brief M3 += M1*M2
 *
 * \param _m1 t_gmat
 * \param _m2 t_gmat
 * \param _m3 t_gmat
 *
 */
#define 	C_peq_A_mul_B_SU2(_m3, _m1, _m2)						\
	(_m3)[0] += (_m1)[0]*(_m2)[0] - (_m1)[1]*CCALL(conj)((_m2)[1]);	\
	(_m3)[1] +=	(_m1)[0]*(_m2)[1] + (_m1)[1]*CCALL(conj)((_m2)[0]);

/*!
 * \brief \f$ M3 = M1^{\dagger}*M2 \f$
 *
 * \param _m1 t_gmat
 * \param _m2 t_gmat
 * \param _m3 t_gmat
 *
 */
#define 	C_eq_Ad_mul_B_SU2(_m3, _m1, _m2)								 \
	(_m3)[0] = CCALL(conj)((_m1)[0])*(_m2)[0]+(_m1)[1]*CCALL(conj)((_m2)[1]);\
	(_m3)[1] = CCALL(conj)((_m1)[0])*(_m2)[1]-(_m1)[1]*CCALL(conj)((_m2)[0]);

/*!
 * \brief \f$ M3 = M1*M2^{\dagger} \f$
 *
 * \param _m1 t_gmat
 * \param _m2 t_gmat
 * \param _m3 t_gmat
 *
 */
#define 	C_eq_A_mul_Bd_SU2(_m3, _m1, _m2)								 \
	(_m3)[0] = (_m1)[0]*CCALL(conj)((_m2)[0])+(_m1)[1]*CCALL(conj)((_m2)[1]);\
	(_m3)[1] =-(_m1)[0]*(_m2)[1]+(_m1)[1]*(_m2)[0];


/*!
 * \brief  \f$ M3 = M1^{\dagger}*M2^{\dagger}\f$
 *
 * \param _m1 t_gmat
 * \param _m2 t_gmat
 * \param _m3 t_gmat
 *
 */
#define 	C_eq_Ad_mul_Bd_SU2(_m3, _m1, _m2)								\
	(_m3)[0] = CCALL(conj)((_m1)[0])*CCALL(conj)((_m2)[0])-					\
				            (_m1)[1]*CCALL(conj)((_m2)[1]);					\
	(_m3)[1] =-CCALL(conj)((_m1)[0])*(_m2)[1] - (_m1)[1]*(_m2)[0]; 

/*!
 * \brief Re(Tr(M1*M2)
 *
 * \param _m1 t_gmat
 * \param _m2 t_gmat
 *
 */
#define 	re_tr_A_mul_B_SU2(_m1, _m2)						\
	(2.0*(CCALL(creal)((_m1)[0])*CCALL(creal)((_m2)[0]) -	\
	CCALL(cimag)((_m1)[0])*CCALL(cimag)((_m2)[0]) -			\
	CCALL(creal)((_m1)[1])*CCALL(creal)((_m2)[1]) -			\
	CCALL(cimag)((_m1)[1])*CCALL(cimag)((_m2)[1]))) 

/*!
 * \brief Tr(M)
 *
 * \param _m t_gmat
 *
 */
#define 	tr_A_SU2(_m)									\
	(2.0*CCALL(creal)((_m)[0]))

/*!
 * \brief normalize(M)
 *
 * \param _m
 *
 */
#define 	norm_A_SU2(_m)									\
{															\
	t_real tmp = sqrt(CCALL(creal)(det_A_SU2((_m))));		\
	(_m)[0] /= tmp;											\
	(_m)[1] /= tmp;											\
}

/*!
 * \brief Hconj(M) 
 *
 * \param _m
 *
 */
#define 	hconj_A_SU2(_m)									\
	(_m)[0] = CCALL(conj)((_m)[0]);							\
	(_m)[1] = -(_m)[1];														

/*!
 * \brief M = inv(M) = M^{-1}
 *
 *	If
 *	\f[ A = \left(\begin{tabular}{cc}a & b \\ c & d \\ \end{tabular}\right)\f]
 *	then
 *	\f[ A^{-1} =\frac{1}{\det A}
 *	\left(\begin{tabular}{cc}d & -b \\ -c & a \\ \end{tabular}\right)\f]
 *	For \f$ SU(2) \f$ case
 *	\f[ U^{-1} = \frac{1}{\det U}
 *	\left(\begin{tabular}{cc}$U_{11}^*$ & $-U_{12}$ \\ $U_{12}^*$ & $U_{11}$\\
 *	\end{tabular}\right)\f]
 * 
 * \param _m
 *
 */
#define 	inv_A_SU2(_m)								\
{														\
	t_real	idet = 1.0/CCALL(creal)(det_A_SU2(_m));		\
	(_m)[0] = CCALL(conj)((_m)[0])*idet;				\
	(_m)[1] = -((_m)[1])*idet;							\
}

/*!
 * \brief Generate random SU2 matrix M
 *
 * \param _m t_gmat
 *
 */
#ifdef SINGLE_PREC
#define		rand_A_SU2(_m)									\
{															\
	t_real o1, o2, o3, tmp, tmp2;							\
	t_real or = M_PI;										\
	t_real om = M_PI/2.0;									\
	o1 = genrand_real2();									\
	o2 = genrand_real2();									\
	o3 = genrand_real2(); 									\
	o1 = or*o1 - om; o2 = or*o2 - om; o3 = or*o3 - om;		\
	tmp = sqrtf(o1*o1 + o2*o2 + o3*o3);						\
	if (tmp != 0.0)											\
		tmp2 = sinf(tmp)/tmp;								\
	else													\
		tmp2 = 0.0;											\
	(_m)[0] = cosf(tmp) + I*o3*tmp2;						\
	(_m)[1] = o2*tmp2  + I*o1*tmp2;		 					\
}																		
#else
#define		rand_A_SU2(_m)									\
{															\
	t_real o1, o2, o3, tmp, tmp2;							\
	t_real or = M_PI;										\
	t_real om = M_PI/2.0;									\
	o1 = genrand_res53();									\
	o2 = genrand_res53();									\
	o3 = genrand_res53(); 									\
	o1 = or*o1 - om; o2 = or*o2 - om; o3 = or*o3 - om;		\
	tmp = sqrt(o1*o1 + o2*o2 + o3*o3);						\
	if (tmp != 0.0)											\
		tmp2 = sin(tmp)/tmp;								\
	else													\
		tmp2 = 0.0;											\
	(_m)[0] = cos(tmp) + I*o3*tmp2;							\
	(_m)[1] = o2*tmp2  + I*o1*tmp2;		 					\
}															
#endif

/*!
 * \brief Det(M)
 *
 * \param _m t_gmat
 *
 */
#define		det_A_SU2(_m)									\
	((_m)[0]*CCALL(conj)((_m)[0]) + (_m)[1]*CCALL(conj)((_m)[1]))

/*****************************************************************************/

void A_mul_CDSV_SU2(t_complex* out, t_gmat A, t_complex* in);
void Ad_mul_CDSV_SU2(t_complex* out, t_gmat A, t_complex* in);

void print_mat_SU2(t_complex * m);

#endif
