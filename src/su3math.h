/*!
 * \file su3math.h
 *
 * \brief SU3 macroses
 *
 *
 * $Id$
 *
 * \author Pavel Buividovich, email: gbuividovich@gmail.com (implemented background magnetic field, chemical potential, SU(3) gauge group in 2008 - 2009)
 */
#ifndef _SU3MATH_H_
#define _SU3MATH_H_

#include <math.h>
#include "defs.h"
#include "types.h"

/*!
 * \brief Initialize SU3 matrix _m with zero matrix
 *
 * \param _m t_gmat
 *
 */
 
#define		LOADZERO_SU3(_m)	\
   {(_m)[0] = 0.0 + I*0.0;		\
    (_m)[1] = 0.0 + I*0.0;		\
    (_m)[2] = 0.0 + I*0.0;		\
    (_m)[3] = 0.0 + I*0.0;		\
    (_m)[4] = 0.0 + I*0.0;		\
    (_m)[5] = 0.0 + I*0.0;		\
    (_m)[6] = 0.0 + I*0.0;		\
    (_m)[7] = 0.0 + I*0.0;		\
    (_m)[8] = 0.0 + I*0.0;};	\

/*!
 * \brief Initialize SU3 matrix _m with UNITY matrix
 *
 * \param _m t_gmat
 *
 */
 
#define		LOADIDENTITY_SU3(_m)    \
   {(_m)[0] = 1.0 + I*0.0;		\
    (_m)[1] = 0.0 + I*0.0;		\
    (_m)[2] = 0.0 + I*0.0;		\
    (_m)[3] = 0.0 + I*0.0;		\
    (_m)[4] = 1.0 + I*0.0;		\
    (_m)[5] = 0.0 + I*0.0;		\
    (_m)[6] = 0.0 + I*0.0;		\
    (_m)[7] = 0.0 + I*0.0;		\
    (_m)[8] = 1.0 + I*0.0;}		


/*!
 * \brief Load complex constant _c to matrix _m
 *
 * \param _m t_gmat
 * \param _c t_complex
 *
 */
#define		LOADCONST_SU3(_m, _c)			\
   {(_m)[0] = (_c);					\
	(_m)[1] = (_c);					\
	(_m)[2] = (_c);					\
    (_m)[3] = (_c);					\
	(_m)[4] = (_c);					\
	(_m)[5] = (_c);					\
	(_m)[6] = (_c);					\
	(_m)[7] = (_c);					\
	(_m)[8] = (_c);}					
/*!
 * \brief M2 = M2 + M1
 *
 * \param _m2 t_gmat
 * \param _m1 t_gmat
 *
 */
#define 	A_add_B_SU3(_m2, _m1)		\
   {(_m2)[0] += (_m1)[0];			\
	(_m2)[1] += (_m1)[1];			\
	(_m2)[2] += (_m1)[2];			\
	(_m2)[3] += (_m1)[3];			\
	(_m2)[4] += (_m1)[4];			\
	(_m2)[5] += (_m1)[5];			\
	(_m2)[6] += (_m1)[6];			\
	(_m2)[7] += (_m1)[7];			\
	(_m2)[8] += (_m1)[8];}			
												

/*!
 * \brief  M = a*M
 *
 * \param _a t_complex
 * \param _m t_gmat
 *
 */
#define 	a_mul_A_SU3(_a, _m)			\
	{(_m)[0] *= (_a);				\
	(_m)[1] *= (_a);				\
	(_m)[2] *= (_a);				\
	(_m)[3] *= (_a);				\
	(_m)[4] *= (_a);				\
	(_m)[5] *= (_a);				\
	(_m)[6] *= (_a);				\
	(_m)[7] *= (_a);				\
	(_m)[8] *= (_a);}				

/*!
 * \brief M3 = M1*M2
 *
 * \param _m1 t_gmat
 * \param _m2 t_gmat
 * \param _m3 t_gmat
 *
 */
#define 	C_eq_A_mul_B_SU3(_m3, _m1, _m2)					 \
    for(tmp_i=0; tmp_i<3; tmp_i++)                                               \
     for(tmp_j=0; tmp_j<3; tmp_j++)                                              \
     {                                                                           \
      (_m3)[3*tmp_i + tmp_j] = 0.0 + I*0.0;                                      \
      for(tmp_k=0; tmp_k<3; tmp_k++)                                             \
       (_m3)[3*tmp_i + tmp_j] += (_m1)[3*tmp_i + tmp_k]*(_m2)[3*tmp_k + tmp_j];  \
     };                                                                          


/*!
 * \brief M3 = M1+M2
 *
 * \param _m1 t_gmat
 * \param _m2 t_gmat
 * \param _m3 t_gmat
 *
 */
#define 	C_eq_A_add_B_SU3(_m3, _m1, _m2)	\
   {(_m3)[0] = (_m1)[0] + (_m2)[0];			\
	(_m3)[1] = (_m1)[1] + (_m2)[1];			\
	(_m3)[2] = (_m1)[2] + (_m2)[2];			\
    (_m3)[3] = (_m1)[3] + (_m2)[3];			\
	(_m3)[4] = (_m1)[4] + (_m2)[4];			\
	(_m3)[5] = (_m1)[5] + (_m2)[5];			\
	(_m3)[6] = (_m1)[6] + (_m2)[6];			\
	(_m3)[7] = (_m1)[7] + (_m2)[7];			\
	(_m3)[8] = (_m1)[8] + (_m2)[8];}
/*!
 * \brief M3 += M1*M2
 *
 * \param _m1 t_gmat
 * \param _m2 t_gmat
 * \param _m3 t_gmat
 *
 */
#define 	C_peq_A_mul_B_SU3(_m3, _m1, _m2)				 \
    for(tmp_i=0; tmp_i<3; tmp_i++)                                               \
     for(tmp_j=0; tmp_j<3; tmp_j++)                                              \
     {                                                                           \
      for(tmp_k=0; tmp_k<3; tmp_k++)                                             \
       (_m3)[3*tmp_i + tmp_j] += (_m1)[3*tmp_i + tmp_k]*(_m2)[3*tmp_k + tmp_j];  \
     };

/*!
 * \brief \f$ M3 = M1^{\dagger}*M2 \f$
 *
 * \param _m1 t_gmat
 * \param _m2 t_gmat
 * \param _m3 t_gmat
 *
 */
#define 	C_eq_Ad_mul_B_SU3(_m3, _m1, _m2)				 \
    for(tmp_i=0; tmp_i<3; tmp_i++)                                               \
     for(tmp_j=0; tmp_j<3; tmp_j++)                                              \
     {                                                                           \
      (_m3)[3*tmp_i + tmp_j] = 0.0 + I*0.0;                                      \
      for(tmp_k=0; tmp_k<3; tmp_k++)                                             \
       (_m3)[3*tmp_i + tmp_j] += CCALL(conj)((_m1)[3*tmp_k + tmp_i])*(_m2)[3*tmp_k + tmp_j];  \
     }; 

/*!
 * \brief \f$ M3 = M1*M2^{\dagger} \f$
 *
 * \param _m1 t_gmat
 * \param _m2 t_gmat
 * \param _m3 t_gmat
 *
 */
#define 	C_eq_A_mul_Bd_SU3(_m3, _m1, _m2)				 \
    for(tmp_i=0; tmp_i<3; tmp_i++)                                               \
     for(tmp_j=0; tmp_j<3; tmp_j++)                                              \
     {                                                                           \
      (_m3)[3*tmp_i + tmp_j] = 0.0 + I*0.0;                                      \
      for(tmp_k=0; tmp_k<3; tmp_k++)                                             \
       (_m3)[3*tmp_i + tmp_j] += (_m1)[3*tmp_i + tmp_k]*CCALL(conj)((_m2)[3*tmp_j + tmp_k]);  \
     };  
     

/*!
 * \brief  \f$ M3 = M1^{\dagger}*M2^{\dagger}\f$
 *
 * \param _m1 t_gmat
 * \param _m2 t_gmat
 * \param _m3 t_gmat
 *
 */
#define 	C_eq_Ad_mul_Bd_SU3(_m3, _m1, _m2)				 \
    for(tmp_i=0; tmp_i<3; tmp_i++)                                               \
     for(tmp_j=0; tmp_j<3; tmp_j++)                                              \
     {                                                                           \
      (_m3)[3*tmp_i + tmp_j] = 0.0 + I*0.0;                                      \
      for(tmp_k=0; tmp_k<3; tmp_k++)                                             \
       (_m3)[3*tmp_i + tmp_j] += CCALL(conj)((_m1)[3*tmp_k + tmp_i])*CCALL(conj)((_m2)[3*tmp_j + tmp_k]);  \
     };  

/*!
 * \brief Re(Tr(M1*M2)
 *
 * \param _m1 t_gmat
 * \param _m2 t_gmat
 *
 */
#define 	re_tr_A_mul_B_SU3(_m1, _m2)		               	\
	CCALL(creal)((_m1)[0]*(_m2)[0] + (_m1)[1]*(_m2)[3] +            \
    (_m1)[2]*(_m2)[6] + (_m1)[3]*(_m2)[1] + (_m1)[4]*(_m2)[4] +         \
    (_m1)[5]*(_m2)[7] + (_m1)[6]*(_m2)[2] + (_m1)[7]*(_m2)[5] +         \
    (_m1)[8]*(_m2)[8]) 

/*!
 * \brief Re Tr(M)
 *
 * \param _m t_gmat
 *
 */
#define 	tr_A_SU3(_m)						\
	CCALL(creal)((_m)[0] + (_m)[4] + (_m)[8])

/*!
 * \brief Hconj(M) 
 *
 * \param _m
 *
 */
#define 	hconj_A_SU3(_m)					\
{   t_complex hcajtmp[9];                                       \
	hcajtmp[0] = CCALL(conj)((_m)[0]);			\
	hcajtmp[1] = CCALL(conj)((_m)[3]);			\
	hcajtmp[2] = CCALL(conj)((_m)[6]);			\
	hcajtmp[3] = CCALL(conj)((_m)[1]);			\
	hcajtmp[4] = CCALL(conj)((_m)[4]);			\
	hcajtmp[5] = CCALL(conj)((_m)[7]);			\
	hcajtmp[6] = CCALL(conj)((_m)[2]);			\
	hcajtmp[7] = CCALL(conj)((_m)[5]);			\
	hcajtmp[8] = CCALL(conj)((_m)[8]);			\
	(_m)[0] = hcajtmp[0];                                   \
	(_m)[1] = hcajtmp[1];                                   \
	(_m)[2] = hcajtmp[2];                                   \
	(_m)[3] = hcajtmp[3];                                   \
	(_m)[4] = hcajtmp[4];                                   \
	(_m)[5] = hcajtmp[5];                                   \
	(_m)[6] = hcajtmp[6];                                   \
	(_m)[7] = hcajtmp[7];                                   \
	(_m)[8] = hcajtmp[8];}
	
/*!
 * \brief Det(M)
 *
 * \param _m t_gmat
 *
 */
#define		det_A_SU3(_m)									\
	(-(_m)[2]*(_m)[4]*(_m)[6] + (_m)[1]*(_m)[5]*(_m)[6] + (_m)[2]*(_m)[3]*(_m)[7] - (_m)[0]*(_m)[5]*(_m)[7] - (_m)[1]*(_m)[3]*(_m)[8] + (_m)[0]*(_m)[4]*(_m)[8])

/*****************************************************************************/

void A_mul_CDSV_SU3(t_complex* out, t_gmat A, t_complex* in);
void Ad_mul_CDSV_SU3(t_complex* out, t_gmat A, t_complex* in);

void print_mat_SU3(t_complex * m);

#endif
