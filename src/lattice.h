/*!
 * \file lattice.h
 *
 * \brief Lattice routines
 *
 *
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \author Pavel Buividovich, email: gbuividovich@gmail.com (implemented background magnetic field, chemical potential, SU(3) gauge group in 2008 - 2009)
 */
#ifndef _LATTICE_H_
#define _LATTICE_H_

#include "defs.h"
#include "types.h"

#define FWD		0	//! < Move forward in mu direction
#define	BWD		1	//! < Move backward in mu direction
#define	FBS		2	//! < Forward_Backward_Size

/*!
 * \brief Movement on the lattice 
 */
typedef long int t_mv[VOL][DIM][FBS];	

extern t_mv * lat_mov;

#ifdef EFIELD
t_complex ev[VOL][DIM];
#endif

#define FORALL_4D(_x)											\
	for ((_x)[3] = 0; (_x)[3] < LT; (_x)[3]++)					\
	for ((_x)[2] = 0; (_x)[2] < LS; (_x)[2]++)					\
	for ((_x)[1] = 0; (_x)[1] < LS; (_x)[1]++)					\
	for ((_x)[0] = 0; (_x)[0] < LS; (_x)[0]++)

#define FORALL_3D(_x)											\
	for ((_x)[2] = 0; (_x)[2] < LT; (_x)[2]++)					\
	for ((_x)[1] = 0; (_x)[1] < LS; (_x)[1]++)					\
	for ((_x)[0] = 0; (_x)[0] < LS; (_x)[0]++)

#define FORALL_2D(_x)											\
	for ((_x)[1] = 0; (_x)[1] < LT; (_x)[1]++)					\
	for ((_x)[0] = 0; (_x)[0] < LS; (_x)[0]++)


#define INDEX_X_4D(_x)											\
	((_x)[0] + LS*(_x)[1] + LS2*(_x)[2] + LS3*(_x)[3])

#define INDEX_X_3D(_x)											\
	((_x)[0] + LS*(_x)[1] + LS2*(_x)[2])

#define INDEX_X_2D(_x)											\
	((_x)[0] + LS*(_x)[1])

#if (DIM == 4)
#define FORALL(_x) FORALL_4D((_x))
#define INDEX_X(_x) INDEX_X_4D((_x))
#elif (DIM == 3)
#define FORALL(_x) FORALL_3D((_x))
#define INDEX_X(_x) INDEX_X_3D((_x))
#elif (DIM == 2)
#define FORALL(_x) FORALL_2D((_x))
#define INDEX_X(_x) INDEX_X_2D((_x))
#endif

#define INDEX_X_MU(_x, _mu)			\
	(INDEX_X(_x)*DIM + _mu)
	

/*!
 * \brief  Move in \f$\mu\f$ direction forward one step
 *
 * \param _x  DIM-dimensional vector. Move from this site
 * \param _mu Direction of move
 *
 */
#define X_FWD(_x, _mu)													\
	((_x)[(_mu)] = (((_x)[(_mu)]+1) % lat_size[(_mu)]))

/*!
 * \brief  Move in \f$\mu\f$ direction backward one step
 *
 * \param _x  DIM-dimensional vector. Move from this site
 * \param _mu Direction of move
 *
 */
#define X_BWD(_x, _mu)													\
	((_x)[(_mu)] = 														\
	 (lat_size[(_mu)] + (_x)[(_mu)]-1) % lat_size[(_mu)])

/*****************************************************************************/

t_gauge_vector * lat_gauge_create(void);
void lat_gauge_destroy(t_gauge_vector * gv);
void lat_gauge_load(t_gauge_vector * gv, char * info_fname);
void lat_gauge_check_unit(t_gauge_vector * gv);
t_real lat_gauge_action(t_gauge_vector * gv);
void lat_mov_init(void);
void lat_gauge_identity(t_gauge_vector * gv);

void   xcoord(int x, int* xc);
int    taxi_distance(int*, int*);
t_real distance4D(int x1, int x2);

void lat_gauge_instanton(t_gauge_vector* gv, t_real rho, t_real* x0);
void lat_gauge_aperiodic(t_gauge_vector* gv);
void init_efield(int E, int H);
#endif
