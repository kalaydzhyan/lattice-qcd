/*!
 * \file defs.h
 *
 * \brief Various definitions and macros
 *
 *
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \author Pavel Buividovich, email: gbuividovich@gmail.com (implemented background magnetic field, chemical potential, SU(3) gauge group in 2008 - 2009)
 */

#ifndef _DEFS_H_
#define _DEFS_H_

#define	DIM			4						//!< Lattice dimension
#define LS			(LAT_S)					//!< Lattice space size 
#define LT 			(LAT_T)					//!< Lattice time size
#define	LS2			(LS*LS)					//!< LS*LS to speedup index calculation
#define LS3			(LS*LS*LS)				//!< to speedup index calculation in 3D

#if(DIM == 4)
#define N_SITES 	(LS*LS*LS*LT)			//!< Number of lattice sites (DIM = 4)
#elif(DIM == 3)
#define N_SITES		(LS*LS*LT)				//!< Number of lattice sites (DIM = 3)
#elif(DIM == 2)
#define N_SITES		(LS*LT)					//!< Number of lattice sites (DIM = 2)
#endif
#define N_LINKS 	(N_SITES*DIM)			//!< Number of lattice links
#define N_PLAQS 	(N_SITES*DIM*(DIM-1)/2)	//!< Number of plaquettes

#define	VOL			N_SITES					//!< Lattice volume = number of sites
#define	VOLH		(VOL/2)					//!< half lattice volume

#ifdef SCALAR
#define NDIRAC		(1)						//!< number of dirac indicies
#else
#define NDIRAC		(4)						//!< number of dirac indicies
#endif

/*!
 * \brief Wrapper for fortran functions
 *
 *	Fortran function without underscore in its name should be called using
 *	this macro. This is because fortran compiler adds one underscore to the
 *	name of function.
 * 
 * \param a Fortran function name 
 */
#define UNDERSC(a) a ## _
/*!
 * \brief Wrapper for fortran functions
 *
 *	Fortran function with underscore in its name should be called using
 *	this macro. This is because fotran compiler adds two underscore to the
 *	name of function already containing one or more underscore in its name.
 * 
 * \param a Fortran function name 
 */

#ifdef  UNSC
#define UNDERSC2(a) a ## _
#elif defined UNSC2
#define UNDERSC2(a) a ## __
#endif

#endif

