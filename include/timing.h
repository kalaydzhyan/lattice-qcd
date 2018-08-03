/*{{{*/
/*!
 * \file timing.h
 *
 * \brief Timing of you code
 *
 *
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \date   Û“ƒ ÔÀ‘ 06 13:56:42 MSD 2004
 */
/*}}}*/
#ifndef _TIMING_H_
#define _TIMING_H_

#define	TIMING_MAX_ITEMS	(8)
#define TIMING_STR_LEN		(40)

#define	TIMING_WDIRAC		(0)
#define	TIMING_OVDIRAC		(1)
#define TIMING_H_WDIRAC		(2)
#define TIMING_H_WDIRAC_SQ	(3)
#define TIMING_OV_DIRAC_LN	(4)
#define TIMING_H_WDIRAC_LN	(5)
#define TIMING_A_MUL_CDSV	(6)
#define TIMING_AD_MUL_CDSV	(7)

extern long tim_ncalls[TIMING_MAX_ITEMS];
extern long long tim_timer[TIMING_MAX_ITEMS];
extern char tim_strs[TIMING_MAX_ITEMS][TIMING_STR_LEN];

/*!
 * \brief Init timing
 *
 * Set tim_ncall[i] = 0 and tim_time[i] = 0, 0 <= i < TIMING_MAX_ITEMS
 *
 * \param void
 *
 * \return void 
 */
void init_time(void);

void print_time(void);

#ifndef I686

#define start_time(_num)				\
	tim_ncalls[(_num)]++;				\
	tim_timer[(_num)] += time(NULL);

#define end_time(_num)					\
	tim_timer[(_num)] -= time(NULL);

#else

/*!
 * \brief Start time measurement
 *
 * \param ti timing item
 *
 * \return void 
 */
void start_time(int ti);

/*!
 * \brief Stom time measurement
 *
 * \param ti timings item
 *
 * \return void 
 */
void end_time(int ti);

#endif /* I686 */

#endif
