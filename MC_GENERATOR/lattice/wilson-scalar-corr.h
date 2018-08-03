#ifndef _WILSON_SCALAR_CORR_H_
#define _WILSON_SCALAR_CORR_H_

#include <MC-SU2.h>
#include <fftw.h>
#include <wilson-loops.h>

/* Possible values:
   0 - uses lot of memory, fast;
   1 - uses less memory, slower;
   2 - minimal memory, slowest */
extern uchar wilson_scalar_lowmem;

/* Wilson loops <-> scalar correlation function.
   Arguments:
   smearing - whether to use smeared Wilson loops, smearing parameters are set
              in willson-loops.c;
   new_run  - if !=0 then all final data files are truncated first; otherwise
              data files will be appended;
   callback - array of callbacks which determine scalar observables to be measured;
              argument - index of point at which the scalar value is requested;
   filename - array of data files names, one per scalar observable;
              routine dumps obtained data to these files;
   N_callback - total number of scalars */
extern void make_wilson_scalar( uchar smearing, uchar new_run,
				double (*callback[])( uint ), char *filename[],
				uchar N_callback );
extern void free_wilson_scalar();
extern void mean_wilson_scalar( char *old, char *new );
#endif
