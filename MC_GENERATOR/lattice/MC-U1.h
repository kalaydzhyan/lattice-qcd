#ifndef _MC_U1_H_
#define _MC_U1_H_

#include <math.h>
#include <fcntl.h>
#include <unistd.h>
#include <complex.h>
#include <geometry.h>

#define NUMBER_OF_OVERRELAX_U1 3

extern double complex *F_U1;
extern int *mono;


extern int  fields_init_U1( const char *filename, uchar T, uchar L, uchar D );
extern void fields_free_U1();
extern int  read_filedata_U1(  const char *filename );
extern int  write_filedata_U1( const char *filename );

/* returns <1/2 Tr U_p> */
extern double MC_U1();
extern double mean_plaq_U1();
extern void   random_gauge_U1();
extern double Landau_gauge_U1();

extern int * get_mono_U1();
extern double mono_density_U1();
#endif
