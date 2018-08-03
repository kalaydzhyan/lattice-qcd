#ifndef _MC_SU2_H_
#define _MC_SU2_H_

#include <math.h>
#include <fcntl.h>
#include <unistd.h>

#include <geometry.h>
#include <SU2-utils.h>

#define NUMBER_OF_OVERRELAX 3
#define COOLING_DELTA       0.05

extern SU2 *F;

extern int fields_init( const char *filename, uchar T, uchar L, uchar D );
extern void fields_free();
extern int read_SU2_filedata(  const char *filename );
extern int write_SU2_filedata( const char *filename );

/* returns <1/2 Tr U_p> */
extern double MC_SU2();
extern double cool_SU2();

extern double mean_plaq_SU2();
extern double get_plaq_trace(int x, int mu, int nu);

extern void   random_gauge_SU2();

#define LANDAU_GAUGE_PRECISION  0.00001
extern double Landau_gauge_SU2();

#define MAA_GAUGE_PRECISION  0.0000001
extern double MaA_gauge_SU2();

#define MAA_CENTER_GAUGE_PRECISION  0.0000001
extern double MaA_center_gauge_SU2();

#endif
