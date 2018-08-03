#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>
#include <float.h>
#include <strings.h>

#define MAX_BUF  1024

#define IMPROVEMENT

typedef unsigned char uchar;
#ifndef RAVEN
typedef unsigned int  uint;
#endif
//#define SYSTEM_RANDOM
#ifdef SYSTEM_RANDOM
extern void   rnd_init(void);
#define RND()    (((double)random())/((double)RAND_MAX+1.0))
#else
#include <random.h>
#define RND()    (rnd_get())
#endif

typedef struct{
  uchar D, *size;
  uint *ipw;
  double beta;
#ifdef IMPROVEMENT
//uzero is just mean plaquette i.e. u0^4;
  double uzero,beta1, beta2, beta3;
#endif
  double gamma, lambda; /* for nAHM only */
  uchar q;              /* for nAHM only */
} Param;
extern Param *param;

/* ----------------------------------------------------------------- */
extern int  param_init( uchar T, uchar L, uchar D );
#ifdef IMPROVEMENT
extern void param_improvement();
#endif
extern void param_free();
extern uint site_index( uchar *x );
extern void site_coordinates( uchar *x, uint m );
extern uint index_up( int m, uchar k );
extern uint index_down( int m, uchar k );
extern uint get_distance2( uint m1, uint m2 );
extern double mod_pi( double x );
extern double mod_one_half( double x );
extern void   abuse( char *module, int line, char *fmt, ... );
extern int parse_geometry( char *pattern, uchar *T, uchar *L, uchar *D );

/* Number of planes D*(D-1)/2, enumeration:
   plane_enum(i,j) + (D*(D-1)/2) * site_index
   0 <=i,j < D */
extern uchar plane_enum( uchar i, uchar j );
extern uint  plaq_index( uchar i, uchar j, uint m );
extern uint  link_index( uchar i, uint m );
extern void  D3dual( uchar i0, uchar *i1, uchar *i2 );
extern void  D4dual( uchar i0, uchar i1, uchar *i2, uchar *i3 );

#ifdef DEBUG
extern void print_link(uint m, uchar i0 );
#endif
/* ----------------------------------------------------------------- */
#define MEMORY_CHECK(x,y)  if( !(x) ){ error("Memory allocation failed"); return y; };

#define error(x)           abuse(__FILE__,__LINE__,x)
#define error2(x1,x2)         abuse(__FILE__,__LINE__,x1,x2)
#define error3(x1,x2,x3)         abuse(__FILE__,__LINE__,x1,x2,x3)
#define error4(x1,x2,x3,x4)         abuse(__FILE__,__LINE__,x1,x2,x3,x4)
#define error5(x1,x2,x3,x4,x5)         abuse(__FILE__,__LINE__,x1,x2,x3,x4,x5)

#define panic(x)           { abuse(__FILE__,__LINE__,x); exit(EXIT_FAILURE); }
#define panic2(x1,x2)         { abuse(__FILE__,__LINE__,x1,x2); exit(EXIT_FAILURE); }
#define panic3(x1,x2,x3)         { abuse(__FILE__,__LINE__,x1,x2,x3); exit(EXIT_FAILURE); }
#define panic4(x1,x2,x3,x4)         { abuse(__FILE__,__LINE__,x1,x2,x3,x4); exit(EXIT_FAILURE); }
#define panic5(x1,x2,x3,x4,x5)         { abuse(__FILE__,__LINE__,x1,x2,x3,x4,x5); exit(EXIT_FAILURE); }
/* ----------------------------------------------------------------- */

/* Generator of gaussian random numbers distributed with probability
     exp{ -\alpha ( x - x_0 )^2   */
extern double gauss_random( double alpha, double x0 );


/* Generator of POSITIVE gaussian random numbers distributed with probability
     exp{ -\alpha ( x - x_0 )^2   */
extern double positive_gauss_random( double alpha, double x0 );


/* Generator of random numbers distributed with probability
       exp[ -alpha (x - x0)^2 - gamma (1 - cos( q x ) ) ]
*/
extern double heatbath_nAHM( double alpha, double x0, double gamma, uchar q );


/* Generator of positive random numbers distributed with probability
   exp[ -lambda (x - 1)^2 - alpha ( sqrt(x) - x0 )^2 -sqrt(x) beta ]  */
extern double heatbath_higgs_nAHM( double lambda, double alpha, double x0, double beta );

/* ----------------------------------------------------------------- */
/*                         SU2 stuff.                                */
/* ----------------------------------------------------------------- */
/*    This defines "critical alpha" for Wilson model: we use Kennedy algorithm
      for alpha>ALPHA_CROSS and standard heatbath overwise */
#define ALPHA_CROSS    1.685

/*  Generates random number x\in[-1;1] with probability \sqrt{1-x^2}exp[alpha*x]
    Creutz algorithm, which works fast only for not so large \alpha (~< 1.685 )
    P(x) ~ \sqrt{1-x^2} e^{\alpha x} dx ~ \sqrt{ 1- (\frac{ \ln z }{ \alpha })^2 } dz
    x \in [-1;1],   z \in [ e^{-\alpha}; e^{\alpha} ]
    We are generating z and then return \ln z / \alpha:
       x = rnd() -- \in [0;1]
       z = e^{-\alpha} + x ( e^{\alpha} - e^{-\alpha} ) = e^{\alpha} ( x + e^{-2 \alpha}( 1-x ) )
       help = \ln z / \alpha = 1.0 + \ln[x + e^{-2 \alpha}( 1-x )] / \alpha
*/
extern double Creutz( double alpha );

/*  Generates random number x\in[-1;1] with probability \sqrt{1-x^2}exp[alpha*x]
    Algorithm of Kennedy-Pendleton (PLB 156,393) sutable at weak coupling regime
    ( for \alpha >~ 1.685 )
*/
extern double Kennedy( double alpha );

#endif
