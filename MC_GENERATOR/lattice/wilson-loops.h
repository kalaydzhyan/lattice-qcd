#ifndef _WILSON_LOOPS_H_
#define _WILSON_LOOPS_H_

#include <MC-SU2.h>
#include <bessel.h>
#include <strings.h>

#define  WILSON_LOOPS_VOLUME_SHIFT  1


extern double smearing_alpha;
extern uint   N_smearing;

extern void SU2_spacial_smearing( SU2 *G, uchar i0 );


/* Wilson loops calculation. Returns averaged on this configuration
   values of wilson loops (see however filename argument below) or
   NULL in case of failure. Arguments:
   - size: two uchars which will recieve the maximal temporal (size[0])
           and spatial (size[1]) extent of loops. If you don't need these
           (see filename below) just pass NULL pointer.
   - smearing: whether we have to smear given gauge configuration before
           measurements. Number of smearing sweeps is determined by this arg.
   - filename: if not NULL the values of wilson loops are appended to this file.
	       Format of each record: uchar size[0], uchar size[1], N_loops double values.
   'Temporal' blocking:
     0 - no blocking at all
     1 - use standard multihit *USE FOR WILSON ACTION ONLY*
     2 - use hypercubic blocking   */

extern double         *wilson_loops( uchar *size, uchar smearing, uchar blocking, char *filename );
extern double *adjoint_wilson_loops( uchar *size, uchar smearing, uchar blocking, char *filename );

extern double *wilson_loops_scan_smearing( char *filename );

/* Call this when done with loops - you can do this only once at very end */
extern void    wilson_loops_free();

/* Calculation of Polyakov lines. Returns averaged on given configuration
   values of <P(0) P(R)> or NULL on exception.
   uint *count will recieve the total number of PLine entries;
   double *line   -  will get |1/V3 \sum L_x| on this configuration (pass NULL if you don't need it)
   char *filename - if not NULL every time correlators are appended to that file
                    Format of each record:
		    uint N_PLine
		    N_PLine records { uint R2, double value0, double value1 }
   argument blocking - the same as for wilson loops  */
struct s_PLine {
  uint R2;
  double value;
};
typedef struct s_PLine  PLine;

extern PLine         *PLine_correlator( uint *count, double *line, uchar blocking, char *filename );
extern PLine *adjoint_PLine_correlator( uint *count, double *line, uchar blocking, char *filename );
extern void   PLine_correlator_free();
#endif
