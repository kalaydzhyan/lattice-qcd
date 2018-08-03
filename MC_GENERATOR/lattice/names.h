#ifndef _NAMES_H_
#define _NAMES_H_

#include <string.h>
#include <geometry.h>

#define COUPLING_LIMIT  -100.0


/* Prepears file name pattern = b*_D*_T*_L* and stores its internally */
extern void  lattice_set_name( double beta, uchar T, uchar L, uchar D );

/* returns the above pattern concatenated with given ext: pattern + '.' + ext */
extern char *lattice_get_name( char *ext );


/* Parses file name of the type 'b*_D*_T*_L*' + '.' + ([^\.]+) + (\.ext)? 
   and returns $1 setting beta, T, L, D */
extern char *lattice_parse_name( char *name, char *ext,
				 double *beta, int *size0, int *size1, int *dim );
#endif
