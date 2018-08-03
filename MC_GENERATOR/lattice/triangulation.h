#ifndef _TRIANGULATION_H_
#define _TRIANGULATION_H_

#include <geometry.h>

#define MAX_PARITY    16

/* Tringulation of four dimensional hypercubical lattice consists of
   cells of various dimensions, which are:
     dim 0 (sites) -- same as on the original lattice;
     dim 1 (links) -- original links plus additional ones required by
                      triangulation;
     dim 2 (triangles) -- various triangles (bisected plaquettes etc...);
     dim 3 (simplices) -- .....
     dim 4 (hypersimplices) -- ....
   Below are the number of various cells per one lattice point */

#define          LINK_PER_SITE    11
#define      TRIANGLE_PER_SITE    34
#define       SIMPLEX_PER_SITE    40
#define  HYPERSIMPLEX_PER_SITE    16

/* Handy defines to enumerate distinct cells in triangulation */
#define    L_ENUM(i,m)       ((i) +         LINK_PER_SITE * (m))
#define   TR_ENUM(i,m)       ((i) +     TRIANGLE_PER_SITE * (m))
#define    S_ENUM(i,m)       ((i) +      SIMPLEX_PER_SITE * (m))
#define   HS_ENUM(i,m)       ((i) + HYPERSIMPLEX_PER_SITE * (m))

/* Differential on triangulated complex. We worked out only differential,
   not its adjoint boundary operator. Name of the function Tr_diff_'cell'
   means that it defines the action of differential on (D-1)-dim cell 'cell'
   with value in D-dim cells. Arguments:
     first (input):  index of target D-dim cell on which differential is
                     to be evaluated;
     second (output): indices of (D-1)-dim cells which enter the differential;
     third (output) : signs with which (D-1)-dim cells enter the differential.

   Typically, (diff (d-1)cell)_{d-cell} = \sum sign_k (d-1)-cell_k

   Number of required (d-1)-cells (hence the dimensionality of *idx and *sign)
   is trivial and is indicated below. */

extern void Tr_diff_simplex( uint hs, uint *s_idx, int *s_sign );   /* 5 simplices */
extern void Tr_diff_triangle( uint s, uint *tr_idx, int *tr_sign ); /* 4 triangles */
extern void Tr_diff_link( uint tr, uint *l_idx, int *l_sign );      /* 3 links */
extern void Tr_diff_site( uint l, uint *idx, int *sign );           /* 2 sites */

/* ------------------------------------------------------------------------ */
/* Following routines return indices of points entering given cell.
   Points orientation is in accord with differential above (see also below) */
/* Link orientation idx[0] --> idx[1], differential is idx[1] - idx[0]
   (which is the usual rule 'end minus start') */
extern void Tr_points_in_link( uint l_idx, uint *idx );

/* Triangle orientation  idx[0]-->idx[1]-->idx[2] */
extern void Tr_points_in_triangle( uint tr_idx, uint *idx );

/* Simplex orientation  idx[0]-->idx[1]-->idx[2]-->idx[3] */
extern void Tr_points_in_simplex( uint s_idx, uint *idx );

/* Hypersimplex orientation  idx[0]-->idx[1]-->idx[2]-->idx[3]-->idx[4] */
extern void Tr_points_in_hypersimplex( uint hs_idx, uint *idx );


/* ------------------------------------------------------------------------ */
/* Dependence of various cells on particular site. Number of dependent cells varies
   with site's parity, below are maximal number of dependent cells */

#define MAX_SITE_DEPENDENCE_LINK             48
#define MAX_SITE_DEPENDENCE_TRIANGLE        240
#define MAX_SITE_DEPENDENCE_SIMPLEX         384
#define MAX_SITE_DEPENDENCE_HYPERSIMPLEX    192

/* Functions returning indices of dependent cells via 'idx' and
   actual number of depening cells as a return value. */
extern uint Tr_link_dependent_on_site( uint m_s, uint *idx );
extern uint Tr_triangle_dependent_on_site( uint m_s, uint *idx );
extern uint Tr_simplex_dependent_on_site( uint m_s, uint *idx );
extern uint Tr_hypersimplex_dependent_on_site( uint m_s, uint *idx );



#define MAX_LINK_DEPENDENCE_TRIANGLE        14
extern uchar link2triangle_get( uint m_l, uint *idx, int *sign );


/* ------------------------------------------------------------------------ */
/* Restoration of normal hypercubical geometry (done only for three and four
   dimensional cells). Crudely, hypecube value is given by sum over some hypersimplices,
   and similar to cube construction from simplices. Routines below take
   (hyper) cube index as first argument and return indices and signs of corresponding
   (hyper) simplices. Number of (hyper) simplices is defined below.
   Orientations are such that if
      hypersimplex = d simplices
   then
      hypercube = diff cube
   where divergence 'diff' is sum of 'outgoing link minus incoming link' */

#define   HYPERSIMPLEX_PER_HYPERCUBE      HYPERSIMPLEX_PER_SITE
#define   SIMPLEX_PER_CUBE                5

extern void Tr_hypercube_from_hypersimplex( uint m, uint *idx, int *sign );
extern void Tr_cube_from_simplex( uint m, uint *idx, int *sign );


/* Performs cumulative checks of various routines */
extern void Tr_check_all();
#endif
