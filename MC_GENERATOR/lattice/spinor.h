#ifndef _SPINOR_H_
#define _SPINOR_H_

/* Compile-time defines:

   DEBUG                - include debugging code (might be slow)
   SPINOR_RANDOM_PHASE  - spinor phase is randomized as much as possible
                          (see 'spinor_rephasing_random' below)

   The following make 'spinor_rephasing_random' and friends (see below) effectively dummy.
   SPINOR_SECTION_PLUS  - all spinors are taken in 'plus'  form [1;z]
   SPINOR_SECTION_MINUS - all spinors are taken in 'minus' form [z;1]
*/
#include <strings.h>
#include <MC-SU2.h>

struct s_Spinor;
struct s_Bond;
struct s_Plaq;
struct s_Triangle;
struct s_Square;
struct s_Cube;
struct s_Site;
struct s_Diamond;
struct s_Hypercube;
struct s_Hypersite;

extern uint N_triangle;
extern uint N_square;
extern uint N_site;
extern uint N_cube;
extern uint N_diamond;
extern uint N_hypersite;
extern uint N_hypercube;
/* ::::::::::::::::::::::::::::::::::::::::::::::::::: */
struct s_Spinor {
  double complex z1, z2; /* z1 (z2) effectively useless
			    if SECTION_PLUS (SECTION_MINUS) is taken */
  struct s_Plaq *plaq;
  uchar self;            /* this = this->plaq->H[this->self] */
};
/* ::::::::::::::::::::::::::::::::::::::::::::::::::: */
struct s_Bond {
  SU2 *F;
  struct s_Plaq **plaq;  /* 2(D-1) plaq which share this bond */
  uchar *self;           /* this = this->plaq[k]->bond[ this->self[k] ]*/

  /* These depend on me, some of them might be NULL, number of entries - N_* above */
  struct s_Triangle  **triangle; 
  struct s_Square    **square;
  struct s_Site      **site;
  struct s_Cube      **cube;
  struct s_Diamond   **diamond;
  struct s_Hypersite **hypersite;
  struct s_Hypercube **hypercube;
};
/* ::::::::::::::::::::::::::::::::::::::::::::::::::: */
struct s_Plaq {           /* Positive orientation iff mu < nu */
  struct s_Spinor *Sp[4];
  struct s_Bond *bond[4];
  uchar  self[4];         /* this = this->bond[k]->plaq[ this->self[k] ] */
  double phase;
};
/* ::::::::::::::::::::::::::::::::::::::::::::::::::: */
struct s_Triangle {
  struct s_Spinor *Sp[3];
  double phase;
};
/* ::::::::::::::::::::::::::::::::::::::::::::::::::: */
struct s_Square {
  struct s_Spinor *Sp[4];
  double phase;
};
/* ::::::::::::::::::::::::::::::::::::::::::::::::::: */
struct s_Site {           /* Positive orientation == outgoing flux */
  char q;
  struct s_Triangle *triangle[8];
  struct s_Square   *square[6];
};
/* ::::::::::::::::::::::::::::::::::::::::::::::::::: */
struct s_Cube {           /* Positive orientation == outgoing flux */
  char q;
  struct s_Triangle *triangle[8];
  struct s_Plaq *plaq[6];
  char plaq_sign[6];
};
/* ::::::::::::::::::::::::::::::::::::::::::::::::::: */
struct s_Diamond {        /* Only for D=4, positive orientation == outgoing flux */
  char group;             /* Two group with different charge difinitions - first bit;
		             Charge sign convention - second bit (if set then sign is 'reversed') */
  char q;
  struct s_Triangle *triangle[8];
};
/* ::::::::::::::::::::::::::::::::::::::::::::::::::: */
struct s_Hypersite {
  struct s_Site    *site[4];
  struct s_Diamond *diamond[24];
};
/* ::::::::::::::::::::::::::::::::::::::::::::::::::: */
struct s_Hypercube {
  struct s_Diamond *diamond[16];
  struct s_Cube    *cube[8];
};
/* ::::::::::::::::::::::::::::::::::::::::::::::::::: */

typedef struct s_Spinor   Spinor;
typedef struct s_Bond     Bond;
typedef struct s_Plaq     Plaq;
typedef struct s_Triangle Triangle;
typedef struct s_Square   Square;
typedef struct s_Site     Site;
typedef struct s_Cube     Cube;
typedef struct s_Diamond  Diamond;
typedef struct s_Hypersite Hypersite;
typedef struct s_Hypercube Hypercube;

/* ********************************************************** */

extern Spinor    *spinor;
extern Bond      *bond;
extern Plaq      *plaq;
extern Triangle  *triangle;
extern Square    *square;
extern Site      *site;
extern Cube      *cube;
/* These are NULL in D=3 */
extern Diamond   *diamond;
extern Hypersite *hypersite;
extern Hypercube *hypercube;

/******************************************************************/
/******************************************************************/
extern uchar spinor_enum( uchar i0, char d0, uchar i1, char d1 );
extern uint  spinor_index( uchar i0, char d0, uchar i1, char d1, uint m );

extern int  spinor_build();
extern void spinor_free();

extern double density_site();
extern double density_cube();
extern double density_diamond();
extern double density_diamond0();
extern double density_diamond1();
extern double density_site_total();

extern void update_plaq( Plaq *p );
extern void update_triangle( Triangle *tr );
extern void update_square( Square *sq );
extern double measure_triangle( Triangle *tr );
extern double measure_square( Square *sq );
extern char site_charge( Site *s );
extern char cube_charge( Cube *c );
extern char diamond_charge( Diamond *diam );
extern char hypersite_charge( Hypersite *hs );
extern char hypercube_charge( Hypercube *hc );

#ifdef DEBUG
void spinor_print( Spinor *sp );
#endif

#define  TOTAL_BINS    60
extern double *triangle_histogram();
extern double *square_histogram();
extern double *plaq_histogram();
#endif
