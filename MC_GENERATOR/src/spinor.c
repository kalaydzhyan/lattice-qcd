#include <spinor.h>

/* Compile-time defines:

   DEBUG                - include debugging code (slow)
   SPINOR_RANDOM_PHASE  - spinor phase is randomized as much as possible
                          (see 'spinor_rephasing_random' below)

   The following make 'spinor_rephasing_random' and friends (see below) effectively dummy.

   SPINOR_SECTION_PLUS  - all spinors are taken in 'plus'  form [1;z]
   SPINOR_SECTION_MINUS - all spinors are taken in 'minus' form [z;1]
*/
#if (defined(SPINOR_RANDOM_PHASE))&&(defined(SPINOR_SECTION_PLUS))
#error "You CANNOT define both SPINOR_RANDOM_PHASE and SPINOR_SECTION_PLUS"
#endif
#if (defined(SPINOR_RANDOM_PHASE))&&(defined(SPINOR_SECTION_MINUS))
#error "You CANNOT define both SPINOR_RANDOM_PHASE and SPINOR_SECTION_MINUS"
#endif
#if (defined(SPINOR_SECTION_PLUS))&&(defined(SPINOR_SECTION_MINUS))
#error "You CANNOT define both SPINOR_SECTION_PLUS and SPINOR_SECTION_MINUS"
#endif
/* *********************************************************************** */
/* ************* Miscellaneous Spinor Stuff ****************************** */
/* *********************************************************************** */
static const char spinor_table_D3[6][6]={
  { -1, -1,  0,  1,  2,  3 }, { -1, -1,  4,  5,  6,  7 },
  {  0,  4, -1, -1,  8,  9 }, {  1,  5, -1, -1, 10, 11 },
  {  2,  6,  8, 10, -1, -1 }, {  3,  7,  9, 11, -1, -1 }
};
static const char spinor_table_D4[8][8]={
  { -1,  -1,  0,  1,  2,  3,  4,  5 }, { -1,  -1,  6,  7,  8,  9, 10, 11 },
  {  0,   6, -1, -1, 12, 13, 14, 15 }, {  1,   7, -1, -1, 16, 17, 18, 19 },
  {  2,   8, 12, 16, -1, -1, 20, 21 }, {  3,   9, 13, 17, -1, -1, 22, 23 },
  {  4,  10, 14, 18, 20, 22, -1, -1 }, {  5,  11, 15, 19, 21, 23, -1, -1 }
};
uchar spinor_enum( uchar i0, char d0, uchar i1, char d1 ){
  char ret = -1;
  i0 *= 2;
  i1 *= 2;
  i0 += ( d0 > 0 ) ? 1 : 0 ;
  i1 += ( d1 > 0 ) ? 1 : 0 ;
  if( param->D == 3 ) ret = spinor_table_D3[i0][i1];
  else                ret = spinor_table_D4[i0][i1];
  if( ret < 0 ) panic("Illegal call to spinor_enum");
  return (uchar) ret;
};
uint spinor_index( uchar i0, char d0, uchar i1, char d1, uint m ){
  return spinor_enum(i0,d0,i1,d1) + 2*param->D*(param->D-1)*m;
};
/* ----------------------------------------------------------------------- */
#ifdef DEBUG
void spinor_print( Spinor *sp ){
  uint m;
  uchar i0, i1, k, found = 0;
  Plaq *p;
  for( m = 0; m < param->ipw[param->D] ; m++  )
    for( i0 = 0; i0 < param->D-1 ; i0++ )
      for( i1 = i0+1; i1 < param->D ; i1++ ){
	p = &plaq[ plaq_index(i0, i1, m) ];
	for( k = 0; k < 4; k++ ) if( p->Sp[k] == sp ){ found = 1; goto final; };
      };
 final:
  if( !found ) panic("Spinor NOT found");
  if( k == 0 )      printf(" (%1d+, %1d+)", i0, i1 );
  else if( k == 1 ) printf(" (%1d-, %1d+)", i0, i1 );
  else if( k == 2 ) printf(" (%1d-, %1d-)", i0, i1 );
  else              printf(" (%1d+, %1d-)", i0, i1 );
  fflush( stdout );
};
#endif
/* ----------------------------------------------------------------------- */
#if (!defined(SPINOR_SECTION_PLUS))&&(!defined(SPINOR_SECTION_MINUS))
static inline void spinor_normalize( Spinor *sp ){
  double help = sqrt(cnorm(sp->z1) + cnorm(sp->z2));
  sp->z1 /= help; sp->z2 /= help;
};
/* ----------------------------------------------------------------------- */
static inline void spinor_random_phase( Spinor *sp ){
  double phase = (2.0*RND()-1.0)*M_PI;
  double complex zu = cos(phase) + I*sin(phase);
  sp->z1 *= zu;
  sp->z2 *= zu;
};
/* ----------------------------------------------------------------------- */
static inline void spinor_section_plus( Spinor *sp ){
  double complex zu = conj(sp->z1) / cabs( sp->z1 );
  sp->z1 *= zu;
  sp->z2 *= zu;
};
/* ----------------------------------------------------------------------- */
static inline void spinor_section_minus( Spinor *sp ){
  double complex zu = conj(sp->z2) / cabs( sp->z2 );
  sp->z1 *= zu;
  sp->z2 *= zu;
};
#endif
/* ----------------------------------------------------------------------- */
static inline double spinor_next( Spinor *sp, SU2 *V, Spinor *next ){
#if defined(SPINOR_SECTION_PLUS)
  double complex z = V->alpha - conj(V->beta) * sp->z2;
  next->z1 = 1.0;
  next->z2 = ( V->beta + conj(V->alpha) * sp->z2 )/z;
  return atan2( cimag(z), creal(z) );
#elif defined(SPINOR_SECTION_MINUS)
  double complex z = V->beta * sp->z1 + conj(V->alpha);
  next->z1 = ( V->alpha * sp->z1 - conj(V->beta) )/ z;
  next->z2 = 1.0;
  return atan2( cimag(z), creal(z) );
#else
  next->z1 = V->alpha * sp->z1 - conj( V->beta  ) * sp->z2;
  next->z2 = V->beta  * sp->z1 + conj( V->alpha ) * sp->z2;
  return 0.0;
#endif
};
/* ----------------------------------------------------------------------- */
static inline double spinor_product( Spinor *sp1, Spinor *sp2 ){
  double complex z;
#if defined(SPINOR_SECTION_PLUS)
  z = 1.0 + sp1->z2 * conj(sp2->z2);
#elif defined(SPINOR_SECTION_MINUS)
  z = sp1->z1 * conj(sp2->z1) + 1.0;
#else
  z = sp1->z1 * conj(sp2->z1) + sp1->z2 * conj(sp2->z2);
#endif
  return atan2( cimag(z), creal(z) );
};
/* ----------------------------------------------------------------------- */
static inline double spinor_scalar_product( Spinor *sp1, Spinor *sp2 ){
  double complex z;
#if defined(SPINOR_SECTION_PLUS)
  z = sp1->z2 - sp2->z2;
  return 1.0 - 2.0 * cnorm(z)/ ( (1.0 + cnorm(sp1->z2)) * (1.0 + cnorm(sp2->z2)) );
#elif defined(SPINOR_SECTION_MINUS)
  z = sp1->z1 - sp2->z1;
  return 1.0 - 2.0 * cnorm(z)/ ( (1.0 + cnorm(sp1->z1)) * (1.0 + cnorm(sp2->z1)) );
#else
  z = sp1->z2 * sp2->z1 - sp1->z1 * sp2->z2;
  return 1 - 2.0 * cnorm( z );
#endif
};
/* ----------------------------------------------------------------------- */
static inline double spinor_matrix_element( Spinor *sp1, SU2 *V, Spinor *sp2 ){
#if defined(SPINOR_SECTION_PLUS)
  double complex z = V->alpha - conj(V->beta) * sp1->z2;
  return atan2( cimag(z), creal(z) );
#elif defined(SPINOR_SECTION_MINUS)
  double complex z = V->beta * sp1->z1 + conj(V->alpha);
  return atan2( cimag(z), creal(z) );
#else
  Spinor sp;
  spinor_next( sp1, V, &sp );
  return spinor_product( &sp, sp2 );
#endif
};
/* ----------------------------------------------------------------------- */
#ifdef DEBUG
static int spinor_same_basepoint( Spinor *sp1, Spinor *sp2 ){
#if defined(SPINOR_SECTION_PLUS)
  if( cabs(sp1->z2 - sp2->z2) > 0.000001 ) return 0;
  return 1;
#elif defined(SPINOR_SECTION_MINUS)
  if( cabs(sp1->z1 - sp2->z1) > 0.000001 ) return 0;
  return 1;
#else
  double complex z1, z2;
  if( cabs(sp1->z2) > 0.5 ){
    z1 = sp1->z1 / sp1->z2;
    z2 = sp2->z1 / sp2->z2;
  }else{
    z1 = sp1->z2 / sp1->z1;
    z2 = sp2->z2 / sp2->z1;
  };
  if( cabs(z1 - z2) > 0.000001 ) return 0;
  return 1;
#endif
};
#endif
/* ----------------------------------------------------------------------- */
static inline void spinor_eigenstate( SU2 *V, Spinor *eigen ){
  double help = creal( V->alpha );
#if defined(SPINOR_SECTION_PLUS)
  eigen->z1 = 1.0;
  eigen->z2 = I*(cimag(V->alpha) - sqrt(1.0-help*help))/conj( V->beta );
#elif defined(SPINOR_SECTION_MINUS)
  eigen->z1 = I*(cimag(V->alpha) + sqrt(1.0-help*help))/V->beta;
  eigen->z2 = 1.0;
#else
  eigen->z1 = conj( V->beta );
  eigen->z2 = I * ( cimag(V->alpha) - sqrt(1.0-help*help) );
  spinor_normalize( eigen );
#ifdef SPINOR_RANDOM_PHASE
  spinor_random_phase( eigen );
#endif
#endif
#ifdef DEBUG
  {
    Spinor next;
#if defined(SPINOR_SECTION_PLUS)||defined(SPINOR_SECTION_MINUS)
    help = spinor_next( eigen, V, &next );
#else
    spinor_next( eigen, V, &next );
    help = spinor_product( &next, eigen );
#endif
    if( fabs(help - acos( creal(V->alpha) )) > 0.0000001 )
      panic("Eigenstate: failed in phase");
    if( !spinor_same_basepoint(eigen, &next) ) panic("Eigenstate: failed in basepoint");
  };
#endif
};
/* ----------------------------------------------------------------------- */
/* Assumes that spinors are already eigenstates (but U(1) phase is arbitrary) */
static inline void update_plaq_phase( Plaq *p ){
  p->phase = mod_pi(
		    spinor_matrix_element( p->Sp[0], p->bond[0]->F, p->Sp[1] ) +
		    spinor_matrix_element( p->Sp[1], p->bond[1]->F, p->Sp[2] ) -
		    spinor_matrix_element( p->Sp[0], p->bond[3]->F, p->Sp[3] ) -
		    spinor_matrix_element( p->Sp[3], p->bond[2]->F, p->Sp[2] ) );
};
/* ----------------------------------------------------------------------- */
void update_plaq( Plaq *p ){
  SU2 total = SU2_mult( *(p->bond[0]->F), *(p->bond[1]->F) );
  total = SU2_mult( total, SU2_conj( SU2_mult( *(p->bond[3]->F), *(p->bond[2]->F) ) ) );
  spinor_eigenstate( &total, p->Sp[0] );
#if defined(SPINOR_SECTION_PLUS)||defined(SPINOR_SECTION_MINUS)
  p->phase  = spinor_next( p->Sp[0], p->bond[0]->F, p->Sp[1] );
  p->phase += spinor_next( p->Sp[1], p->bond[1]->F, p->Sp[2] );
  p->phase -= spinor_next( p->Sp[0], p->bond[3]->F, p->Sp[3] );
  p->phase -= spinor_matrix_element( p->Sp[3], p->bond[2]->F, p->Sp[2] );
  p->phase = mod_pi( p->phase );
#else
  spinor_next( p->Sp[0], p->bond[0]->F, p->Sp[1] );
  spinor_next( p->Sp[1], p->bond[1]->F, p->Sp[2] );
  spinor_next( p->Sp[0], p->bond[3]->F, p->Sp[3] );
#ifdef SPINOR_RANDOM_PHASE
  spinor_random_phase( p->Sp[1] );
  spinor_random_phase( p->Sp[2] );
  spinor_random_phase( p->Sp[3] );
  update_plaq_phase( p );
#else
  p->phase = -spinor_matrix_element( p->Sp[3], p->bond[2]->F, p->Sp[2] );
#endif
#endif
#ifdef DEBUG
  {
    Spinor next;
    if( fabs( acos(creal(total.alpha)) - p->phase) > 0.0000001 ) panic("Wrong plaquette phase!");
    total = SU2_mult( total, *(p->bond[0]->F) );
    total = SU2_mult( SU2_conj(*(p->bond[0]->F)), total );
    spinor_eigenstate( &total, &next );
    if( !spinor_same_basepoint( p->Sp[1], &next) ) panic("Basepoint at 1");
    total = SU2_mult( total, *(p->bond[1]->F) );
    total = SU2_mult( SU2_conj(*(p->bond[1]->F)), total );
    spinor_eigenstate( &total, &next );
    if( !spinor_same_basepoint( p->Sp[2], &next) ) panic("Basepoint at 2");
    total = SU2_mult( total, SU2_conj(*(p->bond[2]->F)) );
    total = SU2_mult( *(p->bond[2]->F), total );
    spinor_eigenstate( &total, &next );
    if( !spinor_same_basepoint( p->Sp[3], &next) ) panic("Basepoint at 3");
  };
#endif
};
/* ----------------------------------------------------------------------- */
void update_triangle( Triangle *tr ){
  tr->phase = mod_pi(
		     spinor_product(tr->Sp[0], tr->Sp[1]) +
		     spinor_product(tr->Sp[1], tr->Sp[2]) +
		     spinor_product(tr->Sp[2], tr->Sp[0]) );
};
/* ----------------------------------------------------------------------- */
double measure_triangle( Triangle *tr ){
  double help, ret = 0.0;
  help = spinor_scalar_product( tr->Sp[0], tr->Sp[1] );
  ret += help * help;
  help = spinor_scalar_product( tr->Sp[1], tr->Sp[2] );
  ret += help * help;
  help = spinor_scalar_product( tr->Sp[2], tr->Sp[0] );
  ret += help * help;
  return ret;
};
/* ----------------------------------------------------------------------- */
void update_square( Square *sq ){
  sq->phase = mod_pi(
		     spinor_product(sq->Sp[0], sq->Sp[1]) +
		     spinor_product(sq->Sp[1], sq->Sp[2]) +
		     spinor_product(sq->Sp[2], sq->Sp[3]) +
		     spinor_product(sq->Sp[3], sq->Sp[0]) );
};
/* ----------------------------------------------------------------------- */
double measure_square( Square *sq ){
  double help, ret = 0.0;
  help = spinor_scalar_product( sq->Sp[0], sq->Sp[1] );
  ret += help * help;
  help = spinor_scalar_product( sq->Sp[1], sq->Sp[2] );
  ret += help * help;
  help = spinor_scalar_product( sq->Sp[2], sq->Sp[3] );
  ret += help * help;
  help = spinor_scalar_product( sq->Sp[3], sq->Sp[0] );
  ret += help * help;
  return ret;
};
/* ----------------------------------------------------------------------- */
static void spinor_build_on_plaq( uint n, uchar i0, uchar i1 ){
  uint m[4], l[4];
  uchar k;
  Plaq *p = &plaq[plaq_index(i0,i1,n)];
  if( i0 > i1 ){ k = i1 ; i1 = i0; i0 = k; };
  m[0] = n; m[1] = index_up(n,i0); m[2] = index_up(m[1],i1); m[3] = index_up(n,i1);
  l[0] = link_index( i0, m[0] ); l[1] = link_index( i1, m[1] );
  l[2] = link_index( i0, m[3] ); l[3] = link_index( i1, m[0] );
  p->Sp[0] = &spinor[ spinor_index(i0, 1,i1, 1,m[0]) ];
  p->Sp[1] = &spinor[ spinor_index(i0,-1,i1, 1,m[1]) ];
  p->Sp[2] = &spinor[ spinor_index(i0,-1,i1,-1,m[2]) ];
  p->Sp[3] = &spinor[ spinor_index(i0, 1,i1,-1,m[3]) ];
  for( n = 0; n < 4; n++ ){
    p->Sp[n]->plaq = p;
    p->Sp[n]->self = n;
    p->bond[n] = &bond[l[n]];
    for( k = 0 ; k < 2*(param->D-1) && p->bond[n]->plaq[k]; k++ );
    if( k == 2*(param->D-1) ) panic("Plaq didn't find empty place in bond->plaq references");
    p->bond[n]->plaq[k] = p;
    p->self[n] = k;
    p->bond[n]->self[k] = n;
  };
  update_plaq( p );
};
/* *********************************************************************** */
/* ***************  Real Stuff - Spinor Construction  ******************** */
/* *********************************************************************** */
/*     BOND REFERENCES:                    |     PLAQ REFERENCES:
                    D=4         D=3        |                   D=4    D=3
   N_triangle       232         24         |  Np_triangle      48     8
   N_square         78          26         |  Np_square        16     8
   N_site           30          10         |  Np_site          8      4
   N_cube           12          4          |  Np_cube          4      2
   N_diamond        102         0          |  Np_diamond       24     0
   N_hypersite      14          0          |  Np_hypercube     4      0
   N_hypercube      8           0          |  Np_hypersite     4      0
-------------------------------------------------------------------------------
   Note1: In D=4 N?_hypersite = N?_hypercube = 0 unless defined DEBUG
   Note2: In D=4 among N?_diamond dependent diamonds:
             bond:  64 belong to group 0, 38 to group 1;
             plaq:  16 belong to group 0,  8 to group 1;
------------------------------------------------------------------------------- */
Spinor    *spinor = NULL;
Bond      *bond  = NULL;
Plaq      *plaq  = NULL;
Triangle  *triangle  = NULL;
Square    *square  = NULL;
Site      *site  = NULL;
Cube      *cube  = NULL;
Diamond   *diamond  = NULL;
Hypersite *hypersite  = NULL;
Hypercube *hypercube  = NULL;

uint N_triangle = 0;
uint N_square = 0;
uint N_site = 0;
uint N_cube = 0;
uint N_diamond = 0;
uint N_hypersite = 0;
uint N_hypercube = 0;

static uint current_triangle = 0;
static uint current_square   = 0;
static uint current_site     = 0;
static uint current_cube     = 0;
static uint current_diamond  = 0;

static void site_init( uint m );
static void cube_init( uint m );
static void hypersite_init( uint m );
static void hypercube_init();
/* ----------------------------------------------------------------------- */
void spinor_free(){
  uint m;
  current_triangle = current_square = current_site = current_cube = current_diamond = 0;
  if( spinor ){ free( spinor ); spinor = NULL; };
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( bond ){
    for( m = 0; m < param->D*param->ipw[param->D]; m++  ){
      if( bond[m].plaq ) free( bond[m].plaq );
      if( bond[m].self ) free( bond[m].self );
      if( bond[m].triangle )  free( bond[m].triangle );
      if( bond[m].square )    free( bond[m].square );
      if( bond[m].site )      free( bond[m].site );
      if( bond[m].cube )      free( bond[m].cube );
      if( bond[m].diamond )   free( bond[m].diamond );
      if( bond[m].hypersite ) free( bond[m].hypersite );
      if( bond[m].hypercube ) free( bond[m].hypercube );
    };
    free( bond ); bond = NULL;
  };
  if( plaq      ){ free( plaq );      plaq      = NULL; };
  if( triangle  ){ free( triangle );  triangle  = NULL; };
  if( square    ){ free( square );    square    = NULL; };
  if( site      ){ free( site );      site      = NULL; };
  if( cube      ){ free( cube );      cube      = NULL; };
  if( diamond   ){ free( diamond );   diamond   = NULL; };
  if( hypersite ){ free( hypersite ); hypersite = NULL; };
  if( hypercube ){ free( hypercube ); hypercube = NULL; };
};
/* ----------------------------------------------------------------------- */
int spinor_build( ){
  uint m;
  uint scale = param->D*(param->D-1);
  uchar i0, i1;
  Bond *b;
  if( !param || !F ) panic("No parameters/gauge data");
  if( param->D != 3 && param->D != 4 ) panic("This works only in D=3,4");
  if( spinor ){ error("Spinors already constructed");  return 0; };
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef VERBOSE
  printf("<<===  Spinor Construction [compiled with ");
#if defined(SPINOR_RANDOM_PHASE)
  printf("SPINOR_RANDOM_PHASE");
#elif defined(SPINOR_SECTION_PLUS)
  printf("SPINOR_SECTION_PLUS");
#elif defined(SPINOR_SECTION_MINUS)
  printf("SPINOR_SECTION_MINUS");
#else
  printf("no special options");
#endif
  printf("]  ===>>\n");
  fflush( stdout );
#endif
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( param->D == 3 ){
    N_triangle = 24;  N_square = 26; N_site = 10; N_cube = 4;
  }else{
    N_triangle = 232; N_square = 78; N_site = 30; N_cube = 12;   N_diamond = 102;
#ifdef DEBUG
    N_hypersite = 14; N_hypercube = 8;
#endif
  };
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( !spinor ) MEMORY_CHECK( spinor=(Spinor *)calloc( 2 * scale * param->ipw[param->D],sizeof(Spinor)),-1);
  if( !bond ){
    MEMORY_CHECK( bond = (Bond *) calloc( param->D * param->ipw[param->D], sizeof(Bond)),-1);
    for( m = 0; m < param->D * param->ipw[param->D]; m++  ){
      b = &bond[m];
      b->F  = &F[m];
      MEMORY_CHECK(b->plaq = (Plaq **)calloc( 2 * (param->D-1), sizeof(Plaq *)),-1);
      MEMORY_CHECK(b->self = (uchar *)calloc( 2 * (param->D-1), sizeof(uchar)) ,-1);
      if( N_triangle ) MEMORY_CHECK(b->triangle =(Triangle **) calloc( N_triangle,  sizeof(Triangle *))  ,-1);
      if( N_square )   MEMORY_CHECK(b->square   =(Square **)   calloc( N_square,    sizeof(Square *))    ,-1);
      if( N_site )     MEMORY_CHECK(b->site     =(Site **)     calloc( N_site,      sizeof(Site *))      ,-1);
      if( N_cube )     MEMORY_CHECK(b->cube     =(Cube **)     calloc( N_cube,      sizeof(Cube *))      ,-1);
      if( N_diamond )  MEMORY_CHECK(b->diamond  =(Diamond **)  calloc( N_diamond,   sizeof(Diamond *))   ,-1);
      if( N_hypersite) MEMORY_CHECK(b->hypersite=(Hypersite **)calloc( N_hypersite, sizeof(Hypersite *)) ,-1);
      if( N_hypercube) MEMORY_CHECK(b->hypercube=(Hypercube **)calloc( N_hypercube, sizeof(Hypercube *)) ,-1);
    };
  };
  scale /= 2;
  if( !plaq   ) MEMORY_CHECK( plaq=(Plaq *) calloc( scale * param->ipw[param->D], sizeof(Plaq)), -1);
  if(!triangle) MEMORY_CHECK(triangle=(Triangle *)calloc(((param->D==3)? 8 : 96) * param->ipw[param->D], sizeof(Triangle)),-1);
  if( !square ) MEMORY_CHECK( square =( Square * )calloc(((param->D==3)? 6 : 24) * param->ipw[param->D], sizeof(Square)),-1);
  if( !site   ) MEMORY_CHECK(  site  =(  Site *  )calloc(((param->D==3)? 1 : 4 ) * param->ipw[param->D],sizeof(Site)),-1);
  if( !cube   ) MEMORY_CHECK(  cube  =(  Cube *  )calloc(((param->D==3)? 1 : 4 ) * param->ipw[param->D],sizeof(Cube)),-1);
  if( param->D == 4 ){
    if( !diamond )   MEMORY_CHECK( diamond =( Diamond * )calloc( 24 * param->ipw[param->D],sizeof( Diamond )),-1);
    if( !hypersite ) MEMORY_CHECK(hypersite=(Hypersite *)calloc(      param->ipw[param->D],sizeof(Hypersite)),-1);
    if( !hypercube ) MEMORY_CHECK(hypercube=(Hypercube *)calloc(      param->ipw[param->D],sizeof(Hypercube)),-1);
  };
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for( m = 0; m < param->ipw[param->D]; m++  )
    for( i0 = 0; i0 < param->D-1; i0++ )
      for( i1 = i0+1; i1 < param->D; i1++ )
	spinor_build_on_plaq(m,i0,i1);
  if( param->D == 3 ){
    for( m = 0; m < param->ipw[param->D]; m++  ) site_init( m );
    for( m = 0; m < param->ipw[param->D]; m++  ) cube_init( m );
  }else{
    for( m = 0; m < param->ipw[param->D]; m++  ) hypersite_init( m );
    hypercube_init();
  };
#ifdef DEBUG
  {
    int total_charge = 0;
    if( param->D == 3 ){
      for( m = 0; m < param->ipw[param->D]; m++  ){
	  total_charge += site[m].q;
	  total_charge += cube[m].q;
      };
      if( total_charge ) panic("Non-zero total charge");
      if( current_triangle != 8*param->ipw[param->D] ) panic("You are wasting memory on triangles");
      if( current_square   != 6*param->ipw[param->D] ) panic("You are wasting memory on squares");
    }else{
      for( m = 0; m <  4*param->ipw[param->D]; m++  ) total_charge += site[m].q;
      for( m = 0; m <  4*param->ipw[param->D]; m++  ) total_charge += cube[m].q;
      for( m = 0; m < 24*param->ipw[param->D]; m++  ) total_charge += diamond[m].q;
      if( total_charge ) panic("Non-zero total charge");
      if( current_triangle != 96*param->ipw[param->D] ) panic("You are wasting memory on triangles");
      if( current_square   != 24*param->ipw[param->D] ) panic("You are wasting memory on squares");
      if( current_site     != 4*param->ipw[param->D]  ) panic("You are wasting memory on sites");
      if( current_cube     != 4*param->ipw[param->D]  ) panic("You are wasting memory on cubes");
      if( current_diamond  != 24*param->ipw[param->D] ) panic("You are wasting memory on diamonds");
    };
    scale = (uint)( RND() * ((double) (param->D * param->ipw[param->D])) );
    for( m = 0; m < N_triangle; m++ ) if( !bond[scale].triangle[m] ) panic("Wasting memory in bond triangle references");
    for( m = 0; m < N_square  ; m++ ) if( !bond[scale].square[m]   ) panic("Wasting memory in bond square references");
    for( m = 0; m < N_site    ; m++ ) if( !bond[scale].site[m]     ) panic("Wasting memory in bond site references");
    for( m = 0; m < N_cube    ; m++ ) if( !bond[scale].cube[m]     ) panic("Wasting memory in bond cube references");
    for( m = 0; m < N_diamond ; m++ ) if( !bond[scale].diamond[m]  ) panic("Wasting memory in bond diamond references");
    for( m = 0; m < N_hypersite;m++ ) if( !bond[scale].hypersite[m]) panic("Wasting memory in bond hypersite references");
    for( m = 0; m < N_hypercube;m++ ) if( !bond[scale].hypercube[m]) panic("Wasting memory in bond hypercube references");
  };
#endif
  return 0;
};
/* ----------------------------------------------------------------------- */
static void insert_triangle( Triangle *tr ){
  uchar i, j, k;
  Bond *b;
  Plaq *p;
  for( i = 0; i < 3; i++ ){
    p = tr->Sp[i]->plaq;
    for( j = 0; j < 4; j++ ){
      b = p->bond[j];
      for( k = 0; k < N_triangle && b->triangle[k] ; k++ ) if( b->triangle[k] == tr ) break;
      if( k == N_triangle ) panic("No place for triangle in bond references");
      if( !b->triangle[k] ) b->triangle[k] = tr;
    };
  };
};
/* ----------------------------------------------------------------------- */
static void insert_square( Square *sq ){
  uchar i, j, k;
  Bond *b;
  Plaq *p;
  for( i = 0; i < 4; i++ ){
    p = sq->Sp[i]->plaq;
    for( j = 0; j < 4; j++ ){
      b = p->bond[j];
      for( k = 0; k < N_square && b->square[k] ; k++ ) if( b->square[k] == sq ) break;
      if( k == N_square ) panic("No place for square in bond references");
      if( !b->square[k] ) b->square[k] = sq;
    };
  };
};
/* ----------------------------------------------------------------------- */
static void insert_site( Site *s ){
  uchar i, j, k, n;
  Bond *b;
  Plaq *p;
  for( i = 0 ; i < 8 ; i++ )
    for( j = 0 ; j < 3 ; j++ ){
      p = s->triangle[i]->Sp[j]->plaq;
      for( k = 0 ; k < 4 ; k++ ){
        b = p->bond[k];
	for( n = 0; n < N_site && b->site[n] ; n++ ) if( b->site[n] == s ) break;
	if( n == N_site ) panic("No place for site in bond references");
	if( !b->site[n] ) b->site[n] = s;
      };
    };
};
/* ----------------------------------------------------------------------- */
static void insert_cube( Cube *c ){
  uchar i, j, k, n;
  Bond *b;
  Plaq *p;
  for( i = 0 ; i < 8 ; i++ )
    for( j = 0 ; j < 3 ; j++ ){
      p = c->triangle[i]->Sp[j]->plaq;
      for( k = 0 ; k < 4 ; k++ ){
        b = p->bond[k];
	for( n = 0; n < N_cube && b->cube[n] ; n++ ) if( b->cube[n] == c ) break;
	if( n == N_cube ) panic("No place for cube in bond references");
	if( !b->cube[n] ) b->cube[n] = c;
      };
    };
};
/* ----------------------------------------------------------------------- */
static void insert_diamond( Diamond *diam ){ /* For D=4 only */
  uchar i, j, k, n;
  Bond *b;
  Plaq *p;
  for( i = 0 ; i < 8 ; i++ )
    for( j = 0 ; j < 3 ; j++ ){
      p = diam->triangle[i]->Sp[j]->plaq;
      for( k = 0 ; k < 4 ; k++ ){
        b = p->bond[k];
        for( n = 0; n < N_diamond && b->diamond[n] ; n++ ) if( b->diamond[n] == diam ) break;
        if( n == N_diamond ) panic("No place for diamond in bond references");
        if( !b->diamond[n] ) b->diamond[n] = diam;
      };
    };
};
/* ----------------------------------------------------------------------- */
static void insert_hypersite( Hypersite *hs ){ /* For D=4 only */
#ifdef DEBUG
  uchar i, j, k, n, m;
  Bond *b;
  Plaq *p;
  for( i = 0 ; i < 24 ; i++ )
    for( j = 0 ; j < 8 ; j++ )
      for( k = 0 ; k < 3 ; k++ ){
	p = hs->diamond[i]->triangle[j]->Sp[k]->plaq;
	for( n = 0 ; n < 4 ; n++ ){
	  b = p->bond[n];
	  for( m = 0; m < N_hypersite && b->hypersite[m] ; m++ ) if( b->hypersite[m] == hs ) break;
	  if( m == N_hypersite ) panic("No place for hypersites in bond references");
	  if( !b->hypersite[m] ) b->hypersite[m] = hs;
	};
      };
#endif
};
/* ----------------------------------------------------------------------- */
static void insert_hypercube( Hypercube *hc ){ /* For D=4 only */
#ifdef DEBUG
  uchar i, j, k, n, m;
  Bond *b;
  Plaq *p;
  for( i = 0 ; i < 16 ; i++ )
    for( j = 0 ; j < 8 ; j++ )
      for( k = 0 ; k < 3 ; k++ ){
	p = hc->diamond[i]->triangle[j]->Sp[k]->plaq;
	for( n = 0 ; n < 4 ; n++ ){
	  b = p->bond[n];
	  for( m = 0; m < N_hypercube && b->hypercube[m] ; m++ ) if( b->hypercube[m] == hc ) break;
	  if( m == N_hypercube ) panic("No place for hypercubes in bond references");
	  if( !b->hypercube[m] ) b->hypercube[m] = hc;
	};
      };
#endif
};
/* ----------------------------------------------------------------------- */
static inline uchar three_shifts( char d1, char d2, char d3 ){
  uchar ret = 0;
  if( d1 > 0 ) ret |= 0x04;
  if( d2 > 0 ) ret |= 0x02;
  if( d3 > 0 ) ret |= 0x01;
  return ret;
};
static inline uchar four_shifts( char d0, char d1, char d2, char d3 ){
  uchar ret = three_shifts( d1, d2, d3 );
  if( d0 > 0 ) ret |= 0x08;
  return ret;
};
/* ----------------------------------------------------------------------- */
char site_charge( Site *s ){
  uchar k;
  double phase = 0;
  if( param->D == 3 ){
    for( k = 0; k < 8 ; k++ ) phase += s->triangle[k]->phase;
    for( k = 0; k < 6 ; k++ ) phase += s->square[k]->phase;
  }else{
    char d[3];
    for( k = 0; k < 6 ; k++ ) phase += s->square[k]->phase;
    for( d[0] = -1; d[0] < 2; d[0] += 2 )
      for( d[1] = -1; d[1] < 2; d[1] += 2 )
	for( d[2] = -1; d[2] < 2; d[2] += 2 )
	  phase -= d[0]*d[1]*d[2]* s->triangle[ three_shifts(d[0], d[1], d[2]) ]->phase;
  };
#ifdef DEBUG
  if(fabs(rint(phase/(2*M_PI))-(phase/(2*M_PI)))>0.000001) panic("SITE Charge is not integer!");
#endif
  s->q = (char)(rint( phase/(2*M_PI) ));
  return s->q;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
char cube_charge( Cube *c ){
  uchar k;
  double phase = 0;
  for( k = 0 ; k < 6 ; k++ ) phase += c->plaq_sign[k] * c->plaq[k]->phase;
  if( param->D == 3 ){
    for( k = 0 ; k < 8 ; k++ ) phase -= c->triangle[k]->phase;
  }else{
    char d[3];
    for( d[0] = -1; d[0] < 2; d[0] += 2 )
      for( d[1] = -1; d[1] < 2; d[1] += 2 )
	for( d[2] = -1; d[2] < 2; d[2] += 2 )
	  phase -= d[0] * d[1] * d[2] * c->triangle[three_shifts(d[0],d[1],d[2])]->phase;
  };
#ifdef DEBUG
  if(fabs(rint(phase/(2*M_PI))-(phase/(2*M_PI)))>0.000001) panic("CUBE Charge is not integer!");
#endif
  c->q = (char)(rint( phase/(2*M_PI) ));
  return c->q;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
char diamond_charge( Diamond *diam ){
  double phase = 0.0;
  if( (diam->group & 0x01) ){
    char d[3];
    for( d[0] = -1; d[0] < 2; d[0] += 2 )
      for( d[1] = -1; d[1] < 2; d[1] += 2 )
	for( d[2] = -1; d[2] < 2; d[2] += 2 )
	  phase += d[0]*d[1]*d[2] * diam->triangle[three_shifts(d[0],d[1],d[2])]->phase;
  }else{
    uchar k;
    for( k = 0 ; k < 8 ; k++ ) phase += diam->triangle[k]->phase;
  };
#ifdef DEBUG
  if(fabs(rint(phase/(2*M_PI))-(phase/(2*M_PI)))>0.000001) panic("DIAMOND Charge is not integer!");
#endif
  diam->q = (char)(rint( phase/(2*M_PI) ));
  if( (diam->group & 0x02) ) diam->q *= -1; /* make diamond charges in hypersite consistent */
  return diam->q;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
char hypersite_charge( Hypersite *hs ){
#ifdef DEBUG
  uchar k;
  char charge = 0;
  for( k = 0; k < 24; k++ ) charge += hs->diamond[k]->q;
  /* YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY */
  //  for( k = 0; k < 4; k++ ) charge += hs->site[k]->q;
  /* YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY */
  if( charge ) panic("Hey, man, your hypersite is CHARGED!");
#endif
  return 0;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
char hypercube_charge( Hypercube *hc ){
#ifdef DEBUG
  uchar k;
  char charge = 0;
  for( k = 0; k < 16; k++ ) charge -= hc->diamond[k]->q;
  for( k = 0; k < 4 ; k++ ) charge -= hc->cube[k]->q;
  for( k = 4; k < 8 ; k++ ) charge += hc->cube[k]->q;
  if( charge ) panic("Hey, man, your hypercube is CHARGED!");
#endif
  return 0;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static void site_init( uint m ){ /* Only in D=3 */
  uchar i0, i1, i2;
  char d0, d1, d2;
  Triangle *tr;
  Square *sq;
  Site *s = &site[m];
  if( param->D != 3 ) panic("Illegal call");
  for( d0 = -1; d0 < 2; d0 += 2 )
    for( d1 = -1; d1 < 2; d1 += 2 )
      for( d2 = -1; d2 < 2; d2 += 2 ){
	i0 = three_shifts( d0, d1, d2 );
#ifdef DEBUG
	if( current_triangle == 8*param->ipw[param->D] ) panic("Wrong triangle count");
#endif
	tr = s->triangle[i0] = &triangle[current_triangle++];
	tr->Sp[0]  = &spinor[spinor_index( 0, d0, 1, d1, m)];
	if( d0 * d1 * d2 > 0 ){
	  tr->Sp[1]  = &spinor[spinor_index( 1, d1, 2, d2, m)];
	  tr->Sp[2]  = &spinor[spinor_index( 2, d2, 0, d0, m)];
	}else{
	  tr->Sp[2]  = &spinor[spinor_index( 1, d1, 2, d2, m)];
	  tr->Sp[1]  = &spinor[spinor_index( 2, d2, 0, d0, m)];
	};
	update_triangle( tr );
	insert_triangle( tr );
      };
  for( i0 = 0; i0 < 3; i0++ ){
    D3dual( i0, &i1, &i2 );
    for( d0 = 0; d0 < 2; d0++ ){
#ifdef DEBUG
      if( current_square == 6 * param->ipw[param->D] ) panic("Wrong square count");
#endif
      sq = s->square[ i0 + 3*d0 ] = &square[current_square++];
      sq->Sp[0] = &spinor[spinor_index( i0, d0, i1, 1, m)];
      sq->Sp[2] = &spinor[spinor_index( i0, d0, i1,-1, m)];
      if( d0 ){
	sq->Sp[1] = &spinor[spinor_index( i0, d0, i2, 1, m)];
	sq->Sp[3] = &spinor[spinor_index( i0, d0, i2,-1, m)];
      }else{
	sq->Sp[3] = &spinor[spinor_index( i0, d0, i2, 1, m)];
	sq->Sp[1] = &spinor[spinor_index( i0, d0, i2,-1, m)];
      };
      update_square( sq );
      insert_square( sq );
    };
  };
  site_charge(s);
  insert_site( s );
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static void hypersite_init( uint m ){ /* Only in D=4 */
  char d[4];
  uchar i0, i1, i2, i3, num;
  Diamond *diam;
  Triangle *tr;
  Site *s;
  Square *sq;
  Hypersite *hs = &hypersite[m];
  if( param->D != 4 ) panic("Illegal call");
  for( d[0] = -1; d[0] < 2; d[0] += 2 )
    for( d[1] = -1; d[1] < 2; d[1] += 2 )
      for( d[2] = -1; d[2] < 2; d[2] += 2 )
	for( d[3] = -1; d[3] < 2; d[3] += 2 ){
	  num = four_shifts( d[0], d[1], d[2], d[3] ); /* Diamond we're working with */
	  if( !hs->diamond[num] ){
#ifdef DEBUG
	    if( current_diamond == 24 * param->ipw[param->D] ) panic("Wrong diamond count");
#endif
	    hs->diamond[num] = &diamond[ current_diamond++ ];
	  };
	  diam = hs->diamond[num];
	  for( i0 = 0; i0 < 4; i0++ ){
	    i1 = (i0) ? 0 : 1 ;
	    D4dual( i0, i1, &i2, &i3 );
	    // - - - Triangle which goes to another group - - - - - - 
#ifdef DEBUG
	    if( current_triangle == 96 * param->ipw[param->D] ) panic("Wrong triangle count");
	    if( diam->triangle[ 2*i0 ] ) panic("Someone already settled MY triangle");
#endif
	    tr = diam->triangle[ 2*i0 ] = &triangle[current_triangle++];
	    tr->Sp[0] = &spinor[spinor_index( i0, d[i0], i1, d[i1], m)];
	    tr->Sp[1] = &spinor[spinor_index( i0, d[i0], i2, d[i2], m)];
	    tr->Sp[2] = &spinor[spinor_index( i0, d[i0], i3, d[i3], m)];
	    update_triangle(tr);
	    insert_triangle( tr );
	    num = 16 + i0; /* Num of diamond and its triangle num in another group */
	    if( d[i0] < 0 ) num += 4;
	    if( !hs->diamond[num] ){
#ifdef DEBUG
	      if( current_diamond == 24 * param->ipw[param->D] ) panic("Wrong diamond count");
#endif
	      hs->diamond[num] = &diamond[ current_diamond++ ];
	    };
#ifdef DEBUG
	    if( hs->diamond[num]->triangle[three_shifts( d[i1], d[i2], d[i3] )] )
	      panic("Someone already settled MY 'slave' triangle");
#endif
	    hs->diamond[num]->triangle[three_shifts( d[i1], d[i2], d[i3] )] = tr;

	    // - - - Triangle which goes to the same group and to site - - - - - - 
	    if( !(tr = diam->triangle[ 2*i0 + 1 ]) ){
#ifdef DEBUG
	      if( current_triangle == 96 * param->ipw[param->D] ) panic("Wrong triangle count");
#endif
	      tr = diam->triangle[ 2*i0 + 1 ] = &triangle[current_triangle++];
	      tr->Sp[0] = &spinor[spinor_index( i1, d[i1], i2, d[i2], m)];
	      tr->Sp[2] = &spinor[spinor_index( i2, d[i2], i3, d[i3], m)];
	      tr->Sp[1] = &spinor[spinor_index( i3, d[i3], i1, d[i1], m)];
	      update_triangle(tr);
	      insert_triangle( tr );
	      d[i0] *= -1; /* Num of other diamond in the same group */
	      num = four_shifts( d[0], d[1], d[2], d[3] );
	      d[i0] *= -1;
	      if( !hs->diamond[num] ){
#ifdef DEBUG
		if( current_diamond == 24 * param->ipw[param->D] ) panic("Wrong diamond count");
#endif
		hs->diamond[num] = &diamond[ current_diamond++ ];
	      };
#ifdef DEBUG
	      if( hs->diamond[num]->triangle[ 2*i0 + 1 ] )
		panic("Someone already settled MY 'slave' triangle");
#endif
	      hs->diamond[num]->triangle[ 2*i0 + 1 ] = tr;
	    };
	    if( !hs->site[i0] ){
#ifdef DEBUG
	      if( current_site == 4 * param->ipw[param->D] ) panic("Wrong site count");
#endif
	      hs->site[i0] = &site[ current_site++ ];
	    };
	    s = hs->site[i0];
	    num = three_shifts( d[i1], d[i2], d[i3] );
	    if( !s->triangle[num] && d[i0] > 0 ) s->triangle[num] = tr;
	  }; /* i0 */
	  diam->group = 0;
	  if( d[0] * d[1] * d[2] * d[3] < 0 ) diam->group |= 0x02;
	  diamond_charge( diam );
	  insert_diamond(diam);
	};   /* d0, d1, d2, d3 */
  // - finally construct sites ---------
  for( i0 = 0; i0 < 4; i0++ ){
    s = hs->site[i0];    num = 0; /* count of site's squares */
    for( i1 = 0; i1 < 4; i1++ ){
      if( i1 == i0 ) continue;
      D4dual( i0, i1, &i2, &i3 );
      for( d[i1] = 0; d[i1] < 2; d[i1]++ ){
#ifdef DEBUG
	if( s->square[num] ) panic("Someone already settled MY square");
	if( current_square == 24 * param->ipw[param->D] ) panic("Wrong square count");
#endif
	sq = s->square[num] = &square[current_square++];
	sq->Sp[0] = &spinor[spinor_index( i1, d[i1], i2, 1, m)];
	sq->Sp[2] = &spinor[spinor_index( i1, d[i1], i2,-1, m)];
	if( d[i1] ){
	  sq->Sp[1] = &spinor[spinor_index( i1, d[i1], i3, 1, m)];
	  sq->Sp[3] = &spinor[spinor_index( i1, d[i1], i3,-1, m)];
	}else{
	  sq->Sp[3] = &spinor[spinor_index( i1, d[i1], i3, 1, m)];
	  sq->Sp[1] = &spinor[spinor_index( i1, d[i1], i3,-1, m)];
	};
	update_square( sq );
	insert_square( sq );
	num++;
      };
    };
    site_charge(s);
    insert_site( s );
  };
#ifdef DEBUG
  for( num = 0; num < 24; num++ )
    for( i0 = 0; i0 < 8; i0++ ){
      if( !(hs->diamond[num]->triangle[i0]) ) panic("NULL triangle");
      i3 = 0;
      for( i1 = 0; i1 < 24; i1++ )
	for( i2 = 0; i2 < 8; i2++ ){
	  if( !(hs->diamond[i1]->triangle[i2]) ) panic("NULL triangle");
	  if( hs->diamond[num]->triangle[i0] == hs->diamond[i1]->triangle[i2] ) i3++;
	};
      if( i3 != 2 ) panic("Hypersite construction failed");
    };
#endif
  /* Second group of diamonds is still not inserted into bonds */
  for( i0 = 0; i0 < 4 ; i0++ )
    for( d[0] = -1; d[0] < 2; d[0] += 2 ){
      num  = 16 + i0;
      if( d[0] < 0 ) num += 4;
      diam = hs->diamond[num];
#ifdef DEBUG
      if( !diam ) panic("NULL diamond in second group");
#endif
      diam->group = (char) 0x01;
      if( d[0] > 0 ) diam->group |= 0x02;
      diamond_charge( diam );
      insert_diamond( diam );
    };
  hypersite_charge(hs);
  insert_hypersite(hs);
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static void cube_init( uint m ){ /* For D=3 only */
  uchar i1 = 0, i2 = 1, i3 = 2;
  uint n[8];
  Cube *c = &cube[m];
  if( param->D != 3 ) panic("Illegal call");
  n[0] = m;
  n[1] = index_up(n[0],i3);  n[2] = index_up(n[0],i2);
  n[3] = index_up(n[2],i3);  n[4] = index_up(n[0],i1);
  n[5] = index_up(n[4],i3);  n[6] = index_up(n[4],i2);
  n[7] = index_up(n[6],i3);
  c->plaq[0] = &plaq[ plaq_index(i2,i3,n[0]) ]; c->plaq_sign[0] = (i3 < i2) ? 1 : -1;
  c->plaq[1] = &plaq[ plaq_index(i2,i3,n[4]) ]; c->plaq_sign[1] = (i2 < i3) ? 1 : -1;
  c->plaq[2] = &plaq[ plaq_index(i3,i1,n[0]) ]; c->plaq_sign[2] = (i1 < i3) ? 1 : -1;
  c->plaq[3] = &plaq[ plaq_index(i3,i1,n[2]) ]; c->plaq_sign[3] = (i3 < i1) ? 1 : -1;
  c->plaq[4] = &plaq[ plaq_index(i1,i2,n[0]) ]; c->plaq_sign[4] = (i2 < i1) ? 1 : -1;
  c->plaq[5] = &plaq[ plaq_index(i1,i2,n[1]) ]; c->plaq_sign[5] = (i1 < i2) ? 1 : -1;
  for( i1 = 0 ; i1 < 8 ; i1++ ) c->triangle[i1] = site[ n[i1] ].triangle[7-i1];
  cube_charge(c);
  insert_cube(c);
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static void hypercube_init(){ /* Only in D=4 */
  uchar i0, i1, i2, i3;
  uchar d_num;
  char d[4];
  uint m, n, n2;
  Hypercube *hc;
  Cube *c;
  if( param->D != 4 ) panic("Illegal call");
  for( m = 0; m < param->ipw[param->D]; m++  ){
    hc = &hypercube[m];
    //-------------------------------------------------
    for( i0 = 0; i0 < 4; i0++ ){
#ifdef DEBUG
      if( hc->cube[i0] ) panic("My own cube was already defined by someone else");
      if( current_cube == 4*param->ipw[param->D] ) panic("Wrong cube count");
#endif
      c = hc->cube[i0] = &cube[ current_cube++ ];
      i1 = (i0) ? 0 : 1 ;
      D4dual( i0, i1, &i2, &i3 );
      n = index_up(m, i0);
      d[i0] = 1;
      for( d[i1] = -1; d[i1] < 2; d[i1] += 2 )
	for( d[i2] = -1; d[i2] < 2; d[i2] += 2 )
	  for( d[i3] = -1; d[i3] < 2; d[i3] += 2 ){
	    n2 = n;
	    if( d[i1] > 0 ) n2 = index_up( n2, i1 );
	    if( d[i2] > 0 ) n2 = index_up( n2, i2 );
	    if( d[i3] > 0 ) n2 = index_up( n2, i3 );
	    d_num = four_shifts( -d[0], -d[1], -d[2], -d[3] );
	    c->triangle[three_shifts(d[i1],d[i2],d[i3])] = hypersite[n2].diamond[ d_num ]->triangle[ 2*i0 + 1 ];
	  };
      c->plaq[0] = &plaq[ plaq_index(i3,i2,n) ];              c->plaq_sign[0] = (i3 < i2) ? 1 : -1;
      c->plaq[1] = &plaq[ plaq_index(i2,i3,index_up(n,i1)) ]; c->plaq_sign[1] = (i2 < i3) ? 1 : -1;
      c->plaq[2] = &plaq[ plaq_index(i1,i3,n) ];              c->plaq_sign[2] = (i1 < i3) ? 1 : -1;
      c->plaq[3] = &plaq[ plaq_index(i3,i1,index_up(n,i2)) ]; c->plaq_sign[3] = (i3 < i1) ? 1 : -1;
      c->plaq[4] = &plaq[ plaq_index(i2,i1,n) ];              c->plaq_sign[4] = (i2 < i1) ? 1 : -1;
      c->plaq[5] = &plaq[ plaq_index(i1,i2,index_up(n,i3)) ]; c->plaq_sign[5] = (i1 < i2) ? 1 : -1;
      cube_charge(c);
      insert_cube(c);
#ifdef DEBUG
      if( hypercube[n].cube[i0 + 4] ) panic("This was MY cube");
#endif
      hypercube[n].cube[i0 + 4] = c;
    };
    //-------------------------------------------------
    for( d[0] = -1; d[0] < 2; d[0] += 2 )
      for( d[1] = -1; d[1] < 2; d[1] += 2 )
	for( d[2] = -1; d[2] < 2; d[2] += 2 )
	  for( d[3] = -1; d[3] < 2; d[3] += 2 ){
	    n = m;
	    if( d[0] > 0 ) n = index_up( n, 0 );
	    if( d[1] > 0 ) n = index_up( n, 1 );
	    if( d[2] > 0 ) n = index_up( n, 2 );
	    if( d[3] > 0 ) n = index_up( n, 3 );
	    d_num = four_shifts( -d[0], -d[1], -d[2], -d[3] );
	    i0    = four_shifts( d[0], d[1], d[2], d[3] );
	    hc->diamond[i0] = hypersite[n].diamond[ d_num ];
	  };
    //-------------------------------------------------
    /* Cannot do hypercube insert/charge here */
  };
  for( m = 0; m < param->ipw[param->D]; m++  ){
    hc = &hypercube[m];
    hypercube_charge( hc );
    insert_hypercube( hc );
  };
};
/* ****************************************************************************** */
/* ****************************************************************************** */
/* ****************************************************************************** */
double density_site(){
  uint m, scale = (param->D == 3) ? 1 : 4 ;
  double value = 0.0;
  scale *= param->ipw[param->D];
  for( m = 0; m < scale; m++ ) value += abs( site[m].q );
  value /= (double) scale;
  return value;
};
double density_cube(){
  uint m, scale = (param->D == 3) ? 1 : 4 ;
  double value = 0.0;
  scale *= param->ipw[param->D];
  for( m = 0; m < scale; m++ ) value += abs( cube[m].q );
  value /= (double) scale;
  return value;
};
double density_diamond(){
  uint m, scale = (param->D == 3) ? 0 : 24 ;
  double value = 0.0;
  scale *= param->ipw[param->D];
  for( m = 0; m < scale; m++ ) value += abs( diamond[m].q );
  if( scale ) value /= (double) scale;
  return value;
};
double density_diamond0(){
  uint m, scale = (param->D == 3) ? 0 : 16 ;
  double value = 0.0;
  scale *= param->ipw[param->D];
  for( m = 0; m < scale; m++ )
    if( !(diamond[m].group & 0x01) )
      value += abs( diamond[m].q );
  if( scale ) value /= (double) scale;
  return value;
};
double density_diamond1(){
  uint m, scale = (param->D == 3) ? 0 : 8 ;
  double value = 0.0;
  scale *= param->ipw[param->D];
  for( m = 0; m < scale; m++ )
    if( diamond[m].group & 0x01 )
      value += abs( diamond[m].q );
  if( scale ) value /= (double) scale;
  return value;
};
double density_site_total(){
  uint m, scale, total = 0;
  double value = 0.0;
  //----------------------
  scale  = (param->D == 3) ? 1 : 4 ;
  scale *= param->ipw[param->D];
  total += scale;
  for( m = 0; m < scale; m++ ) value += abs( site[m].q );
  //----------------------
  scale = (param->D == 3) ? 0 : 24 ;
  scale *= param->ipw[param->D];
  total += scale;
  for( m = 0; m < scale; m++ ) value += abs( diamond[m].q );
  //----------------------
  value /= (double) total;
  return value;
};
/* ************************************************************************** */
double *triangle_histogram(){
  uint m;
  uint limit = ( param->D == 3 ) ? 8 : 96 ;
  double *ret = NULL;
  double help;
  if( !triangle ) return NULL;
  MEMORY_CHECK( ret = (double *) calloc( TOTAL_BINS, sizeof(double)), NULL);
  help = 1.0/( (double)(limit * param->ipw[param->D]) );
  for( m = 0 ; m < limit * param->ipw[param->D]; m++ )
    ret[ (int) (0.5 * TOTAL_BINS * (1.0 + triangle[m].phase / M_PI))  ] += help;
  return ret;
};
double *square_histogram(){
  uint m;
  uint limit = ( param->D == 3 ) ? 6 : 24 ;
  double *ret = NULL;
  double help;
  if( !square ) return NULL;
  MEMORY_CHECK( ret = (double *) calloc( TOTAL_BINS, sizeof(double)), NULL);
  help = 1.0/( (double)(limit * param->ipw[param->D]) );
  for( m = 0 ; m < limit * param->ipw[param->D]; m++ )
    ret[ (int) (0.5 * TOTAL_BINS * (1.0 + square[m].phase / M_PI))  ] += help;
  return ret;
};
double *plaq_histogram(){
  uint m;
  uint limit = ( param->D == 3 ) ? 3 : 6 ;
  double *ret = NULL;
  double help;
  if( !plaq ) return NULL;
  MEMORY_CHECK( ret = (double *) calloc( TOTAL_BINS, sizeof(double)), NULL);
  help = 1.0/( (double)(limit * param->ipw[param->D]) );
  for( m = 0 ; m < limit * param->ipw[param->D]; m++ )
    ret[ (int) (0.5 * TOTAL_BINS * (1.0 + plaq[m].phase / M_PI))  ] += help;
  return ret;
};
