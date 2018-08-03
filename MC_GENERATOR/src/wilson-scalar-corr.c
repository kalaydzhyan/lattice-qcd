#include <wilson-scalar-corr.h>

/* ****************************************************************************************** */
/* *********   Useful defines    ************************************************************ */
/* ****************************************************************************************** */

/* If set, header section in data files is written in separate file
   with additional extension '.header' */
#define SEPARATE_HEADER_FILE  0

/* If set, turns on various time-conserving simplifications (not all loops are considered) */
#define SIMPLIFIED_VERSION    0

/* This is effective ONLY if SIMPLIFIED_VERSION is set. Possible values and meaning:
   0  -  loop time direction is resctricted to be 0;
   1  -  loop time direction is resctricted to be 1 (meaningfull only for symmetric lattices);
   2  -  loop time direction is resctricted to be 2 (meaningfull only for symmetric lattices);
   3  -  loop time direction is resctricted to be 3 (meaningfull only for symmetric lattices); */
#define SIMPLIFIED_VARIANT    0
/* ****************************************************************************************** */

/* Aim is to calculate for given scalar quantity Y(x) its correlation with Wilson loops
  <Y(h,R) W(T,R)> together with <Y> and <W(T,R)>.
  Generic geometry:  T = (2 Lt), R = (2 Ls), Lt is not shown below, R is radial distance.
               /
              /Ls
             /
            * center of the loop
          h/
          /___R____*
         /                                     */
/* ****************************************************************************************** */
/* Global static vars used throughout the file in various functions */

/* Possible values:
   0 - uses lot of memory, fast;
   1 - uses less memory, slower;
   2 - minimal memory, most slow */
uchar wilson_scalar_lowmem = 0;

/* Maximal loop size to consider, temporal and spacial extents.
   Since loops sizes are always even, these are essentially half of the actual lengths */
static uchar  max_loop_size[2] = { 0, 0 };

/* ****************************************************************************************** */
/* ******   Cylindrical coordinates   ******************************************************* */
/* ****************************************************************************************** */
/* Maximal total number of radial points (136 is for 32^4 lattice) */
#define  N_SLICE_POINTS   136

/* Geometry of two-dimensional slice orthogonal to Wilson plane.
   N - number of points with different R; R2 - *squared* distance to origin;
   count - how many points have the same R2; */
typedef struct {
  uint R2[N_SLICE_POINTS], count[N_SLICE_POINTS];
} Slice;
static Slice *slice = NULL;
static uint slice_N = 0;
/* ------------------------------------------------------------------------------------------ */
static int insert_into_slice( uint R2 ){
  uint k;
  for( k = 0; k < slice_N; k++ )
    if( R2 == slice->R2[k] ){ slice->count[k]++; return 0; };
  if( slice_N == N_SLICE_POINTS ) panic("N_SLICE_POINTS is too small");
  slice->R2[slice_N] = R2;
  slice->count[slice_N++] = 1;
  return 0;
};
/* ------------------------------------------------------------------------------------------ */
static void make_slice(){
  uchar x[4];
  uint i,j, tmp;
  if( slice ) return;
  if( !(slice = (Slice *) calloc( 1, sizeof(Slice)) ) ) panic("Memory");
  slice_N = 0;
  bzero( x, 4 * sizeof(uchar) );
  for( x[1] = 0; x[1] < param->size[1]; x[1]++ )
    for( x[2] = 0; x[2] < param->size[2]; x[2]++ )
      if( insert_into_slice( get_distance2( site_index( x ), 0 ) ) )
	panic("Slice insertion failed");
  for( i = 0; i < slice_N - 1; i++ )
    for( j = i+1; j < slice_N; j++ ){
      if( slice->R2[i] == slice->R2[j] ) panic("Doubled entry in slice");
      if( slice->R2[i] > slice->R2[j] ){
	tmp = slice->R2[i];
	slice->R2[i] = slice->R2[j];
	slice->R2[j] = tmp;
	tmp = slice->count[i];
	slice->count[i] = slice->count[j];
	slice->count[j] = tmp;
      };
    };
};
/* ------------------------------------------------------------------------------------------ */
static uint real_get_R2idx( uint start, uint end, uint R2 ){
  uint mid = (uint)(0.5 * ((double)(start + end)) );
  if( R2 == slice->R2[start] ) return start;
  if( R2 == slice->R2[end] ) return end;
  if( start == end ) panic("Zero interval");
  if( R2 >= slice->R2[mid] ) return real_get_R2idx( mid, end, R2 );
  return real_get_R2idx( start, mid, R2 );
};
/* ------------------------------------------------------------------------------------------ */
static uint get_R2idx( uint R2 ){
  return real_get_R2idx( 0, slice_N - 1, R2 );
};
/* ------------------------------------------------------------------------------------------ */
static void free_slice(){
  if( slice ) free(slice );
  slice = NULL;
  slice_N = 0;
};
/* ****************************************************************************************** */
/* ********   Wilson loops construction   *************************************************** */
/* ****************************************************************************************** */
/* Purpose: for given time and spacial directions mu[0], mu[1] and for given temporal and
   spacial extents T[0], T[1] construct all Wilson loops centered at every lattice point. */
/* ------------------------------------------------------------------------------------------ */
/* To speed up Wilson loops construction we first make all possible segments (open Wilson lines)
   of given (even) length.
   !! THIS MIGHT BE MEMORY PROHIBITIVE FOR LARGE LATTICES !! */
static SU2 *segment[4] = { NULL, NULL, NULL, NULL };
/* ------------------------------------------------------------------------------------------ */
/* Segments indexation (separately for time and space dirs). Pointer args unchanged!
   m - segment base point; l - half of segment length (length is always even)  */
static inline uint segment_index( uint m, uchar l ){ return m + param->ipw[param->D] * l; };
/* ------------------------------------------------------------------------------------------ */
/* Segments construction. Pointer args unchanged! G - gauge fields.
   Aimed for unlimited memory ONLY  */
static void make_segments( SU2 *G ){
  uint m, n;
  uchar dir, l;
  SU2 *seg, *seg_next;
  if( wilson_scalar_lowmem ) panic("Low memory - unable to make ALL segments");
  for( dir = 0; dir < 4; dir++ ){
    n = (dir) ? max_loop_size[1] : max_loop_size[0] ;
    n *= param->ipw[param->D];
    if( !segment[dir] && !(segment[dir] = (SU2 *) calloc( n, sizeof(SU2)) ) ) panic("Memory");
    bzero( segment[dir], n * sizeof(SU2) );
  };
#ifdef VERBOSE
  printf("Full segments "); fflush(stdout);
#endif
  for( m = 0; m < param->ipw[param->D]; m++ ){
    for( dir = 0; dir < 4; dir++ ){
      seg = &(segment[dir][segment_index( m, 0 )]);
      *seg = G[ link_index(dir, m) ];
      n = index_up( m, dir );
      *seg = SU2_mult( *seg, G[ link_index( dir, n) ] );
      n = index_up( n, dir );
      for( l = 1; l < max_loop_size[ (dir) ? 1 : 0 ]; l++ ){
	seg_next = &(segment[dir][segment_index( m, l )]);
	*seg_next = SU2_mult( *seg, G[ link_index( dir, n) ] );
	n = index_up( n, dir );
	*seg_next = SU2_mult( *seg_next, G[ link_index( dir, n) ] );
	n = index_up( n, dir );
	seg = seg_next;
      };
    };
#ifdef VERBOSE
    if( !(m % 800) ){ printf("."); fflush(stdout); }
#endif
  };
#ifdef VERBOSE
  printf(" done\n"); fflush(stdout);
#endif
};
/* ------------------------------------------------------------------------------------------ */
static inline uint segment_index_lowmem( uchar *mu, uchar *x, uchar l ){
  return x[mu[0]] + param->size[mu[0]] * ( x[mu[1]] + param->size[mu[1]] * l );
};
/* ------------------------------------------------------------------------------------------ */
/* Segments construction. Pointer args unchanged! G - gauge fields.
   Aimed for limited memory ONLY  */
static void make_segments_lowmem( uchar *mu, uchar *x, SU2 *G ){
  uint m, n;
  uchar dir, l, y[4];
  SU2 *seg, *seg_next;
  if( wilson_scalar_lowmem != 1 ) panic("Low memory - unable to make segments");
  for( dir = 0; dir < 2; dir++ ){
    n = param->size[mu[0]] * param->size[mu[1]] * max_loop_size[dir];
    if( !segment[dir] && !(segment[dir] = (SU2 *) calloc( n, sizeof(SU2)) ) ) panic("Memory");
    bzero( segment[dir], n * sizeof(SU2) );
  };
  bcopy( x, y, 4 * sizeof(uchar) );
  for( y[mu[0]] = 0; y[mu[0]] < param->size[mu[0]]; y[mu[0]]++ ){
    for( y[mu[1]] = 0; y[mu[1]] < param->size[mu[1]]; y[mu[1]]++ ){
      m = site_index( y );
      for( dir = 0; dir < 2; dir++ ){
	seg = &(segment[dir][segment_index_lowmem( mu, y, 0 )]);
	*seg = G[ link_index( mu[dir], m) ];
	n = index_up( m, mu[dir] );
	*seg = SU2_mult( *seg, G[ link_index( mu[dir], n) ] );
	n = index_up( n, mu[dir] );
	for( l = 1; l < max_loop_size[dir]; l++ ){
	  seg_next = &(segment[dir][segment_index_lowmem( mu, y, l )]);
	  *seg_next = SU2_mult( *seg, G[ link_index( mu[dir], n) ] );
	  n = index_up( n, mu[dir] );
	  *seg_next = SU2_mult( *seg_next, G[ link_index( mu[dir], n) ] );
	  n = index_up( n, mu[dir] );
	  seg = seg_next;
	};
      };
    };
  };
};
/* ------------------------------------------------------------------------------------------ */
/* Pair of handy routines (essentially the same as in geometry.c, but might be faster (?)) */
static inline uchar mod_minus(uchar x, uchar dir ){ if(!x ) x = param->size[dir]; return x-1;};
static inline uchar mod_plus( uchar x, uchar dir ){ return (x+1) % param->size[dir]; };
/* ------------------------------------------------------------------------------------------ */
/* From the open lines (segments) constructed above makes the closed Wilson loop.
   mu[0], mu[1] - time/space plane; x - loop center;
   t,r - (half of) temporal/space dimensions of the loop.
   For the case of unrestricted memory ONLY  */
static double get_wilson_loop_value( uchar *mu, uchar *x, uchar t, uchar r ){
  SU2 loop[2];
  uchar i, y0, y1, y[4];
  if( wilson_scalar_lowmem ) panic("Low memory - unable to get loop values");
  y0 = x[mu[0]]; y1 = x[mu[1]];
  for( i = 0; i <= t ; i++ ) y0 = mod_minus( y0, mu[0] );
  for( i = 0; i <= r ; i++ ) y1 = mod_minus( y1, mu[1] );
  bcopy( x, y, 4 * sizeof(uchar));
  y[mu[0]] = y0;   y[mu[1]] = y1;
  loop[0] = segment[mu[0]][segment_index( site_index(y), t )];
  loop[1] = segment[mu[1]][segment_index( site_index(y), r )];
  y[mu[1]] = (y1 + 2* (r + 1)) % param->size[mu[1]];
  loop[1] = SU2_mult( loop[1], segment[mu[0]][segment_index( site_index(y), t )] );
  y[mu[1]] = y1;
  y[mu[0]] = (y0 + 2*(t + 1)) % param->size[mu[0]];
  loop[0] = SU2_mult( loop[0], segment[mu[1]][segment_index( site_index(y), r )] );
  loop[1] = SU2_mult( loop[1], SU2_conj(loop[0]) );
  return creal(loop[1].alpha);
};
/* ------------------------------------------------------------------------------------------ */
/* From the open lines (segments) constructed above makes the closed Wilson loop.
   mu[0], mu[1] - time/space plane; x - loop center;
   t,r - (half of) temporal/space dimensions of the loop
   For restricted memory usage ONLY  */
static double get_wilson_loop_value_lowmem( uchar *mu, uchar *x, uchar t, uchar r ){
  SU2 loop[2];
  uchar i, y0, y1, y[4];
  if( wilson_scalar_lowmem != 1 ) panic("Low memory - unable to get loop values");
  y0 = x[mu[0]]; y1 = x[mu[1]];
  for( i = 0; i <= t ; i++ ) y0 = mod_minus( y0, mu[0] );
  for( i = 0; i <= r ; i++ ) y1 = mod_minus( y1, mu[1] );
  bcopy( x, y, 4 * sizeof(uchar));
  y[mu[0]] = y0;   y[mu[1]] = y1;
  loop[0] = segment[0][segment_index_lowmem( mu, y, t )];
  loop[1] = segment[1][segment_index_lowmem( mu, y, r )];
  y[mu[1]] = (y1 + 2* (r + 1)) % param->size[mu[1]];
  loop[1] = SU2_mult( loop[1], segment[0][segment_index_lowmem( mu, y, t )] );
  y[mu[1]] = y1;
  y[mu[0]] = (y0 + 2*(t + 1)) % param->size[mu[0]];
  loop[0] = SU2_mult( loop[0], segment[1][segment_index_lowmem( mu, y, r )] );
  loop[1] = SU2_mult( loop[1], SU2_conj(loop[0]) );
  return creal(loop[1].alpha);
};
/* ------------------------------------------------------------------------------------------ */
/* Same as above but for lowest memory environment */
static double get_wilson_loop_value_lowestmem( uchar *mu, uchar *x, uchar t, uchar r, SU2 *G ){
  SU2 seg_l, seg_r, seg_b, seg_t;
  uchar i, y_l[4], y_r[4], y_t[4];
  uint n_l, n_r, n_t, n_b;
  bcopy( x, y_l, 4 * sizeof(uchar));
  bcopy( x, y_r, 4 * sizeof(uchar));
  bcopy( x, y_t, 4 * sizeof(uchar));
  for( i = 0; i <= t ; i++ ){
    y_l[mu[0]] = y_r[mu[0]] = mod_minus( y_l[mu[0]], mu[0] );
    y_t[mu[0]] = mod_plus( y_t[mu[0]], mu[0] );
  };
  for( i = 0; i <= r ; i++ ){
    y_l[mu[1]] = y_t[mu[1]] = mod_minus( y_l[mu[1]], mu[1] );
    y_r[mu[1]] = mod_plus( y_r[mu[1]], mu[1] );
  };
  n_l = n_b = site_index( y_l );
  n_t = site_index( y_t );
  n_r = site_index( y_r );
  seg_l.alpha = seg_r.alpha = seg_t.alpha = seg_b.alpha = 1.0;
  seg_l.beta  = seg_r.beta  = seg_t.beta  = seg_b.beta  = 0.0;
  for( i = 0; i < 2*(t+1); i++, n_l = index_up(n_l,mu[0]), n_r = index_up(n_r,mu[0]) ){
    seg_l = SU2_mult( seg_l, G[link_index( mu[0], n_l) ] );
    seg_r = SU2_mult( seg_r, G[link_index( mu[0], n_r) ] );
  };
  for( i = 0; i < 2*(r+1); i++, n_t = index_up(n_t,mu[1]), n_b = index_up(n_b,mu[1]) ){
    seg_t = SU2_mult( seg_t, G[link_index( mu[1], n_t) ] );
    seg_b = SU2_mult( seg_b, G[link_index( mu[1], n_b) ] );
  };
  seg_b = SU2_mult( seg_b , seg_r );
  seg_l = SU2_mult( seg_l , seg_t );
  seg_b = SU2_mult( seg_b , SU2_conj(seg_l) );
  return creal( seg_b.alpha );
};
/* ------------------------------------------------------------------------------------------ */
static void free_loops(){
  uchar i;
  for( i = 0; i < 4; i++ ){
    if( segment[i] ) free( segment[i] );
    segment[i] = NULL;
  };
};
/* ****************************************************************************************** */
/* *********    Fourier  transform  ********************************************************* */
/* ****************************************************************************************** */
static fftwnd_plan fourier_plan = NULL;
static int *fourier_size = NULL;
static fftw_complex *fourier_scalar = NULL; /* fourier image of scalar observables */
static fftw_complex *fourier_loops = NULL;  /* fourier image of loops at fixed plane/geometry */
static fftw_complex *fourier_data = NULL;   /* used to store final correlator in both fourier
					       and eventually x-space */
/* ------------------------------------------------------------------------------------------ */
static void fourier_free(){
  if( fourier_scalar ){ free(fourier_scalar); fourier_scalar = NULL; };
  if( fourier_loops ){ free(fourier_loops); fourier_loops = NULL; };
  if( fourier_data ){ free(fourier_data); fourier_data = NULL; };
  if( fourier_size ){ free(fourier_size); fourier_size = NULL; };
  if( fourier_plan ){ fftwnd_destroy_plan(fourier_plan); fourier_plan = NULL; };
};
/* ------------------------------------------------------------------------------------------ */
/* Allocates the storage and makes fourier transform for scalar observables */
static void fourier_init( uchar N_callback, double (*callback[])( uint ), double *scalar ){
  uchar i;
  uint m, n;
  if( !fourier_scalar
      && !(fourier_scalar = (fftw_complex *) calloc( N_callback * param->ipw[param->D], sizeof(fftw_complex)) ) )
    panic("Memory");
  if( !fourier_loops
      && !(fourier_loops = (fftw_complex *) calloc( param->ipw[param->D], sizeof(fftw_complex)) ) ) panic("Memory");
  if( !fourier_data
      && !(fourier_data = (fftw_complex *) calloc( param->ipw[param->D], sizeof(fftw_complex)) ) ) panic("Memory");
  if( !fourier_size ){
    if( !(fourier_size = (int *) calloc( param->D, sizeof(int))) ) panic("Memory");
    for( i = 0; i < param->D; i++ ) fourier_size[i] = (int) param->size[param->D - i - 1];
  };
  if( !fourier_plan
      && !(fourier_plan = fftwnd_create_plan( param->D, fourier_size, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE )) )
    panic("FFTW plan creation failed");
  bzero( fourier_scalar, N_callback * param->ipw[param->D] * sizeof(fftw_complex) );
#ifdef VERBOSE
  printf("Fourier scalars "); fflush(stdout);
#endif
  for( i = 0; i < N_callback; i++ ){
    n = i * param->ipw[param->D];
    for( m = 0; m < param->ipw[param->D]; m++ )
      fourier_scalar[m + n].re = (*callback[i])(m);
    fftwnd_one( fourier_plan, &(fourier_scalar[n]), NULL );
#ifdef VERBOSE
    printf("."); fflush(stdout);
#endif
    for( m = 0; m < param->ipw[param->D]; m++ ){
      fourier_scalar[m + n].re /= (double) param->ipw[param->D];
      fourier_scalar[m + n].im /= (double) param->ipw[param->D];
    };
    if( scalar ) scalar[i] = fourier_scalar[n].re;
  };
#ifdef VERBOSE
    printf(" done\n"); fflush(stdout);
#endif
};
/* ------------------------------------------------------------------------------------------ */
/* For given loops plane/geometry fourier transform loops array */
static double fourier_make_loops( uchar *mu, uchar *T, SU2 *G ){
  uchar x[4], nu[2];
  uint m;
  if( !fourier_plan || !fourier_loops || !fourier_size ) panic("Fourier is not inited");
  bzero( fourier_loops, param->ipw[param->D] * sizeof(fftw_complex) );
  D4dual( mu[0], mu[1], &nu[0], &nu[1] );
  x[mu[0]] = x[mu[1]] = 0;
  for( x[nu[0]] = 0; x[nu[0]] < param->size[nu[0]]; x[nu[0]]++ )
    for( x[nu[1]] = 0; x[nu[1]] < param->size[nu[1]]; x[nu[1]]++ ){
      if( wilson_scalar_lowmem == 1  && !T[0] && !T[1] )
	make_segments_lowmem( mu, x, G );
      for( x[mu[0]] = 0; x[mu[0]] < param->size[mu[0]]; x[mu[0]]++ )
	for( x[mu[1]] = 0; x[mu[1]] < param->size[mu[1]]; x[mu[1]]++ ){
	  m = site_index(x);
	  if( wilson_scalar_lowmem == 2 ){
	    fourier_loops[m].re = get_wilson_loop_value_lowestmem( mu, x, T[0], T[1], G );
	  }else if( wilson_scalar_lowmem == 1 ){
	    fourier_loops[m].re = get_wilson_loop_value_lowmem( mu, x, T[0], T[1] );
	  }else{
	    fourier_loops[m].re = get_wilson_loop_value( mu, x, T[0], T[1] );
	  };
	};
    };
  fftwnd_one( fourier_plan, fourier_loops, NULL );
#ifdef VERBOSE
  printf("f"); fflush(stdout);
#endif
  for( m = 0; m < param->ipw[param->D]; m++ ){
    fourier_loops[m].re /= (double) param->ipw[param->D];
    fourier_loops[m].im /= (double) param->ipw[param->D];
  };
  return fourier_loops[0].re;
};
/* ------------------------------------------------------------------------------------------ */
/* For given loops plane/geometry construct correlator in fourier space and transforms it back
   to x-space. Fourier images of loops and scalars MUST be stored already */
static void fourier_make_data( uchar cb ){
  uint m, n;
  if( !fourier_plan || !fourier_loops || !fourier_scalar || !fourier_size ) panic("Fourier is not inited");
  for( m = 0; m < param->ipw[param->D]; m++ ){
    n = m + cb * param->ipw[param->D];
    fourier_data[m].re  = fourier_loops[m].re * fourier_scalar[n].re;
    fourier_data[m].re += fourier_loops[m].im * fourier_scalar[n].im;
    fourier_data[m].im  = fourier_loops[m].im * fourier_scalar[n].re;
    fourier_data[m].im -= fourier_loops[m].re * fourier_scalar[n].im;
  };
  fftwnd_one( fourier_plan, fourier_data, NULL );
#ifdef DEBUG
  for( m = 0; m < param->ipw[param->D]; m++ )
    if( fabs(fourier_data[m].im) > FLT_EPSILON )
      panic2("Non zero imaginary part %.10f", fourier_data[m].im );
#endif
};
/* ****************************************************************************************** */
/* ***** Correlation function construction for fixed plane/geometry   *********************** */
/* ****************************************************************************************** */
/* Storage for Wilson loops <--> scalar correlator at given plane/geometry */
static double *wilson_scalar = NULL;
static uint wilson_scalar_count = 0;
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
/* Loops are parameterized by their time/spacial extents  */
static inline uint loop_index( uchar t, uchar r ){ /* Actual loop size is (2 t) X (2 r) */
  return t + max_loop_size[0] * r;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
/* Total number of h-points, it is convenient to set once below */
static uchar H_total = 0;
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
/* R2 - squared distance to the loop center line; h - distance along center line to loop center;
   t,r - Wilson loop geometry  */
static uint wilson_scalar_index( uchar h, uint R2idx, uchar t, uchar r ){
  return h + H_total * ( R2idx + slice_N * loop_index(t,r) );
};
/* ------------------------------------------------------------------------------------------ */
/* Resets the correlation function above and allocates storage if needed */
static void wilson_scalar_reset( uchar N_callback ){
  uint N_total = N_callback * slice_N * H_total * max_loop_size[0] * max_loop_size[1];
  if( !wilson_scalar 
      && !(wilson_scalar = (double *) calloc( N_total, sizeof(double) )) )
    panic("Memory");
  bzero( wilson_scalar, N_total * sizeof(double) );
  wilson_scalar_count = 0;
};
/* ------------------------------------------------------------------------------------------ */
static inline uchar get_h( int i, int L ){  while( i >  L/2 ) i -= L; return (uchar) (abs(i)); };
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
/* Adds together correlator's values keeping track of count number */
static void wilson_scalar_add( uchar *mu, uchar *T , uchar cb  ){
  uchar r, h, nu[2], x[4];
  uint R2, N_total = slice_N * H_total * max_loop_size[0] * max_loop_size[1];
  D4dual( mu[0], mu[1], &nu[0], &nu[1] );
  for( x[mu[0]] = 0; x[mu[0]] < param->size[mu[0]]; x[mu[0]]++ ){
    for( x[mu[1]] = 0; x[mu[1]] < param->size[mu[1]]; x[mu[1]]++ ){
      h = get_h( x[mu[1]], param->size[mu[1]] );
#ifdef DEBUG
      if( h >= H_total ) panic("h is out of range");
#endif
      for( x[nu[0]] = 0; x[nu[0]] < param->size[nu[0]]; x[nu[0]]++ )
	for( x[nu[1]] = 0; x[nu[1]] < param->size[nu[1]]; x[nu[1]]++ ){
	  r = get_h( x[nu[0]], param->size[nu[0]] );
	  R2 = r * r;
	  r = get_h( x[nu[1]], param->size[nu[1]] );
	  R2 += r * r;
	  R2 = get_R2idx( R2 );
	  wilson_scalar[ cb * N_total + wilson_scalar_index( h, R2, T[0], T[1] ) ] +=
	    fourier_data[site_index(x)].re/((double) slice->count[R2] );
	};
    };
    if( !cb && !T[0] && !T[1] ) wilson_scalar_count++;
  };
#ifdef VERBOSE
  if(!cb){ printf("a"); fflush(stdout);};
#endif
};
/* ------------------------------------------------------------------------------------------ */
static void wilson_scalar_average( uchar N_callback ){
  uchar cb, h, t[2];
  uint R2, N_total = slice_N * H_total * max_loop_size[0] * max_loop_size[1];;
  double norm;
  if( !wilson_scalar_count ) panic("Zero count");
  for( cb = 0; cb < N_callback; cb++ )
    for( h = 0; h < H_total; h++ ){
      norm = wilson_scalar_count;
      if( h && (h+1) != H_total ) norm *= 2.0;
      for( R2 = 0; R2 < slice_N; R2++ )
	for( t[0] = 0; t[0] < max_loop_size[0]; t[0]++ )
	  for( t[1] = 0; t[1] < max_loop_size[1]; t[1]++ )
	    wilson_scalar[ cb * N_total + wilson_scalar_index( h, R2, t[0], t[1] ) ] /= norm;
  };
};
/* ------------------------------------------------------------------------------------------ */
static void wilson_scalar_free(){
  if( wilson_scalar ){ free(wilson_scalar);  wilson_scalar = NULL; };
  wilson_scalar_count = 0;
  H_total = 0;
};
/* ****************************************************************************************** */
/* ******  Driver for all the above stuff  ************************************************** */
/* ****************************************************************************************** */
static SU2 *Fuzzy = NULL;     /* Do we need smearing? */
static double *scalar = NULL; /* Mean scalar values */
static double *wilson = NULL; /* Mean loops values */
/* ------------------------------------------------------------------------------------------ */
void make_wilson_scalar( uchar smearing, uchar new_run,
			 double (*callback[])( uint ), char *filename[], uchar N_callback ){
  uchar cb, limit, mu[2], T[2];
  uint n, wilson_count;
  FILE *f;
#if SEPARATE_HEADER_FILE
  char buf[MAX_BUF];
#endif
  /* --- Sanity checks --- */
  if( !param || !F ) panic("Need input data");
  if( param->D != 4 ) panic("Works only in D=4");
  if( param->size[1] != param->size[2] || param->size[1] != param->size[3] )
    panic("Lattice MUST be spatially symmetric");
  if( param->size[0] % 2 || param->size[1] % 2 ) panic("Odd lattice extension");
  if( param->size[0] < 4 || param->size[1] < 4 ) panic("Too small lattice - make no sence");
  for( cb = 0; cb < N_callback; cb++ )
    if( !callback[cb] || !filename[cb] || !*filename[cb] ) panic("Illegal usage");

  /* --- 'Predefined' parameters and actions --- */
  max_loop_size[0] = param->size[0]/4;
  max_loop_size[1] = param->size[1]/4;
  limit = ( param->size[0] == param->size[1] ) ? param->D : 1 ;
  H_total = param->size[1]/2 + 1;
  make_slice();
  wilson_scalar_reset( N_callback );
  if( !scalar && !(scalar=(double *)calloc(N_callback,sizeof(double))) ) panic("Memory");
  if( !wilson &&
      !(wilson = (double *)calloc( max_loop_size[0] * max_loop_size[1],sizeof(double))) ) panic("Memory");
  bzero( wilson, max_loop_size[0] * max_loop_size[1] * sizeof(double) );
  wilson_count = 0;
  if( new_run ){ /* Header/data files are truncated for new run ONLY */
    for( cb = 0; cb < N_callback; cb++ ){
#if SEPARATE_HEADER_FILE
      if( !(f=fopen( filename[cb], "w")) ) panic("Failed to write open");
      fclose(f);
      sprintf( &buf[0], "%s.header", filename[cb] );
      if( !(f=fopen( buf, "w")) ) panic("Failed to write open");
#else
      if( !(f=fopen( filename[cb], "w")) ) panic("Failed to write open");
#endif
      if( fwrite( &max_loop_size[0], sizeof(uchar), 2, f) != 2
	  || fwrite( &H_total, sizeof(uchar), 1, f) != 1
	  || fwrite( &slice_N, sizeof(uint), 1, f) != 1
	  || fwrite( slice->R2, sizeof(uint), slice_N, f) != slice_N ) panic("Failed to write metadata");
      fclose(f);
    };
  };
  if( smearing ){
    if( !Fuzzy 
	&& !(Fuzzy = (SU2 *) calloc( param->D * param->ipw[param->D], sizeof(SU2))) )
      panic("Memory");
  }else  Fuzzy = F;

  /* --------- Let's go --------- */
  fourier_init( N_callback, callback, scalar ); /* Got scalars in fourier space */

  if( !smearing   /* Without smearing segments are taken ONLY once (if memory allows) */
      && !wilson_scalar_lowmem )  make_segments( Fuzzy );

  for( mu[0] = 0; mu[0] < limit; mu[0]++ ){ /* Loops time direction */
#if SIMPLIFIED_VERSION
    if( mu[0] != SIMPLIFIED_VARIANT ){
#ifdef VERBOSE
      printf("SIMPLIFIED_VERSION, VARIANT %d: skiping time dir %d\n",
	     SIMPLIFIED_VARIANT , mu[0]);
      fflush(stdout);
#endif
      continue;
    };
#endif
#ifdef VERBOSE
    printf("Time dir %d:\n", mu[0]); fflush(stdout);
#endif
    if( smearing ){ /* Smearing if needed - and segments afterwards (if memory allows) */
      bcopy( F, Fuzzy, param->D * param->ipw[param->D] * sizeof(SU2) );
#ifdef VERBOSE
      printf("\t Smearing "); fflush(stdout);
#endif
      SU2_spacial_smearing( Fuzzy, mu[0] );
#ifdef VERBOSE
      printf(" done\n"); fflush(stdout);
#endif
      if( !wilson_scalar_lowmem ){
#ifdef VERBOSE
	printf("\t "); fflush(stdout);
#endif
	make_segments( Fuzzy );
      };
    };
    for( mu[1] = 0; mu[1] < param->D; mu[1]++ ){ /* Loops space direction */
      if( mu[1] == mu[0] ) continue;
#ifdef VERBOSE
      printf("\t Space dir %d: ", mu[1]); fflush(stdout);
#endif
      for( T[0] = 0; T[0] < max_loop_size[0]; T[0]++ )  /* Loops temporal extent */
	for( T[1] = 0; T[1] < max_loop_size[1]; T[1]++ ){ /* Loops spacial extent */
	  wilson[loop_index(T[0],T[1])] += /* Got fourier image of loops at fixed geometry */
	    fourier_make_loops( mu, T, Fuzzy);
	  if( !T[0] && !T[1] ) wilson_count++;
	  for( cb = 0; cb < N_callback; cb++ ){ /* Take average for each callback */
	    fourier_make_data( cb ); /* Got x-space correlator for fixed loops geometry */
	    wilson_scalar_add( mu, T , cb  ); /* Took values in cylindrical coordinates */
	  }; /* cb */
	}; /* T[1], T[0] */
#ifdef VERBOSE
      printf("\n"); fflush(stdout);
#endif
    }; /* mu[1] */
  }; /* mu[0] */
  wilson_scalar_average( N_callback ); /* Obtain averaged values ready for dumping */
  for( n = 0; n < max_loop_size[0] * max_loop_size[1]; n++ ) wilson[n] /= (double) wilson_count;

  /* --- Dump what had been obtained --- */
  n = slice_N * H_total * max_loop_size[0] * max_loop_size[1];
  for( cb = 0; cb < N_callback; cb++ ){
    if( !(f=fopen( filename[cb], "a")) ) panic("Failed to append open");
    if( fwrite( &wilson_scalar[cb * n], sizeof(double), n, f) != n
	|| fwrite( wilson, sizeof(double), max_loop_size[0]*max_loop_size[1], f) != max_loop_size[0]*max_loop_size[1]
	|| fwrite( &scalar[cb], sizeof(double), 1, f) != 1 ) panic("Failed to append data");
    fclose(f);
  };
};
/* ------------------------------------------------------------------------------------------ */
void free_wilson_scalar(){
  if( Fuzzy ){ free(Fuzzy); Fuzzy = NULL; };
  if( scalar ){ free(scalar); scalar = NULL; };
  if( wilson ){ free(wilson); wilson = NULL; };
  free_slice();
  free_loops();
  fourier_free();
  wilson_scalar_free();
};
/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */
void mean_wilson_scalar( char *old, char *new ){
  uchar t, r, h;
  uint  i, k, count = 0, WS_total, W_total, *R2 = NULL;
  FILE *f;
  double *mean_ws = NULL, *err_ws = NULL;
  double *mean_w = NULL, *err_w = NULL;
  double *ws = NULL, *w = NULL;
  double mean_s = 0.0, err_s = 0.0, s = 0.0, help, tmp;
  if( !old || !*old || !new || !*new ) panic("Illegal usage");
  if( !(f=fopen(old, "r")) ) panic("Read open");
  if( fread( &max_loop_size[0], sizeof(uchar), 2, f) != 2
      || fread( &H_total, sizeof(uchar), 1, f) != 1
      || fread( &slice_N, sizeof(uint), 1, f) != 1 ) panic("Metadata read");
  if( !(R2 = (uint *) calloc( slice_N, sizeof(uint)) ) ) panic("Memory");
  if( fread( R2, sizeof(uint), slice_N, f) != slice_N ) panic("R2 values read");

  W_total = max_loop_size[0] * max_loop_size[1];
  WS_total = slice_N * H_total * W_total;

  if( !(ws = (double *) calloc( WS_total, sizeof(double)) )) panic("Memory");
  if( !(mean_ws = (double *) calloc( WS_total, sizeof(double)) )) panic("Memory");
  if( !(err_ws = (double *) calloc( WS_total, sizeof(double)) )) panic("Memory");
  if( !(w  = (double *) calloc( W_total, sizeof(double)) )) panic("Memory");
  if( !(mean_w  = (double *) calloc( W_total, sizeof(double)) )) panic("Memory");
  if( !(err_w  = (double *) calloc( W_total, sizeof(double)) )) panic("Memory");
  for( ; ; ){
    k = (uint) fread( ws, sizeof(double), WS_total, f);
    if( !k ) break;
    if( k != WS_total
	|| fread( w, sizeof(double), W_total, f) != W_total
	|| fread( &s, sizeof(double), 1, f) != 1 ) panic("Trancated data");
    for( k = 0; k < WS_total; k++ ) mean_ws[k] += ws[k];
    for( k = 0; k < W_total; k++ )  mean_w[k] += w[k];
    mean_s += s;
    count++;
  };
  if( !count ) panic("No one valid record found");
#ifdef VERBOSE
  printf("Mean values: got %d records\n", count);
  fflush(stdout);
#endif
  for( k = 0; k < WS_total; k++ ) mean_ws[k] /= (double) count;
  for( k = 0; k < W_total; k++ ) mean_w[k] /= (double) count;
  mean_s /= (double) count;
  rewind(f);
  if( fread( &max_loop_size[0], sizeof(uchar), 2, f) != 2
      || fread( &H_total, sizeof(uchar), 1, f) != 1
      || fread( &slice_N, sizeof(uint), 1, f) != 1
      || fread( R2, sizeof(uint), slice_N, f) != slice_N ) panic("Reread metadata");
  for( ; ; ){
    k = (uint) fread( ws, sizeof(double), WS_total, f);
    if( !k ) break;
    if( k != WS_total
	|| fread( w, sizeof(double), W_total, f) != W_total
	|| fread( &s, sizeof(double), 1, f) != 1 ) panic("Trancated data");
    for( k = 0; k < WS_total; k++ ) err_ws[k] += (mean_ws[k] - ws[k]) * (mean_ws[k] - ws[k]);
    for( k = 0; k < W_total; k++ )  err_w[k] += (mean_w[k] - w[k]) * (mean_w[k] - w[k]);
    err_s += (mean_s - s) * (mean_s - s);
  };
  fclose(f);
  for( k = 0; k < WS_total; k++ ) err_ws[k] = sqrt(err_ws[k])/((double) count);
  for( k = 0; k < W_total; k++ ) err_w[k] = sqrt(err_w[k])/((double) count);
  err_s = sqrt(err_s)/((double) count);
  /* --------------------------------------------------- */
  if( !(f=fopen( new, "w")) ) panic("Write open");
  fprintf(f, "# Wilson loops <--> Scalar correlation\n");
  fprintf(f, "# Scalar observable mean value: \t %.8f +/- %.8f\n", mean_s, err_s );
  fprintf(f, "#\n");
  fprintf(f, "# Index 0 - <W(T,R)>\n");
  fprintf(f, "#      Format:  T   R   W(T,R)   W_err(T,R)\n");
  fprintf(f, "#\n");
  fprintf(f, "# Index 1 - correlation functions at various T, R, h, rho2\n");
  fprintf(f, "# h    : displacement from center point along W-plane in W-spatial direction\n");
  fprintf(f, "# rho2 : squared distance to W center axis (cylindrical coordinates)\n");
  fprintf(f, "#\n");
  fprintf(f, "#  WS: normalized dimensionless correlator <W(T,R) s(h,rho)>/(<s> <W(T,R)>)\n");
  fprintf(f, "#   S: scalar correlator <W(T,R) s(h,rho)>/<W(T,R)> - <s>\n");
  fprintf(f, "#\n");
  fprintf(f, "#      Format:  rho2   h   T   R    WS WS_err   S S_err\n");
  fprintf(f, "#\n");
  fprintf(f, "#\n");
  fprintf(f, "# T \t R \t\t W(T,R) W_err(T,R)\n");
  for( t = 0; t < max_loop_size[0]; t++ )
    for( r = 0; r < max_loop_size[1]; r++ ){
      k = loop_index(t,r);
      fprintf(f, "  %d \t %d \t\t %.8f %.8f\n", 2 * (t+1), 2 * (r+1), mean_w[k], err_w[k] );
    };
  fprintf(f, "\n\n");
  fprintf(f, "# rho2 \t h \t T \t R \t\t WS WS_err \t\t S S_err\n");
  for( t = 0; t < max_loop_size[0]; t++ )
    for( r = 0; r < max_loop_size[1]; r++ ){
      i = loop_index(t,r);
      for( h = 0; h < H_total; h++ )
	for( count = 0; count < slice_N; count++ ){
	  k = wilson_scalar_index( h, count, t, r );
	  help = mean_ws[k]/(mean_s * mean_w[i]);
	  fprintf( f, "%d \t %d \t %d \t %d \t\t %.10f ", R2[count], h, 2*(t+1), 2*(r+1), help);
	  help = err_s/mean_s;
	  help *= help;
	  tmp = err_w[i]/mean_w[i];
	  help += tmp * tmp;
	  tmp = err_ws[k]/mean_ws[k];
	  help += tmp * tmp;
	  fprintf(f, "%.10f \t\t ", sqrt(help) * mean_ws[k]/(mean_s * mean_w[i]) );
	  help = mean_ws[k]/mean_w[i] - mean_s;
	  fprintf(f, "%.10f ", help );
	  help = err_ws[k]/mean_ws[k];
	  help *= help;
	  tmp = err_w[i]/mean_w[i];
	  help += tmp * tmp;
	  tmp = mean_ws[k]/mean_w[i];
	  tmp *= tmp;
	  help *= tmp;
	  help += err_s * err_s;
	  fprintf(f, "%.10f\n", sqrt(help) );
	};
    };
  fclose(f);
  /* --------------------------------------------------- */
  slice_N = 0;
  max_loop_size[0] = max_loop_size[1] = 0;
  H_total = 0;
  free( ws );
  free( w );
  free( R2 );
  free( mean_ws );
  free( err_ws );
  free( mean_w );
  free( err_w );
};
