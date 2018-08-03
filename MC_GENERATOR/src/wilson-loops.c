#include <wilson-loops.h>

/* **************************************************************************** */
/* makes standard smearing treating i0 as 'time' direction */
double smearing_alpha = 0.5;
uint   N_smearing     = 100;

void SU2_spacial_smearing( SU2 *G, uchar i0 ){
  uchar i1, i2;
  uint step, m, m1, m2;
  SU2  g, staple;
  if( !G ) panic("NULL input");
  for( step = 0; step < N_smearing; step++ ){
    for( m = 0; m < param->ipw[param->D] ; m++ )
      for( i1 = 0; i1 < param->D ; i1++){
	if( i1 == i0 ) continue;
	m1 = index_up(m,i1);
        g.alpha = g.beta = 0.0;
	for( i2 = 0; i2 < param->D ; i2++){
	  if( i2 == i1 || i2 == i0 ) continue;
	  m2 = index_down(m,i2);
	  staple = SU2_mult( G[link_index(i2,m)], G[link_index(i1,index_up(m,i2))]);
	  staple = SU2_mult( staple , SU2_conj( G[link_index(i2,m1)] ));
	  g.alpha += staple.alpha;
	  g.beta += staple.beta;
	  staple = SU2_mult( G[link_index(i1,m2)], G[link_index(i2,index_up(m2,i1))] );
	  staple = SU2_mult(SU2_conj( G[link_index(i2,m2)] ), staple );
	  g.alpha += staple.alpha;
	  g.beta += staple.beta;
	};
	m1 = link_index(i1,m);
	G[m1].alpha += smearing_alpha * g.alpha;
	G[m1].beta  += smearing_alpha * g.beta;
        (void) SU2_normalize( &G[m1] );
      };
#ifdef VERBOSE
    if( !(step % 20) ){
      printf("."); fflush( stdout );
    };
#endif
  };
};
/* **************************************************************************** */
/* Keep the last entry zero! */
static double hypercubic_alpha_D3[3] = { 3.0,  1.0,       0.0 };
static double hypercubic_alpha_D4[4] = { 3.0,  2.0,  1.0, 0.0 };

static SU2 link_hypercubic_block( SU2 *G, uint m, uchar i0, uchar no_go ){
  uchar i1, depth = 0;
  SU2 ret, g;
  double *alpha = (param->D == 3) ? &hypercubic_alpha_D3[0] : &hypercubic_alpha_D4[0] ;
  ret.alpha = 0.0;
  ret.beta  = 0.0;
  for( i1 = 0; i1 < param->D; i1++ )
    if( no_go & (1 << i1) ){
      depth++;
    }else{
      no_go |= 1 << i1;
      g = link_hypercubic_block( G, m, i1, no_go);
      g = SU2_mult( g, link_hypercubic_block( G, index_up(m,i1) , i0, no_go) );
      g = SU2_mult( g, SU2_conj( link_hypercubic_block( G, index_up(m,i0) , i1, no_go) ) );
      ret.alpha += g.alpha;
      ret.beta  += g.beta;
      m = index_down(m, i1);
      g = SU2_conj( link_hypercubic_block( G, m, i1, no_go) );
      g = SU2_mult( g, link_hypercubic_block( G, m, i0, no_go) );
      g = SU2_mult( g, link_hypercubic_block( G, index_up(m,i0) , i1, no_go) );
      ret.alpha += g.alpha;
      ret.beta  += g.beta;
      m = index_up(m, i1);
    };
  m = link_index( i0, m);
  ret.alpha = G[m].alpha + alpha[depth - 1] * ret.alpha;
  ret.beta  = G[m].beta  + alpha[depth - 1] * ret.beta;
  SU2_normalize( &ret );
  return ret;
};

static SU2 SU2_hypercubic_block( SU2 *G, uint m, uchar i0 ){
  uchar no_go = 1 << i0;
  return link_hypercubic_block( G, m, i0, no_go );
};
/* **************************************************************************** */
/*  makes standard multihit 'integration' of 'time' link (m,i0). Use ONLY for Wilson action */
static SU2 SU2_multihit( SU2 *G, uint m, uchar i0 ){
  uint m1;
  uchar i1;
  SU2 g, tmp;
  double help;
  g.alpha = g.beta = 0.0;
  for( i1 = 0; i1 < param->D ; i1++){
    if( i1 == i0 ) continue;
    tmp = SU2_mult( G[link_index(i1,m)], G[link_index(i0,index_up(m,i1))] );
    tmp = SU2_mult( tmp, SU2_conj(G[link_index(i1,index_up(m,i0))]) );
    g.alpha += tmp.alpha;
    g.beta += tmp.beta;
    m1 = index_down( m, i1 );
    tmp = SU2_mult( SU2_conj(G[link_index(i1,m1)]), G[link_index(i0,m1)] );
    tmp = SU2_mult( tmp, G[link_index(i1,index_up(m1,i0))] );
    g.alpha += tmp.alpha;
    g.beta += tmp.beta;
  };
  help = param->beta * SU2_normalize( &g );
  help = bessel( help );
  tmp.alpha = help * g.alpha;
  tmp.beta  = help * g.beta;
  return tmp;
};
/* **************************************************************************** */
/* Calculates wilson loops for given gauge configuration F (unchanged) in the plane
   (i0,i1) starting from given point m. Loop size [1:1]...[T0:T1].  Loops are *added*
   to double *loops of size T0*T1 which is parameterized as
      k(t,r) = t + T0*r,   t \in [0:T0-1],  r \in [0:T1-1]
   Actual loop size is (t+1) X (r+1)
   It is assumed that SU2 *Fuzzy contains gauge fields already spacially smeared 
   treating i0 as 'time'. Blocking : 
    1 - use standard multihit for links in i0 direction
    2 - use hypercubic blocking for links in i0 direction
    everything else - no blocking at all
*/
static SU2 *wl_build_g_v = NULL;
static SU2 *wl_build_g_h = NULL;
static uint wl_build_count = 0;

static void wilson_loops_build( uint m,  uchar adjoint,
				uchar i0, uchar T0, uchar i1, uchar T1,
				SU2 *G, double *loops ){
  uchar t, r;
  uint m0, n;
  SU2 tmp;
  double help;
  if( !G ) panic("NULL input");
  if( !wl_build_g_v || !wl_build_g_h || !wl_build_count ){
    if( wl_build_g_v ) free(wl_build_g_v);
    if( wl_build_g_h ) free(wl_build_g_h);
    wl_build_count = 0;
    for( t = 0; t < param->D; t++ )
      if( param->size[t] > wl_build_count ) wl_build_count = param->size[t];
    wl_build_count *= wl_build_count;
    wl_build_g_v = (SU2 *) calloc( wl_build_count, sizeof(SU2) );
    wl_build_g_h = (SU2 *) calloc( wl_build_count, sizeof(SU2) );
    if( !wl_build_g_v || !wl_build_g_h ) panic("Memory");
  };
  bzero( wl_build_g_v, wl_build_count * sizeof(SU2) );
  bzero( wl_build_g_h, wl_build_count * sizeof(SU2) );
  T0++;
  T1++;
  for( m0 = m, t = 0; t < T0 ; t++, m0 = index_up(m0,i0) ){
    wl_build_g_h[t].alpha = 1;
    wl_build_g_h[t].beta = 0;
    for( n = m0, r = 1; r < T1 ; r++, n = index_up(n,i1) )
      wl_build_g_h[ t+T0*r ] = SU2_mult( wl_build_g_h[ t+T0*(r-1) ], G[link_index(i1,n)] );
  };
  for( m0 = m, r = 0; r < T1 ; r++, m0 = index_up(m0,i1) ){
    wl_build_g_v[r].alpha = 1;
    wl_build_g_v[r].beta = 0;
    for( n = m0, t = 1; t < T0 ; t++, n = index_up(n,i0) ){
      wl_build_g_v[ r+T1*t ] = SU2_mult( wl_build_g_v[ r+T1*(t-1) ], G[link_index(i0,n)] );
    };
  };
  for( t = 1; t < T0 ; t++ )
    for( r = 1; r < T1 ; r++ ){
      tmp = SU2_mult( wl_build_g_v[  T1*t], wl_build_g_h[t+T0*r] );
      tmp = SU2_mult( wl_build_g_v[r+T1*t], SU2_conj(tmp));
      tmp = SU2_mult( wl_build_g_h[  T0*r], tmp);
      help = creal(tmp.alpha);
      if( adjoint ){
	help *= 2.0;
	help *= help;
	help -= 1.0;
	help /= 3.0;
      };
      loops[ t-1 + (T0-1)*(r-1) ] += help;
    };
};
/* **************************************************************************** */
static uint  N_loops = 0;
static double *loops = NULL;
static SU2 *Fuzzy = NULL;

void wilson_loops_free(){
  N_loops = 0;
  if( loops ){ free(loops); loops = NULL;  };
};

static double *real_wilson_loops( uchar *T, uchar smearing, uchar blocking, uchar adjoint,
				  char *filename ){
  uint m, count;
  uchar i0, i1;
  FILE *f;
  uchar limit = ( param->size[0] == param->size[1] ) ? param->D : 1 ;
  if( !F ) panic("NULL gauge data");
  if( blocking && blocking != 1 && blocking != 2 ) panic("Illegal BLOCKING parameter");
  if( !T ) panic("Loops geometry NOT given");
  if( N_loops != T[0] * T[1] ) wilson_loops_free();
  if( loops ){
    bzero( &loops[0], N_loops * sizeof(double));
  }else{
    N_loops = T[0] * T[1];
    MEMORY_CHECK(loops=(double *)calloc( N_loops, sizeof(double)),NULL);
  };
  if( (smearing || blocking) && !Fuzzy )
    MEMORY_CHECK( Fuzzy = (SU2 *) calloc( param->D * param->ipw[param->D], sizeof(SU2)) , NULL );
#ifdef VERBOSE
  printf("Wilson loops " ); fflush( stdout );
#endif
  for( count = 0, i0 = 0; i0 < limit; i0++ ){
    if( smearing || blocking ){
      bcopy( F, Fuzzy, param->D * param->ipw[param->D] * sizeof(SU2) );
      if( smearing ){
#ifdef VERBOSE
	printf(" Smearing " ); fflush( stdout );
#endif
	SU2_spacial_smearing( Fuzzy, i0 );
      };
      if( blocking ){
#ifdef VERBOSE
	printf(" Blocking " ); fflush( stdout );
#endif
	for( m = 0; m < param->ipw[param->D]; m++ )
	  Fuzzy[ link_index(i0,m) ] = ( blocking == 1 ) ? SU2_multihit(F, m, i0) : SU2_hypercubic_block(F, m, i0) ;
      };
    }else Fuzzy = F;
    for( m = 0 ; m < param->ipw[param->D]; m += WILSON_LOOPS_VOLUME_SHIFT )
      for( i1 = 0; i1 < param->D; i1++ ){
	if( i1 == i0 ) continue;
	wilson_loops_build( m, adjoint, i0, T[0], i1,T[1], Fuzzy, loops );
	count++;
#ifdef VERBOSE
	if( !(count % 400) ){ printf("."); fflush(stdout); };
#endif
      };
  };
  if( !smearing && !blocking ) Fuzzy = NULL;
  for( m = 0 ; m < N_loops; m++ ) loops[m] /= (double) count;
  if( filename && *filename ){
    if( !(f=fopen( filename, "a")) ) panic("File open");
    if( fwrite( T, sizeof(uchar), 2, f) != 2 || fwrite( &loops[0], sizeof(double), N_loops, f) != N_loops )
      panic("File write");
    fclose(f);
  };
#ifdef VERBOSE
  printf("Done\n"); fflush( stdout );
#endif
  return loops;
};
/* **************************************************************************** */
double *wilson_loops_scan_smearing( char *filename ){
  uchar T[2];
  if( !param ) panic("Geometry undefined");
  T[0] = 1;
  T[1] = param->size[1]/2;
  return real_wilson_loops( T, 1, 0, 0, filename );
};
/* **************************************************************************** */
double *wilson_loops( uchar *size, uchar smearing, uchar blocking, char *filename ){
  uchar T[2];
  if( !param ) panic("Geometry undefined");
  T[0] = param->size[0]/2;
  T[1] = param->size[1]/2;
  if( size ){ size[0] = T[0]; size[1] = T[1]; };
  return real_wilson_loops( T, smearing, blocking, 0, filename );
};
/* **************************************************************************** */
double *adjoint_wilson_loops( uchar *size, uchar smearing, uchar blocking, char *filename ){
  uchar T[2];
  if( !param ) panic("Geometry undefined");
  T[0] = param->size[0]/2;
  T[1] = param->size[1]/2;
  if( size ){ size[0] = T[0]; size[1] = T[1]; };
  return real_wilson_loops( size, smearing, blocking, 1, filename );
};
/* **************************************************************************** */
static uint  N_PLine = 0;
static uint  *PLine_mult = NULL;
static PLine *pline = NULL;

static inline uint PLine_get_distance2( uint m1, uint m2 ){
  uchar x1[4], x2[4], k;
  int i;
  uint ret = 0;
  site_coordinates( x1, param->size[0] * m1 );
  site_coordinates( x2, param->size[0] * m2 );
  for( ret = 0, k = 1; k < param->D; k++ ){
    i = x1[k] -  x2[k];
    if( i >  param->size[k]/2 ) i -= param->size[k];
    if( i < -param->size[k]/2 ) i += param->size[k];
    ret += i*i;
  };
  return ret;
};

static int PLine_get_index( uint R2, uchar mult ){
  int k;
  for( k = 0; k < N_PLine && pline[k].R2 != R2; k++ );
  if( k == N_PLine ){
    if( !mult ) panic("MULTI flag is zero, but I see new element");
    N_PLine++;
    MEMORY_CHECK(pline = (PLine *) realloc( pline, N_PLine * sizeof(PLine)), -1);
    MEMORY_CHECK(PLine_mult = (uint *)realloc(PLine_mult, N_PLine * sizeof(uint)), -1);
    pline[k].R2 = R2;
    pline[k].value = 0.0;
    PLine_mult[k] = 0;
  };
  if( mult ) PLine_mult[k]++ ;
  return k;
};
static void PLine_sort(){
  uint n,m, tmp;
  for( n = 0  ; n < N_PLine - 1 ; n++ )
    for( m = n+1; m < N_PLine     ; m++ )
      if( pline[n].R2 > pline[m].R2 ){
	tmp = pline[n].R2;
	pline[n].R2 = pline[m].R2;
	pline[m].R2 = tmp;
	tmp = PLine_mult[n];
	PLine_mult[n] = PLine_mult[m];
	PLine_mult[m] = tmp;
      };
};

void PLine_correlator_free(){
  N_PLine = 0;
  if( pline ){ free( pline ); pline = NULL; };
  if( PLine_mult ){ free( PLine_mult ); PLine_mult = NULL; };
};

static PLine *real_PLine_correlator( uint *count, double *line, uchar blocking, uchar adjoint, char *filename ){
  uint n, m, V3 = param->ipw[param->D] / param->size[0];
  int m1;
  double *values = NULL, line_value;
  double matrix[9];
  SU2 g;
  FILE *f;
  uchar i0, limit = ( param->size[0] == param->size[1] ) ? param->D : 1 ;
  uchar x[4], x_limit[4] = { 1, 1, 1, 1 };
#ifdef VERBOSE
  printf("PL correlator");
  fflush( stdout );
#endif
  for( i0 = 0 ; i0 < param->D; i0++ ) x_limit[i0] = param->size[i0];
  if( N_PLine ){
    if( !pline || !PLine_mult ) panic("NULL references but not NULL count");
    for( m = 0; m < N_PLine; m++ ) pline[m].value = 0.0;
  }else{
    if( pline || PLine_mult) panic("NULL count but not NULL references");
    for( n = 0; n < V3; n++ )
      for( m = n; m < V3; m++ )
	if( PLine_get_index( PLine_get_distance2(n, m), 1) < 0 ) return NULL;
    if( limit > 1 ) for( n = 0; n < N_PLine; n++ ) PLine_mult[n] *= limit;
    PLine_sort();
  };
  if( count ) *count = N_PLine;
  if( line ) *line = 0.0;
  MEMORY_CHECK( values=(double *)calloc( V3, sizeof(double)), NULL);
  for( i0 = 0 ; i0 < limit; i0++ ){
    for( x[0] = 0; x[0] < x_limit[0]; x[0]++ )
      for( x[1] = 0; x[1] < x_limit[1]; x[1]++ )
	for( x[2] = 0; x[2] < x_limit[2]; x[2]++ )
	  for( x[3] = 0; x[3] < x_limit[3]; x[3]++ ){
	    if( x[i0] ) continue;
	    m = site_index( x );
	    if( blocking == 1 ){
	      g = SU2_multihit( F, m, i0 );
	    }else if( blocking == 2 ){
	      g = SU2_hypercubic_block( F, m, i0 );
	    }else{
	      g = F[link_index(i0,m)];
	    };
	    for( n = index_up(m,i0); n != m ; n = index_up(n,i0) ){
	      if( blocking == 1 ){
		g = SU2_mult( g, SU2_multihit( F, n, i0 ) );
	      }else if( blocking == 2 ){
		g = SU2_mult( g, SU2_hypercubic_block( F, n, i0 ) );
	      }else{
		g = SU2_mult( g, F[link_index(i0,n)] );
	      };
	    };
	    for( m = 0, m1 = 1, n = 0; n < param->D; n++ ){
	      if( n == i0 ) continue;
	      m += x[n] * param->ipw[m1++];
	    };
	    m /= param->size[0];
	    if( adjoint ){
	      SO3_from_SU2( &matrix[0], &g );
	      for( g.alpha = 0.0, adjoint = 1; adjoint < 4; adjoint++ )
		g.alpha += matrix[ MATRIX_ELEMENT( adjoint-1, adjoint-1, 3 )];
	      g.alpha /= 3.0;
	    };
	    values[m] = creal(g.alpha);
	  }; /* all x'es */
    line_value = 0.0;
    for( n = 0; n < V3; n++ ){
#ifdef VERBOSE
      if( !(n % 100) ){ printf("."); fflush(stdout); };
#endif
      line_value += values[n];
      for( m = n; m < V3; m++ ){
	if( (m1 = PLine_get_index( PLine_get_distance2( n, m ), 0 )) < 0 ) return NULL;
	pline[m1].value += values[n]*values[m];
      };
    };

    //if( line ) *line += fabs( line_value ) /( (double) V3 );
    if( line ) *line += line_value /( (double) V3 );

  };/* i0 */
  free( values );
  if( line ) *line /= (double) limit;
  for( n = 0; n < N_PLine; n++ ) pline[n].value /= (double) PLine_mult[n];
#ifdef VERBOSE
  printf("Done\n");
  fflush( stdout );
#endif
  if( filename && *filename ){
    if( !(f=fopen( filename, "a")) ) panic("File open");
    fwrite( &N_PLine, sizeof(uint), 1, f);
    for( n = 0; n < N_PLine; n++ ){
      fwrite( &pline[n].R2, sizeof(uint), 1, f);
      fwrite( &pline[n].value, sizeof(double), 1, f);
    };
    if( fclose(f) ) panic("File close");
  };
  return pline;
};

PLine *PLine_correlator( uint *count, double *line, uchar blocking, char *filename ){
  return real_PLine_correlator( count, line, blocking, 0, filename);
};

PLine *adjoint_PLine_correlator( uint *count, double *line, uchar blocking, char *filename ){
  return real_PLine_correlator( count, line, blocking, 1, filename);
};
