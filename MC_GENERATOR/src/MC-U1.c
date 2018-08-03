#include <MC-U1.h>

double complex *F_U1 = NULL;
int *mono = NULL;

/* ***************************************************************************** */
static inline double complex heatbath( double a ){
  double help;
  do{
    help = M_PI * ( 2.0 * RND() - 1.0 );
  }while( cos( help ) < 1.0 + log( RND() )/a );
  return cos(help) + I * sin(help);
};
/* ***************************************************************************** */
static inline double heatbath_one_link( uchar i0, uint m, uchar relax ){
  double complex tmp = 0.0;
  uint m_up, m_left;
  uchar i1;
  double alpha;
  m_up = index_up( m, i0 );
  for( i1 = 0; i1 < param->D ; i1++)
    if( i1 != i0 ){
      m_left = index_down( m, i1 );
      tmp += F_U1[link_index(i1,m_up)] * conj( F_U1[link_index(i1,m)] * F_U1[link_index(i0,index_up(m,i1))] );
      tmp += F_U1[link_index(i1,m_left)] * conj(F_U1[link_index(i0,m_left)] * F_U1[link_index(i1,index_up(m_left,i0))]);
    };
  alpha = cabs(tmp);
  tmp /= alpha;
  m_left = link_index(i0,m);
  if ( relax )   F_U1[m_left] = conj( F_U1[m_left] * tmp );
  else           F_U1[m_left] = heatbath( param->beta * alpha );
  alpha *= creal( F_U1[m_left] );
  F_U1[m_left] *= conj(tmp);
  return alpha;
};
/* ================================ */
double MC_U1(){
  uchar whole, i0;
  uint m;
  double mean = 0.0;
  for( whole = 0 ; whole <= NUMBER_OF_OVERRELAX_U1 ; whole++ )
    for( m = 0; m < param->ipw[param->D]; m++)
      for( i0 = 0; i0 < param->D ; i0++)
	mean += heatbath_one_link( i0, m, whole );
  return mean/( 2 * (NUMBER_OF_OVERRELAX_U1 + 1) * param->D * (param->D-1) * param->ipw[param->D]);
};
/* ================================================================================================ */
double mean_plaq_U1(){
  uint m,m1,l;
  uchar i0, i1;
  double mean = 0.0;
  if( !param || !F_U1 ) panic("Parameters or fields are not inited");
  for( m = 0; m < param->ipw[ param->D ]; m++)
    for( i0 = 0 ; i0 < param->D-1 ; i0++ ){
      l = link_index(i0,m);
      m1 = index_up( m, i0 );
      for( i1 = i0+1 ; i1 < param->D  ; i1++ )
	mean += creal( F_U1[l] * F_U1[link_index(i1,m1)] *
		       conj(F_U1[link_index(i0,index_up(m,i1))] * F_U1[link_index(i1,m)]));
    };
  return 2.0 * mean/( param->D * (param->D-1) * param->ipw[param->D]);
};
/* ================================================================================================ */
void random_gauge_U1(){
  double phase;
  double complex Omega_plus, Omega_minus;
  uint m, m1;
  uchar i0;
  if( !param || !F_U1 )  panic("Parameters or fields are not inited");
  for( m = 0; m < param->ipw[ param->D ]; m++){
    phase = M_PI * ( 2.0 * RND() - 1.0 );
    Omega_plus = cos(phase) + I * sin(phase);
    Omega_minus = conj(Omega_plus);
    for( i0 = 0; i0 < param->D ; i0++){
      m1 = link_index(i0,m);
      F_U1[m1] *= Omega_plus;
      m1 = link_index(i0,index_down(m,i0));
      F_U1[m1] *= Omega_minus;
    };
  };
};
/* ================================================================================================ */
static double optimal_site( uint m ){
  uchar i0;
  double complex optim = 0.0;
  double norm;
  for( i0 = 0; i0 < param->D; i0++ ){
    optim += conj(F_U1[ link_index(i0,m) ]);
    optim += F_U1[ link_index(i0, index_down(m,i0)) ];
  };
  norm = creal(optim) * creal(optim) + cimag(optim) * cimag(optim);
  norm = sqrt(norm);
  optim /= norm;
  for( i0 = 0; i0 < param->D; i0++ ){
    F_U1[ link_index(i0,m) ] *= optim;
    F_U1[ link_index(i0, index_down(m,i0)) ] *= conj(optim);
  };
  return norm;
};
double Landau_gauge_U1(){
  uint m;
  double new = -1.0, old;
  do{
    old = new;
    for( new = 0.0, m = 0; m < param->ipw[param->D]; m++ )
      new += optimal_site( m );
    new /= (double)( 2 * param->D * param->ipw[ param->D ]);
  }while( new > old || old - new > 1.0e-7 );
  return new;
};
/* ================================================================================================ */
static int mono_on_cube( uint m, uchar i0 ){
  uchar i1, i2, i3;
  double phase = 0.0;
  double complex tmp;
  int n;
  if( param->D == 3 ) i0 = 10;
  for( i1 = 0; i1 < param->D; i1++ )
    if( i1 != i0 ){
      if( param->D == 3 ) D3dual( i1, &i2, &i3 );
      else                D4dual( i0, i1, &i2, &i3 );
      tmp = F_U1[link_index(i2,m)] * F_U1[link_index(i3,index_up(m,i2))] *
	conj( F_U1[link_index(i2,index_up(m,i3))] * F_U1[link_index(i3,m)] );
      phase += atan2( cimag(tmp), creal(tmp) );
      n = index_up(m,i1);
      tmp = F_U1[link_index(i2,n)] * F_U1[link_index(i3,index_up(n,i2))] *
	conj( F_U1[link_index(i2,index_up(n,i3))] * F_U1[link_index(i3,n)] );
      phase -= atan2( cimag(tmp), creal(tmp) );
    };
  return (int)rint(phase/(2*M_PI));
};
/* ================================================================================================ */
int * get_mono_U1(){
  int i, m;
  if( !param || !F_U1 )  panic("Parameters or fields are not inited");
  if( !mono ){
    m = param->ipw[ param->D ];
    if( param->D == 4 ) m *= param->D;
    MEMORY_CHECK( mono = (int *) calloc( m, sizeof(int)) , NULL);
  };
  for( m = 0; m < param->ipw[ param->D ]; m++){
    if( param->D == 3 ){
      mono[m] = mono_on_cube( m, 0 );
    }else{
      for( i = 0; i < param->D; i++ )
	mono[ link_index(i, index_down(m,i)) ] = mono_on_cube(m, i);
    };
  };  
  if( param->D == 4 ){
    int q;
    for( m = 0; m < param->ipw[ param->D ]; m++){
      q = 0;
      for( i = 0; i < param->D; i++ ){
	q += mono[ link_index(i,m) ];
	q -= mono[ link_index(i,index_down(m,i)) ];
      };
      if( q ) panic("Current conservation");
    };
  };
  return mono;
};
/* ================================================================================================ */
double mono_density_U1(){
  double ret = 0.0;
  int m, i;
  if( !param || !F_U1 || !mono )  panic("Parameters or fields are not inited");
  for( m = 0; m < param->ipw[ param->D ]; m++){
    if( param->D == 3 ){
      ret += abs( mono[m] );
    }else{
      for( i = 0; i < param->D; i++ )
	ret += abs( mono[link_index(i,m)] );
    };
  };
  if( param->D == 4 ) ret /= (double) param->D;
  return ret/(double) param->ipw[param->D];
};
/* ***************************************************************************** */
/* ***************************************************************************** */
/* ***************************************************************************** */
int fields_init_U1( const char *filename, uchar T, uchar L, uchar D ){
  uint k;
  if( !param ){
    if( param_init(T,L,D) ){
      error("Parameters initialization failed");
      return -1;
    };
    if( !param ) panic("FATAL: parameters still undefined");
  }else{
    if( param->D != D || param->size[0] != T || param->size[1] != L ){
      error("Inconsistent parameters definition");
      return -1;
    };
  };
  if( !F_U1 )
    MEMORY_CHECK( F_U1 = (double complex *)calloc( param->D*param->ipw[param->D], sizeof(double complex)), -1);
  if( !filename ){
    double phase;
    for( k = 0 ; k < param->D*param->ipw[param->D]; k++ ){
      phase = M_PI * (2.0 * RND() - 1.0);
      F_U1[k] = cos(phase) + I * sin(phase);
    };
    return 0;
  };
  return read_filedata_U1( filename );
};
/* ***************************************************************************** */
void fields_free_U1(){
  if( F_U1 ){ free(F_U1); F_U1 = NULL; };
  param_free();
};
/* ***************************************************************************** */
static int fd = -1;
static const char *Name = NULL;

static void cleanup(){
  if( fd >= 0 ) close(fd);
  fd = -1;
  if( F_U1 ) free(F_U1);
  F_U1 = NULL;
};

static int write_something( void *buf, size_t nbytes ){
  int wrote = write( fd, buf, nbytes );  
  if( wrote != nbytes ){
    error("Write failed");
    close(fd); fd = -1;
    if( Name && *Name ){
      unlink( Name );
      Name = NULL;
    };
    return -1;
  };
  return 0;
};

static int driver_read( int fd, void * ptr, size_t nbytes, void (*callback)() ){
  int readbytes = read(fd, ptr, nbytes);
  if( readbytes != nbytes ){
    error("Unexpected EOF");
    (*callback)();
    return -1;
  };
  return 0;
};
/* ***************************************************************************** */
int read_filedata_U1( const char *filename ){
  float value;
  int j, k, idx;
  uchar i0, x[4];
  if( !param || !param->size ){
    error("Illegal parameters given");
    return -1;
  };
  if( !filename || !*filename || (fd = open( filename, O_RDONLY )) < 0 ){
    error("Cannot read open");
    return -1;
  };
#ifdef VERBOSE
  printf("Reading U1 configuration from %s\n", filename );
  printf("Lattice geometry: D = %d, T = %d, L = %d\n", param->D, param->size[0], param->size[1]);;
#endif
  if( !F_U1 && !( F_U1 = (double complex *)calloc( param->D*param->ipw[param->D], sizeof(double complex)) ) ){
    error("Memory allocation failed");
    cleanup();
    return -1;
  };
  /* First four junk bytes */
  if( driver_read( fd, &value, (size_t)4 , &cleanup ) ) return -1;
  for( k = 0; k < param->ipw[param->D] ; k++ )
    for( i0 = 0; i0 < param->D ; i0++ ){
      for( idx = k, j = param->D - 1; j >= 0 ; j-- , idx /= param->size[1] )
      	x[j] = idx % param->size[1];
      idx = param->D * site_index(x) + param->D - 1 - i0;
      if( driver_read(fd, &value, sizeof(float), &cleanup) ) return -1;
      if( fabsf(value) > M_PI  ) error2("WARNING: illegal Abelian gauge field on input: %f", value);
      F_U1[idx] = cos(value) + I * sin(value);
    };
  /* skip either 4 or 16 junk bytes */
  if( driver_read(fd, &value, (size_t) 4, &cleanup) ) return -1;
  if( (k=read(fd, &value, (size_t) 4 )) == 4 ){
    if( driver_read(fd, &value, (size_t) 4, &cleanup) ) return -1;
    if( driver_read(fd, &value, (size_t) 4, &cleanup) ) return -1;
  }else if( k ){
    error("Unexpected EOF");
    cleanup();
    return -1;
  };
  if( read(fd, &value, 1 ) > 0 ){
    error("Pending bytes");
    cleanup();
    return -1;
  };
  close(fd);
  fd = -1;
#ifdef VERBOSE
  printf("Reading U1 configuration finished\n");
#endif
#ifdef DEBUG_U1
  error2("Mean plaquette:  %f\n", mean_plaq_U1() );
#endif
  return 0;
};
/* ***************************************************************************** */
int write_filedata_U1( const char *filename ){
  float value = 0.0;
  int j, k, idx;
  uchar i0, x[4];
  if( !F_U1 ){ error("Fields are not defined"); return -1; };
  if( !param || !param->size ){ error("Illegal parameters given"); return -1; };
  if( !filename || !*filename || (fd = open( filename, O_WRONLY|O_CREAT, 0664)) < 0 ){
    error("Cannot write open"); return -1;
  };
  Name = filename;
#ifdef VERBOSE
  printf("Writing U1 configuration to %s\n", filename );
  printf("Lattice geometry: D = %d, T = %d, L = %d\n", param->D, param->size[0], param->size[1]);
#endif
  if( write_something( &value, (size_t)4 ) ) return -1;
  for( k = 0; k < param->ipw[param->D] ; k++ )
    for( i0 = 0; i0 < param->D ; i0++ ){
      for( idx = k, j = param->D - 1; j >= 0 ; j-- , idx /= param->size[1] )
      	x[j] = idx % param->size[1];
      idx = param->D * site_index(x) + param->D - 1 - i0;
      value = (float) atan2( cimag(F_U1[idx]), creal(F_U1[idx]));
      if( write_something( &value, sizeof(float) ) ) return -1;
    };
  if( write_something(&value, (size_t)4 ) ) return -1;
  close(fd);
  fd = -1;
#ifdef VERBOSE
  printf("Writing U1 configuration finished\n");
#endif
#ifdef DEBUG_U1
  error2("Mean plaquette:  %f\n", mean_plaq_U1() );
#endif
  return 0;
};

