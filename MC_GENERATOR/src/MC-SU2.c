#include <MC-SU2.h>

SU2 *F = NULL;

/* ================================ */
int fields_init(const char *filename, uchar T, uchar L, uchar D)
{
 uint k;
 if(!param)
 {
  if(param_init(T,L,D))
  {
   error("Parameters initialization failed");
   return -1;
  };
  if(!param)
   panic("FATAL: parameters still undefined");
 }
 else
 {
  if( param->D != D || param->size[0] != T || param->size[1] != L )
  {
   error("Inconsistent parameters definition");
   return -1;
  };
 };
 if(!F)
  MEMORY_CHECK(F = (SU2 *)calloc( param->D*param->ipw[param->D], sizeof(SU2)), -1);
 if(!filename)
 {
  for( k = 0 ; k < param->D*param->ipw[param->D]; k++ )
  {
   F[k].alpha = 1.0;
   F[k].beta  = 0.0;
   SU2_normalize(&F[k]);
  };
  return 0;
 };
 return read_SU2_filedata(filename);
};

void fields_free(){
  if( F ){ free(F); F = NULL; };
  param_free();
};
/* ================================ */
static int fd = -1;
static const char *Name = NULL;
static void cleanup(){
  if( fd >= 0 ) close(fd);
  fd = -1;
  if( F ) free(F);
  F = NULL;
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

static int driver_read( int fd, void * ptr, size_t nbytes, void (*callback)() )
{
  int readbytes = read(fd, ptr, nbytes);
  if( readbytes != nbytes ){
    error("Unexpected EOF");
    (*callback)();
    return -1;
  };
  return 0;
};

int read_SU2_filedata(const char *filename)
{
 float value[4];
 int j, k, idx;
 uchar i0, x[4];
 if(!param || !param->size)
 {
  error("Illegal parameters given");
  return -1;
 };
 if( !filename || !*filename || (fd = open( filename, O_RDONLY )) < 0 )
 {
  error("Cannot read open");
  return -1;
 };
#ifdef VERBOSE
 printf("Reading SU2 configuration from %s\n", filename );
 printf("Lattice geometry: D = %d, T = %d, L = %d\n", param->D, param->size[0], param->size[1]);;
#endif
 if(!F && !( F = (SU2 *)calloc( param->D*param->ipw[param->D], sizeof(SU2)) ) )
 {
  error("Memory allocation failed");
  cleanup();
  return -1;
 };
 /* First four junk bytes */
 if(driver_read(fd, value, (size_t)4 , &cleanup ))
  return -1;
 for(k = 0; k < param->ipw[param->D]; k++)
  for(i0 = 0; i0 < param->D; i0++)
  {
   for(idx = k, j = param->D - 1; j >= 0; j--, idx /= param->size[1])
    if(j)
     x[j] = idx % param->size[1];
    else
     x[j] = idx;
    idx = link_index(param->D - 1 - i0, site_index(x));
    if(driver_read(fd, value, 4 * sizeof(float), &cleanup))
     return -1;
    F[idx].alpha = ((double) value[0]) + I*((double) value[3]);
    F[idx].beta = ((double) value[2]) + I*((double) value[1]);
    if(fabs(SU2_normalize( &F[idx] ) - 1.0) > 0.01)
     error("WARNING: Not unitary matrix on input");
  };
  /* skip either 4 or 16 junk bytes */
  if(driver_read(fd, value, (size_t) 4, &cleanup))
   return -1;
  if((k=read(fd, value, (size_t) 4 )) == 4)
  {
   if(driver_read(fd, value, (size_t) 8, &cleanup))
    return -1;
  }
  else
   if(k)
   {
    error("Unexpected EOF");
    cleanup();
    return -1;
   };
  if(read(fd, value, 1) > 0)
  {
   error("Pending bytes");
   cleanup();
   return -1;
  };
  close(fd);
  fd = -1;
#ifdef VERBOSE
  printf("Reading SU2 configuration finished, mean plaq: %f\n", mean_plaq_SU2() );
#endif
  return 0;
};

int write_SU2_filedata( const char *filename ){
  float value[4];
  int j, k, idx;
  uchar i0, x[4];
  if( !F ){ error("Fields are not defined"); return -1; };
  if( !param || !param->size ){ error("Illegal parameters given"); return -1; };
  if( !filename || !*filename || (fd = open( filename, O_WRONLY|O_CREAT, 0664)) < 0 ){
    error("Cannot write open"); return -1;
  };
  Name = filename;
#ifdef VERBOSE
  printf("Writing SU2 configuration to %s\n", filename );
  printf("Lattice geometry: D = %d, T = %d, L = %d\n", param->D, param->size[0], param->size[1]);
#endif
  if( write_something(value, (size_t)4 ) ) return -1;
  for( k = 0; k < param->ipw[param->D] ; k++ )
    for( i0 = 0; i0 < param->D ; i0++ ){
      for( idx = k, j = param->D - 1; j >= 0 ; j-- , idx /= param->size[1] )
        if( j ) x[j] = idx % param->size[1];
    else    x[j] = idx;
      idx = link_index( param->D - 1 - i0, site_index(x) );
      value[0] = (float) creal( F[idx].alpha );
      value[1] = (float) cimag( F[idx].beta );
      value[2] = (float) creal( F[idx].beta );
      value[3] = (float) cimag( F[idx].alpha );
      if( write_something( value, (size_t)(4*sizeof(float)) ) ) return -1;
    };
  if( write_something(value, (size_t)4 ) ) return -1;
  close(fd);
  fd = -1;
#ifdef VERBOSE
  printf("Writing SU2 configuration finished\n");
#endif
#ifdef DEBUG_SU2
  error2("Mean plaquette:  %f\n", mean_plaq_SU2() );
#endif
  return 0;
};
/* ================================ */
static inline double heatbath_one_link( uchar i0, uint m, uchar relax ){
  SU2 Staple, tmp;
  uint m_up, m_left;
  uchar i1;
  double alpha;
  Staple.alpha = Staple.beta = 0.0;
  m_up = index_up( m, i0 );
  /*       m_up
      * --- * --- *
      |     |     |
      |     |i0   |
      * --- * --- *
  m_left    m  i1     */
  for( i1 = 0; i1 < param->D ; i1++)
    if( i1 != i0 ){
      m_left = index_down( m, i1 );
      tmp = SU2_mult( F[link_index(i1,m)], F[link_index(i0,index_up(m,i1))] );
      tmp = SU2_mult( F[ link_index(i1,m_up)], SU2_conj(tmp) );
      Staple.alpha += tmp.alpha;
      Staple.beta += tmp.beta;
      tmp = SU2_mult( F[link_index(i0,m_left)], F[link_index(i1,index_down(m_up,i1))] );
      tmp = SU2_mult( SU2_conj(tmp), F[link_index(i1,m_left)] );
      Staple.alpha += tmp.alpha;
      Staple.beta += tmp.beta;
    }; /* i1 */
  alpha = SU2_normalize( &Staple ) * param->beta;
  m_left = link_index(i0,m);
  if ( relax ){
    F[m_left] = SU2_conj( SU2_mult( Staple, SU2_mult(F[m_left], Staple)));
  }else{
    double trace = (alpha < ALPHA_CROSS ) ? Creutz(alpha) : Kennedy(alpha);
    double help = sqrt( 1.0 - trace * trace );
    double phi = 2.0 * M_PI * RND();
    double cos_theta = 2.0 * RND() - 1.0;
    double sin_theta = sqrt( 1.0 - cos_theta * cos_theta );
    tmp.alpha = trace + I * help * cos_theta;
    tmp.beta = help * sin_theta * cos(phi) + I * help * sin_theta * sin(phi);
    F[m_left] = SU2_mult( tmp, SU2_conj( Staple ) );
  };
  (void) SU2_normalize( &F[m_left] );
  tmp = SU2_mult( F[m_left], Staple );
  return creal( tmp.alpha ) * alpha/param->beta;
};
/* ================================ */
static inline double cool_one_link( uchar i0, uint m ){
  SU2 Staple, tmp;
  uint m_up, m_left;
  uchar i1;
  double alpha;
  Staple.alpha = Staple.beta = 0.0;
  m_up = index_up( m, i0 );
  /*       m_up
      * --- * --- *
      |     |     |
      |     |i0   |
      * --- * --- *
  m_left    m  i1     */
  for( i1 = 0; i1 < param->D ; i1++)
    if( i1 != i0 ){
      m_left = index_down( m, i1 );
      tmp = SU2_mult( F[link_index(i1,m)], F[link_index(i0,index_up(m,i1))] );
      tmp = SU2_mult( F[ link_index(i1,m_up)], SU2_conj(tmp) );
      Staple.alpha += tmp.alpha;
      Staple.beta += tmp.beta;
      tmp = SU2_mult( F[link_index(i0,m_left)], F[link_index(i1,index_down(m_up,i1))] );
      tmp = SU2_mult( SU2_conj(tmp), F[link_index(i1,m_left)] );
      Staple.alpha += tmp.alpha;
      Staple.beta += tmp.beta;
    }; /* i1 */
  alpha = SU2_normalize( &Staple ) * param->beta;
  m_left = link_index(i0,m);
  F[m_left].alpha += COOLING_DELTA * (conj(Staple.alpha) - F[m_left].alpha);
  F[m_left].beta  -= COOLING_DELTA * (     Staple.beta   + F[m_left].beta );
  (void) SU2_normalize( &F[m_left] );
  tmp = SU2_mult( F[m_left], Staple );
  return creal( tmp.alpha ) * alpha/param->beta;
};
/* ================================================================================================ */
double MC_SU2(){
  uchar whole, i0;
  uint m;
  double mean = 0.0;
  for( whole = 0 ; whole <= NUMBER_OF_OVERRELAX ; whole++ )
    for( m = 0; m < param->ipw[param->D]; m++)
      for( i0 = 0; i0 < param->D ; i0++)
    mean += heatbath_one_link( i0, m, whole );
  return mean/( 2 * (NUMBER_OF_OVERRELAX + 1) * param->D * (param->D-1) * param->ipw[param->D]);
};
/* ================================================================================================ */
double cool_SU2(){
  uchar i0;
  uint m;
  double mean = 0.0;
  for( m = 0; m < param->ipw[param->D]; m++)
    for( i0 = 0; i0 < param->D ; i0++)
      mean += cool_one_link( i0, m );
  return mean/( 2 * param->D * (param->D-1) * param->ipw[param->D]);
};
/* ================================================================================================ */
double mean_plaq_SU2(){
  uint m,m1,l;
  uchar i0, i1;
  double mean = 0.0;
  SU2 tmp1, tmp2;
  if( !param || !F ) panic("Parameters or fields are not inited");
  for( m = 0; m < param->ipw[ param->D ]; m++)
    for( i0 = 0 ; i0 < param->D-1 ; i0++ ){
      l = link_index(i0,m);
      m1 = index_up( m, i0 );
      for( i1 = i0+1 ; i1 < param->D  ; i1++ ){
    tmp1 = SU2_mult(F[l], F[link_index(i1,m1)] );
    tmp2 = SU2_mult(F[link_index(i1,m)], F[link_index(i0,index_up(m,i1))] );
    tmp1 = SU2_mult( tmp1, SU2_conj(tmp2) );
    mean += creal(tmp1.alpha);
      };
    };
  return 2.0 * mean/( param->D * (param->D-1) * param->ipw[param->D]);
};

double get_plaq_trace(int x, int mu, int nu)
{
 SU2 tmp1, tmp2;
 int x1 = index_up(x,mu);
 int x2 = index_up(x,nu);
 tmp1 = SU2_mult(F[link_index(mu,x)],F[link_index(nu,x1)]);
 tmp2 = SU2_mult(F[link_index(nu,x)],F[link_index(mu,x2)]);
 tmp1 = SU2_mult(tmp1, SU2_conj(tmp2));
 return 2.0*creal(tmp1.alpha);
};

/* ================================================================================================ */
double Landau_gauge_SU2(){
  SU2 g_left, g_right;
  double error, help;
  uint m, l;
  uchar i0;
#ifdef VERBOSE
  uint count = 0;
#endif
  if( !param || !F )  panic("Parameters or fields are not inited");
#ifdef VERBOSE
  printf("Landau gauge fixing ");
  fflush(stdout);
#endif
  do{
#ifdef VERBOSE
    if( !((count++) % 20) ){
      printf(".");
      fflush(stdout);
    };
#endif
    error = 0.0;
    for( m = 0; m < param->ipw[ param->D ]; m++){
      g_left.alpha = g_left.beta = g_right.alpha = g_right.beta = 0.0;
      for( i0 = 0; i0 < param->D ; i0++){
    l = link_index( i0, m );
    g_left.alpha += conj(F[l].alpha);
    g_left.beta  -= F[l].beta;
    l = link_index( i0, index_down(m,i0) );
    g_left.alpha += F[l].alpha;
    g_left.beta  += F[l].beta;
      };
      (void) SU2_normalize( &g_left );
      g_right = SU2_conj( g_left );
      if( (help=(1.0 - creal(g_left.alpha))) > error ) error = help;
      for( i0 = 0; i0 < param->D ; i0++){
    l = link_index( i0, m );
    F[l] = SU2_mult( g_left, F[l] );
    l = link_index( i0, index_down(m,i0) );
    F[l] = SU2_mult( F[l], g_right );
      };
    };
  }while( error > LANDAU_GAUGE_PRECISION );
  error = 0.0;
  for( m = 0; m < param->D * param->ipw[ param->D ]; m++){
    help = cimag(F[m].alpha);
    error += help * help + cnorm( F[m].beta );
  };
  error /= ((double)(param->D * param->ipw[ param->D ]));
  error *= 4.0;
#ifdef VERBOSE
  printf(" Done, value : %f, iterations: %d\n", error, count);
  fflush(stdout);
#endif
  return error;
};
/* ================================================================================================ */
double MaA_gauge_SU2(){
  SU2 Omega, g, sigma3;
  uint m, l;
  uchar i0;
  double error, help;
#ifdef VERBOSE
  uint count = 0;
#endif
  if( !param || !F )  panic("Parameters or fields are not inited");
  sigma3.alpha = I;
  sigma3.beta  = 0.0;
#ifdef VERBOSE
  printf("MaA gauge fixing ");
  fflush(stdout);
#endif
  do{
    error = 0.0;
#ifdef VERBOSE
    if( !((count++) % 20) ){
      printf(".");
      fflush(stdout);
    };
#endif
    for( m = 0; m < param->ipw[ param->D ]; m++){
      Omega.alpha = Omega.beta = 0.0;
      for( i0 = 0; i0 < param->D ; i0++){
    l = link_index( i0, m );
    g = SU2_mult( F[l], SU2_mult( sigma3, SU2_conj(F[l])) );
    Omega.alpha += g.alpha;
    Omega.beta += g.beta;
    l = link_index( i0, index_down(m,i0) );
    g = SU2_mult( SU2_conj(F[l]), SU2_mult( sigma3, F[l]) );
    Omega.alpha += g.alpha;
    Omega.beta += g.beta;
      };
      SU2_normalize( &Omega );
      Omega.alpha = 1.0 - I * Omega.alpha;
      Omega.beta  = I * Omega.beta;
      SU2_normalize( &Omega );
      if( (help=(1.0 - creal(Omega.alpha))) > error ) error = help;
      for( i0 = 0; i0 < param->D ; i0++){
    l = link_index( i0, m );
    F[l] = SU2_mult( SU2_conj(Omega), F[l]);
    l = link_index( i0, index_down(m,i0) );
    F[l] = SU2_mult( F[l], Omega );
      };
    };
  }while( error > MAA_GAUGE_PRECISION);
  for( help = 0.0, m = 0; m < param->D * param->ipw[ param->D ]; m++){
    g = SU2_mult( F[m], sigma3 );
    g = SU2_mult( g, SU2_conj(F[m]) );
    g = SU2_mult( g, sigma3 );
    help -= creal( g.alpha );
  };
  help /= ((double) (param->D * param->ipw[ param->D ]));
#ifdef VERBOSE
  printf(" Done, value : %f, iterations: %d\n", help, count);
  fflush(stdout);
#endif
  return help;
};
/* ================================================================================================ */
double MaA_center_gauge_SU2(){
  double complex z, zeta;
  SU2 Omega;
  uint m, l;
  uchar i0;
  double error, help;
#ifdef VERBOSE
  uint count = 0;
#endif
  if( !param || !F )  panic("Parameters or fields are not inited");
#ifdef VERBOSE
  printf("MaA center gauge fixing ");
  fflush(stdout);
#endif
  do{
    error = 0.0;
#ifdef VERBOSE
    if( !((count++) % 20) ){ printf("."); fflush(stdout); };
#endif
    for( m = 0; m < param->ipw[ param->D ]; m++){
      z = 0.0;
      for( i0 = 0; i0 < param->D ; i0++){
    zeta = F[link_index( i0, m )].alpha;
    zeta /= cabs(zeta);
    z += zeta * zeta;
    zeta = conj(F[link_index( i0, index_down(m,i0) )].alpha);
    zeta /= cabs(zeta);
    z += zeta * zeta;
      };
      help = 0.5 * atan2( cimag(z), creal(z));
      z = cos( help ) + I * sin( help );
      Omega.alpha = z;
      Omega.beta  = 0.0;
      if( (help=(1.0 - creal(z))) > error ) error = help;
      for( i0 = 0; i0 < param->D ; i0++){
    l = link_index( i0, m );
    F[l] = SU2_mult( SU2_conj(Omega), F[l]);
    l = link_index( i0, index_down(m,i0) );
    F[l] = SU2_mult( F[l], Omega );
      };
    };
  }while( error > MAA_CENTER_GAUGE_PRECISION );
  for( help = 0.0, m = 0; m < param->D * param->ipw[ param->D ]; m++){
    zeta = F[m].alpha;
    zeta /= cabs(zeta);
    help += creal(zeta) * creal(zeta);
  };
  help /= ((double) (param->D * param->ipw[ param->D ]));
#ifdef VERBOSE
  printf(" Done, value : %f, iterations: %d\n", help, count );
  fflush(stdout);
#endif
  return help;
};
/* ================================================================================================ */
void random_gauge_SU2(){
  SU2 Omega;
  uint m, m1;
  uchar i0;
  if( !param || !F )  panic("Parameters or fields are not inited");
  for( m = 0; m < param->ipw[ param->D ]; m++){
    Omega.alpha = RND() + I * RND();
    Omega.beta = RND() + I * RND();
    (void) SU2_normalize( &Omega );
    for( i0 = 0; i0 < param->D ; i0++){
      m1 = link_index(i0,m);
      F[m1] = SU2_mult( SU2_conj(Omega), F[m1] );
      m1 = link_index(i0,index_down(m,i0));
      F[m1] = SU2_mult( F[m1], Omega );
    };
  };
};
