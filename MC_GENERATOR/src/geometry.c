#include <geometry.h>

Param *param = NULL;

/* ================================================================== */
int param_init( uchar T, uchar L, uchar D ){
  if( param ) panic("Parameters already defined");
  //if( D != 3 && D != 4 ) panic("Dimension must be 3 or 4");
  MEMORY_CHECK( param       = (Param *) calloc( 1, sizeof(Param)), -1);
  MEMORY_CHECK( param->size = (uchar *) calloc( D, sizeof(uchar)), -1);
  MEMORY_CHECK( param->ipw  = (uint *)  calloc( D+1, sizeof(uint)),-1);
  param->D = D;
  param->q = 1;
  param->size[0] = T;
  param->ipw[0] = 1;
  for( D = 1; D < param->D ; D++ ){
    param->size[D] = L;
    param->ipw[D] = param->ipw[D-1] * param->size[D-1];
  };
  param->ipw[D] = param->ipw[D-1] * param->size[D-1];
  rnd_init();
  return 0;
};
#ifdef IMPROVEMENT
void param_improvement(){
  if(!param->uzero ) panic("uzero is undefined");
  param->beta1 = 1.;
  param->beta2 = -(1. - 0.4805*log(param->uzero)/3.06839)/20./sqrt(param->uzero);
  param->beta3 = 0.0332*log(param->uzero)/3.06839/sqrt(param->uzero);
};
#endif
void param_free(){
  if( param ){
    free( param->ipw );
    free( param->size );
    free(param);
    param = NULL;
  };
};
/* ================================================================== */
#ifdef DEBUG
void print_link(uint m, uchar i0 ){
  uchar k, x[4];
  site_coordinates( x, m );
  printf("[( ");
  for( k = 0; k < param->D; k++ ) printf("%d ", x[k]);
  printf(") %d ]", i0 );
};
#endif

uint site_index( uchar *x ){
  uint k, ret = 0;
  if( !param ) panic("Parameters undefined");
  for( k = 0; k < param->D; k++ ) ret += (x[k] % param->size[k]) * param->ipw[k];
  return ret;
};
void site_coordinates( uchar *x, uint m ){
  uchar k;
  if( !param ) panic("Parameters undefined");
  for( k = 0; k < param->D; m /= param->size[k++] )
    x[k] = (uchar) ( m % param->size[k] );
};
uint index_up( int m, uchar k ){
  uint tmp;
  if( !param ) panic("Parameters undefined");
  tmp = m / param->ipw[k];
  m += param->ipw[k];
  if( tmp % param->size[k] +1 == param->size[k] )
    m -= param->ipw[k+1];
  return (uint) m;
};
uint index_down( int m, uchar k ){
  uint tmp;
  if( !param ) panic("Parameters undefined");
  tmp = m / param->ipw[k];
  m -= param->ipw[k];
  if( tmp % param->size[k] == 0 )
    m += param->ipw[k+1];
  return (uint)m;
};
uint get_distance2( uint m1, uint m2 ){
  uchar x1[4] = { 0, 0, 0, 0 };
  uchar x2[4] = { 0, 0, 0, 0 };
  uchar k;
  int i, ret = 0;
  site_coordinates( x1, m1 );
  site_coordinates( x2, m2 );
  for( ret = 0, k = 0; k < param->D; k++ ){
    i = x1[k] -  x2[k];
    if( i >  param->size[k]/2 ) i -= param->size[k];
    if( i < -param->size[k]/2 ) i += param->size[k];
    ret += i*i;
  };
  return ret;
};
/* ================================================================== */
/* Parse geometry specification: either TxL^D or L^D, where T,L,D - int */
int parse_geometry( char *pattern, uchar *T, uchar *L, uchar *D ){
  uint t, l, d;
  if( !pattern || !*pattern ) return -1;
  if( sscanf( pattern, "%dx%d^%d", &t, &l, &d ) == 3
      || sscanf( pattern, "%dX%d^%d", &t, &l, &d ) == 3 ){
    *T = (uchar) t;
    *L = (uchar) l;
    *D = (uchar) (d+1);
    return 0;
  };
  if( sscanf( pattern, "%d^%d", &l, &d ) == 2 ){
    *T = (uchar) l;
    *L = (uchar) l;
    *D = (uchar) d;
    return 0;
  };
  return -1;
};
/* ================================================================== */
void abuse( char *module, int line, char *fmt, ... ){
  char err[MAX_BUF];
  va_list ap;
  snprintf(err, MAX_BUF, "Module: %s, line %d:\n\t", module, line);
  if( errno ){
    perror( err );
    errno = 0;
  }else{
    fprintf( stderr, "%s", err );
    va_start( ap, fmt );
    vfprintf( stderr, fmt, ap );
    va_end(ap);
    fprintf(stderr, "\n" );
  };
};
/* -------------------------------------------param->beta2 ---------------------- */
#ifdef SYSTEM_RANDOM
void rnd_init(void){
  struct timeval tm;
  if( gettimeofday( &tm, NULL) ){
    error("gettimeofday() failed, will use time() for seed");    
    tm.tv_usec = (long) time(NULL);
  };
  srandom( tm.tv_usec );
};
#endif
/* ------------------------------------------------------------------------- */
double mod_pi( double x ){
  while( x >  M_PI ) x -= 2 * M_PI;
  while( x < -M_PI ) x += 2 * M_PI;
  return x;
};
/* ------------------------------------------------------------------------- */
double mod_one_half( double x ){
  while( x >  0.5 ) x -= 1.0;
  while( x < -0.5 ) x += 1.0;
  return x;
};
/* ------------------------------------------------------------------------- */
/* Generator of gaussian random numbers distributed with probability 
     exp{ -\alpha ( x - x_0 )^2 } */
static uchar  gauss_flag = 0;
static double gauss_value = 0.0;
double gauss_random( double alpha, double shift ){
  double x1, x2, radius;
  if( !gauss_flag ){
    do{
      x1 = 2.0*RND()-1.0;
      x2 = 2.0*RND()-1.0;
      radius = x1 * x1 + x2 * x2;
    }while( radius > 1.0 );
    radius = sqrt( -log(radius)/(alpha * radius) );
    gauss_value = shift + x1 * radius;
    return shift + x2 * radius;
  };
  gauss_flag = 0;
  return gauss_value;
};
/* ------------------------------------------------------------------------- */
/* Generator of random numbers distributed with probability 
   exp[ -alpha (x - x0)^2 - gamma (1 - cos( q x ) ) ] */
double heatbath_nAHM( double alpha, double x0, double gamma, uchar q ){
  double ret;
  do{
    ret = gauss_random( alpha, x0 );
  }while( gamma * (1-cos( mod_pi(q * ret) )) >  -log(RND()) );
  return ret;
};
/* ------------------------------------------------------------------------- */
/* Generator of gaussian positive randoms with probability exp( -alpha (x-x0)^2 ) */
double positive_gauss_random( double alpha , double x0 ){
  double R, limit;
  do{ R = -x0 * sqrt( alpha/(-log(RND())) ); }while( R > 1 );
  if( R < -1 ) limit = M_PI;
  else         limit = acos(R);
  return x0*(1 - cos(limit * (2*RND()-1))/R);

};
/* Generator of positive random numbers distributed with probability 
   exp[ -lambda (x - 1)^2 - alpha ( sqrt(x) - x0 )^2 -sqrt(x) beta ]
*/
double heatbath_higgs_nAHM( double lambda, double alpha, double x0, double beta ){
  double ret, help, tmp;
  do{
    ret = positive_gauss_random( lambda, 1 );
    tmp = help = sqrt(ret);
    tmp -= x0;
  }while( (alpha * tmp * tmp + help * beta) > -log(RND()) );
  return help;
};
/* =========================================================== */
/*  Generates random number x\in[-1;1] with probability \sqrt{1-x^2}exp[alpha*x]
    Creutz algorithm, which works fast only for not so large \alpha (~< 1.685 )
    P(x) ~ \sqrt{1-x^2} e^{\alpha x} dx ~ \sqrt{ 1- (\frac{ \ln z }{ \alpha })^2 } dz
    x \in [-1;1],   z \in [ e^{-\alpha}; e^{\alpha} ]
    We are generating z and then return \ln z / \alpha:
       x = rnd() -- \in [0;1]
       z = e^{-\alpha} + x ( e^{\alpha} - e^{-\alpha} ) = e^{\alpha} ( x + e^{-2 \alpha}( 1-x ) )
       help = \ln z / \alpha = 1.0 + \ln[x + e^{-2 \alpha}( 1-x )] / \alpha
*/
double Creutz( double alpha ){
  double trial;
  do{
    trial = RND();
    trial += exp(-2.0*alpha) * (1.0-trial);
    trial=1.0+log(trial)/alpha;
  }while( sqrt(1.0-trial*trial) < RND() );
  return trial;
};
  /*
    Generates random number x\in[-1;1] with probability \sqrt{1-x^2}exp[alpha*x]
    Algorithm of Kennedy-Pendleton (PLB 156,393) sutable at weak coupling regime
    ( for \alpha >~ 1.685 )
   */
double Kennedy( double alpha ){
  double help = 0.0 ,rand = 0.0;
  do{
    help=cos(2*M_PI*RND());
    help=(log(RND())+log(RND())*help*help)/alpha;
    rand=RND();
  }while( (2.0+help) < 2.0*rand*rand);
  return 1.0+help;
};
/* =========================================================== */

/* ****************************************************************************** */
/*     D==3            D=4
  i    0 0 1 1 2 2     0 0 0 1 1 1 2 2 2 3 3 3
  j    1 2 0 2 0 1     1 2 3 0 2 3 0 1 3 0 1 2
  ret  0 1 0 2 1 2     0 1 2 0 3 4 1 3 5 2 4 5
  
  Total number of planes :  D*(D-1)/2
*/
static const char planes_D3[3][3] = {
  { -1,  0,  1 }, {  0, -1,  2 }, {  1,  2, -1 }
};
static const char planes_D4[4][4] = {
  { -1, 0, 1, 2 }, { 0, -1, 3, 4 }, { 1, 3, -1, 5 }, { 2, 4, 5, -1 }
};
uchar plane_enum( uchar i, uchar j ){
  char ret = -1;
  if( !param ) panic("NULL parameters");
  if( param->D == 3 ) ret = planes_D3[i][j];
  else                ret = planes_D4[i][j];
  if( ret < 0 ) panic("Illegal call to plane_enum");
  return (uchar) ret;
};
uint plaq_index( uchar i, uchar j, uint m ){
  if( !param ) panic("NULL parameters");
  return plane_enum(i,j) + param->D*(param->D-1)*m/2;
};
uint link_index( uchar i, uint m ){
  if( !param ) panic("NULL parameters");
  return i + param->D * m;
};
/* ****************************************************************************** */
void D4dual( uchar i0, uchar i1, uchar *i2, uchar *i3 ){
  switch( i0 ){
  case 0:
    switch( i1 ){
    case 1: *i2 = 2; *i3 = 3; break;
    case 2: *i2 = 3; *i3 = 1; break;
    case 3: *i2 = 1; *i3 = 2; break;
    default: panic("dual 0");
    };
    break;
  case 1:
    switch( i1 ){
    case 0: *i2 = 3; *i3 = 2; break;
    case 2: *i2 = 0; *i3 = 3; break;
    case 3: *i2 = 2; *i3 = 0; break;
    default: panic("dual 1");
    };
    break;
  case 2:
    switch( i1 ){
    case 0: *i2 = 1; *i3 = 3; break;
    case 1: *i2 = 3; *i3 = 0; break;
    case 3: *i2 = 0; *i3 = 1; break;
    default: panic("dual 2");
    };
    break;
  case 3:
    switch( i1 ){
    case 0: *i2 = 2; *i3 = 1; break;
    case 1: *i2 = 0; *i3 = 2; break;
    case 2: *i2 = 1; *i3 = 0; break;
    default: panic("dual 3");
    };
    break;
  default: panic("dual");
  };
};
/* ****************************************************************************** */
void D3dual( uchar i0, uchar *i1, uchar *i2 ){
  switch( i0 ){
  case 0:   *i1 = 1; *i2 = 2; break;
  case 1:   *i1 = 2; *i2 = 0; break;
  case 2:   *i1 = 0; *i2 = 1; break;
  default: panic("dual");
  };
};
