#include <MC-SU_N.h>
#include <stdio.h>
SU_N *F_N = NULL;

/* ******************************************************************* */
/* ****   SU(N)/U(N) init and reading/writing   **************************** */
/* ******************************************************************* */
#define U_N_DRIVERS  1
#define SU_N_DRIVERS 3

static int U_N_read(  const char *filename, uchar T, uchar L, uchar D );
static int U_N_write( const char *filename );
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int SU2_read( const char *filename );
static int SU2_write( const char *filename );
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int SU3_bqcd_read( const char *filename );
static int SU3_bqcd_write( const char *filename );
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int SU_N_read( const char *filename );
static int SU_N_write( const char *filename );
/* ----------------------------------------------------------- */
static const char * U_N_read_write_driver_desc[U_N_DRIVERS] = {
  "U(N) data for any N and D=3,4. Links/site enumeration is the same as in SU3 bqcd format.\n"
  "Header goes first:\n"
  "   float beta, uchar N, uchar T, uchar L, uchar D\n"
  "Record for each link U[i][j]\n"
  "   re(U[0][0]),im(U[0][0]),re(U[0][1]),im(U[0][1]), ..., re(U[0][N-1]),im(U[0][N-1])\n"
  "   re(U[1][0]),im(U[1][0]),re(U[1][1]),im(U[1][1]), ..., re(U[1][N-1]),im(U[1][N-1])\n"
  "      ...         ...         ...         ...       ...     ...        ...\n"
  "   re(U[N-1][0]),im(U[N-1][0]),                     ..., re(U[N-1][N-1]),im(U[N-1][N-1])\n"
  "where each particular chunk is in 'float' little endian\n"
};

static const char * SU_N_read_write_driver_desc[SU_N_DRIVERS] = {
  "SU2 data in V.Bornyakov format, only D=3,4 are supported\n"
  "(details could be found elsewhere)\n",

  "SU3 data in BQCD format without time slicing, only D=3,4 are supported\n"
  "(details could be found elsewhere)\n",

  "SU(N) data for any N and D=3,4. For (N=2 || N=3) is the same as the above drivers.\n"
  "For N > 3: links/site enumeration is the same as in SU3 bqcd format.\n"
  "Record for each link U[i][j]\n"
  "   re(U[0][0]),im(U[0][0]),re(U[0][1]),im(U[0][1]), ..., re(U[0][N-1]),im(U[0][N-1])\n"
  "   re(U[1][0]),im(U[1][0]),re(U[1][1]),im(U[1][1]), ..., re(U[1][N-1]),im(U[1][N-1])\n"
  "      ...         ...         ...         ...       ...     ...        ...\n"
  "   re(U[N-1][0]),im(U[N-1][0]),                     ..., re(U[N-1][N-1]),im(U[N-1][N-1])\n"
  "where each particular chunk is in 'float' little endian\n"
};
/* --------------------------------------------------------------------- */
void SU_N_list_drivers(){
  uchar i;
  printf("Available read/write data formats for SU(N) gauge fields:\n");
  for( i = 0; i < SU_N_DRIVERS; i++ ){
    printf("---------------- Driver %d -----------------------------------------\n", i );
    printf("%s", SU_N_read_write_driver_desc[i] );
  };
  fflush(stdout);
};

/* --------------------------------------------------------------------- */
void U_N_list_drivers(){
  uchar i;
  printf("Available read/write data formats for U(N) gauge fields:\n");
  for( i = 0; i < U_N_DRIVERS; i++ ){
    printf("---------------- Driver %d -----------------------------------------\n", i );
    printf("%s", U_N_read_write_driver_desc[i] );
  };
  fflush(stdout);
};
/* --------------------------------------------------------------------- */
int SU_N_init( const char *filename, uchar T, uchar L, uchar D ){
  uint k;
  if( D != 3 && D != 4 ) panic("I'm too lazy, only D=3,4 is implemented");
  if( !param ){
    if( param_init(T,L,D) ) panic("Parameters initialization failed");
  }else{
    if( param->D != D || param->size[0] != T || param->size[1] != L )
      panic("Floating geomety?!");
  };
  if( !F_N )
    MEMORY_CHECK( F_N = (SU_N *)calloc( param->D*param->ipw[param->D], sizeof(SU_N)), -1);
  if( !filename ){
    for( k = 0 ; k < param->D*param->ipw[param->D]; k++ )
      SU_N_random( &F_N[k] );
    return 0;
  };
  if( SU_N_RANK == 2 ) return SU2_read( filename );
  if( SU_N_RANK == 3 ) return SU3_bqcd_read( filename );
  return SU_N_read( filename );
};
/* --------------------------------------------------------------------- */
int U_N_init( const char *filename, uchar T, uchar L, uchar D ){
  uint k;
  if( !filename ){
    if( D != 3 && D != 4 ) panic("I'm too lazy, only D=3,4 is implemented");
    if( !param ){
      if( param_init(T,L,D) ) panic("Parameters initialization failed");
    }else{
      if( param->D != D || param->size[0] != T || param->size[1] != L )
    panic("Floating geomety?!");
    };
    if( !F_N )
      MEMORY_CHECK( F_N = (SU_N *)calloc( param->D*param->ipw[param->D], sizeof(SU_N)), -1);
    for( k = 0 ; k < param->D*param->ipw[param->D]; k++ ) SU_N_random_U_N( &F_N[k] );
    return 0;
  };
  return U_N_read( filename , T, L, D );
};
/* --------------------------------------------------------------------- */
void SU_N_free(){
  if( F_N ){ free(F_N); F_N = NULL; };
  param_free();
};
/* --------------------------------------------------------------------- */
void U_N_free(){
  SU_N_free();
};
/* --------------------------------------------------------------------- */
int SU_N_dump( const char *filename ){
  if( !param || !F_N ) panic("Geometry/gauge data undefined");
  if( SU_N_RANK == 2 ) return SU2_write( filename );
  if( SU_N_RANK == 3 ) return SU3_bqcd_write( filename );
  return SU_N_write( filename );
};
/* --------------------------------------------------------------------- */
int  U_N_dump( const char *filename ){
  if( !param || !F_N ) panic("Geometry/gauge data undefined");
  return U_N_write( filename );
};
/* --------------------------------------------------------------------- */
static void make_SU_N( SU_N *A, float *buf ){
  uchar i, j;
  uint n;
  for( i = 0; i < SU_N_RANK; i++ )
    for( j = 0; j < SU_N_RANK; j++ ){
      n = 2 * (SU_N_RANK * i + j);
      A->U[j][i] = ((double)buf[n]) + I * ((double)buf[n + 1]);
    };
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static void make_buf( SU_N *A, float *buf ){
  uchar i, j;
  uint n;
  for( i = 0; i < SU_N_RANK; i++ )
    for( j = 0; j < SU_N_RANK; j++ ){
      n = 2 * ( SU_N_RANK * i + j );
      buf[n]   = (float) creal(A->U[j][i]);
      buf[n+1] = (float) cimag(A->U[j][i]);
    };
};
/* --------------------------------------------------------------------- */
static int U_N_read( const char *filename, uchar T, uchar L, uchar D ){
  FILE *f;
  uchar i0, x[4], lim3 = (param->D == 4) ? param->size[3] : 1 ;
  uint N;
  float buffer[2 * SU_N_RANK * SU_N_RANK];
  if( !(f= fopen( filename, "r")) ) panic("Cannot read open");
  if( fseek(f, 0, SEEK_END) ) panic("fseek to END failed");
  N = 2 * SU_N_RANK * SU_N_RANK;
  N *= param->ipw[param->D] * param->D * sizeof(float);
  N +=  sizeof(float) + 4*sizeof(uchar);
  if( ftell(f) != N ) panic("Filesize is NOT as expected");
  if( fseek(f, 0, SEEK_SET) ) panic("fseek to START failed");
  if( fread(&buffer[0], sizeof(float), 1, f) != 1 ) panic("Beta read");
  if( param && fabs(buffer[0] - param->beta) > 0.0001 ) panic("Data is for different beta value");
  if( fread( &x[0], sizeof(uchar), 4, f) != 4 ) panic("Header read");
  if( x[0] != SU_N_RANK ) panic3("Incompatible data, found N_c=%d, compilled for %d", x[0], SU_N_RANK);
  if( param ){
    if( x[1] !=  param->size[0] ) panic3("Data is for different T, found %d, have %d", x[1], param->size[0]);
    if( x[2] !=  param->size[1] ) panic3("Data is for different L, found %d, have %d", x[2], param->size[1]);
    if( x[3] !=  param->D ) panic3("Data is for different D, found %d, have %d", x[3], param->D);
  }else{
    if( x[3] != 3 && x[3] != 4 ) panic2("Only D=3,4 implemented (got %d)", x[3] );
    U_N_free();
    if( param_init( x[1], x[2], x[3] ) ) panic("Parameters initialization failed");
    MEMORY_CHECK( F_N = (SU_N *)calloc( param->D*param->ipw[param->D], sizeof(SU_N)), -1);
    param->beta = buffer[0];
  };
#ifdef VERBOSE
  printf("[geometry %dx%d^%d] : reading U(%d) configuration from %s\n",
     param->size[0], param->size[1], param->D-1, SU_N_RANK, filename );
  fflush(stdout);
#endif
  N = 2 * SU_N_RANK * SU_N_RANK;
  for( x[0] = 0; x[0] < param->size[0] ;  x[0]++ )
    for( x[3] = 0; x[3] < lim3;  x[3]++ )
      for( x[2] = 0; x[2] < param->size[2] ;  x[2]++ )
    for( x[1] = 0; x[1] < param->size[1] ;  x[1]++ )
      for( i0 = 0; i0 < param->D; i0++ ){
        if( fread( &buffer[0], sizeof(float), N, f) != N ){
          fclose(f);
          panic("Failed to read chunk");
        };
        make_SU_N( &F_N[link_index( (i0 + 1) % param->D, site_index(x))], buffer);
      };
  fclose(f);
#ifdef VERBOSE
  printf("[geometry %dx%d^%d] : U(%d) gauge data read, mean plaq. %f\n",
     param->size[0], param->size[1], param->D-1, SU_N_RANK, SU_N_mean_plaq() );
  fflush(stdout);
#endif
  return 0;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int U_N_write( const char *filename ){
  FILE *f;
  uchar i0, x[4], lim3 = (param->D == 4) ? param->size[3] : 1 ;
  uint N = 2 * SU_N_RANK * SU_N_RANK;
  float buffer[2 * SU_N_RANK * SU_N_RANK];
  if( !(f= fopen( filename, "w")) ) panic("Cannot write open");
  buffer[0] = (float) param->beta;
  x[0] = SU_N_RANK;
  x[1] = param->size[0];
  x[2] = param->size[1];
  x[3] = param->D;
  if( fwrite( &buffer[0], sizeof(float), 1, f) != 1
      || fwrite( &x[0], sizeof(uchar), 4, f) != 4 ) panic("Header write");
  for( x[0] = 0; x[0] < param->size[0] ;  x[0]++ )
    for( x[3] = 0; x[3] < lim3;  x[3]++ )
      for( x[2] = 0; x[2] < param->size[2] ;  x[2]++ )
    for( x[1] = 0; x[1] < param->size[1] ;  x[1]++ )
      for( i0 = 0; i0 < param->D; i0++ ){
        make_buf( &F_N[link_index( (i0 + 1) % param->D, site_index(x))], buffer);
        if( fwrite( &buffer[0], sizeof(float), N, f) != N ){
          fclose(f);
          panic("Failed to dump gauge matrix");
        };
      };
  fclose(f);
  return 0;
};
/* --------------------------------------------------------------------- */
static int SU2_read( const char *filename ){
  int m, flag = (F == NULL) ? 1 : 0;
  if( (m = read_SU2_filedata( filename )) ) return m;
  if( !F_N )
    MEMORY_CHECK( F_N = (SU_N *)calloc( param->D*param->ipw[param->D], sizeof(SU_N)), -1);
  for( m = 0; m < param->D  * param->ipw[param->D]; m++ ){
    F_N[m].U[0][0] = F[m].alpha;
    F_N[m].U[0][1] = F[m].beta;
    F_N[m].U[1][1] = conj(F[m].alpha);
    F_N[m].U[1][0] = -conj(F[m].beta);
  };
  if( flag ){ free( F ); F = NULL; };
  return 0;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int SU2_write( const char *filename ){
  int m, flag = (F == NULL) ? 1 : 0, ret = 0;
  if( !F )
    MEMORY_CHECK( F = (SU2 *)calloc( param->D*param->ipw[param->D], sizeof(SU2)), -1);
  for( m = 0; m < param->D  * param->ipw[param->D]; m++ ){
    F[m].alpha = F_N[m].U[0][0];
    F[m].beta  = F_N[m].U[0][1];
  };
  ret = write_SU2_filedata( filename );
  if( flag ){ free( F ); F = NULL; };
  return ret;
};
/* --------------------------------------------------------------------- */
static int make_bqcd_SU_N( SU_N *A, float *buf ){
  uchar i, j;
  const uchar m[4] = {1, 2, 0, 1};
  uint n;
  for( i = 0; i < 2; i++ )
    for( j = 0; j < 3; j++ ){
      n = 2 * (j + 3*i);
      A->U[j][i] = ((double)buf[n]) + I * ((double)buf[n + 1]);
    };
  /* 3-d column restoration: A = || a b c || , c = conj(a x b),
     where a,b,c are 3-vectors and "x" is vector product */
  for( i = 0; i < 3; i++ )
    A->U[i][2] = conj(A->U[m[i]][0] * A->U[m[i+1]][1] - A->U[m[i+1]][0] * A->U[m[i]][1]);
  return 0;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int SU3_bqcd_read( const char *filename ){
  FILE *f;
  uchar i0, x[4], lim3 = (param->D == 4) ? param->size[3] : 1 ;
  float buffer[12];
  if( !(f= fopen( filename, "r")) ) panic("Cannot read open");
  if( fseek(f, 0, SEEK_END) ) panic("fseek to END failed");
  if( ftell(f) != param->ipw[param->D] * param->D * 12 * sizeof(float) )
    panic("Filesize is NOT as expected");
  if( fseek(f, 0, SEEK_SET) ) panic("fseek to START failed");
#ifdef VERBOSE
  printf("[geometry %dx%d^%d] : reading SU(3) configuration from %s\n",
     param->size[0], param->size[1], param->D-1, filename );
  fflush(stdout);
#endif
  if( !F_N )
    MEMORY_CHECK( F_N = (SU_N *)calloc( param->D*param->ipw[param->D], sizeof(SU_N)), -1);
  for( x[0] = 0; x[0] < param->size[0] ;  x[0]++ )
    for( x[3] = 0; x[3] < lim3; x[3]++ )
      for( x[2] = 0; x[2] < param->size[2] ;  x[2]++ )
    for( x[1] = 0; x[1] < param->size[1] ;  x[1]++ )
      for( i0 = 0; i0 < param->D; i0++ )
        if( fread( &buffer[0], sizeof(float), 12, f) != 12
        || make_bqcd_SU_N( &F_N[link_index( (i0 + 1) % param->D, site_index(x))], buffer) ){
          fclose(f);
          panic("Failed to reconstruct gauge matrix");
        };
  fclose(f);
#ifdef VERBOSE
  printf("[geometry %dx%d^%d] : SU(3) gauge data read, mean plaq. %f\n",
     param->size[0], param->size[1], param->D-1, SU_N_mean_plaq() );
  fflush(stdout);
#endif
  return 0;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static void make_bqcd_buf( SU_N *A, float *buf ){
  uchar i, j;
  uint n;
  for( i = 0; i < 2; i++ )
    for( j = 0; j < 3; j++ ){
      n = 2 * (j + 3*i);
      buf[n]   = (float) creal(A->U[j][i]);
      buf[n+1] = (float) cimag(A->U[j][i]);
    };
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int SU3_bqcd_write( const char *filename ){
  FILE *f;
  uchar i0, x[4], lim3 = (param->D == 4) ? param->size[3] : 1 ;
  float buffer[12];
  if( !(f= fopen( filename, "w")) ) panic("Cannot write open");
  for( x[0] = 0; x[0] < param->size[0] ;  x[0]++ )
    for( x[3] = 0; x[3] < lim3;  x[3]++ )
      for( x[2] = 0; x[2] < param->size[2] ;  x[2]++ )
    for( x[1] = 0; x[1] < param->size[1] ;  x[1]++ )
      for( i0 = 0; i0 < param->D; i0++ ){
        make_bqcd_buf( &F_N[link_index( (i0 + 1) % param->D, site_index(x))], buffer);
        if( fwrite( &buffer[0], sizeof(float), 12, f) != 12 ){
          fclose(f);
          panic("Failed to dump gauge matrix");
        };
      };
  fclose(f);
  return 0;
};
/* --------------------------------------------------------------------- */
static int SU_N_read( const char *filename ){
  FILE *f;
  uchar i0, x[4], lim3 = (param->D == 4) ? param->size[3] : 1 ;
  uint N = 2 * SU_N_RANK * SU_N_RANK;
  float buffer[2 * SU_N_RANK * SU_N_RANK];
  if( !(f= fopen( filename, "r")) ) panic("Cannot read open");
  if( fseek(f, 0, SEEK_END) ) panic("fseek to END failed");
  if( ftell(f) != param->ipw[param->D] * param->D * N * sizeof(float) )
    panic("Filesize is NOT as expected");
  if( fseek(f, 0, SEEK_SET) ) panic("fseek to START failed");
#ifdef VERBOSE
  printf("[geometry %dx%d^%d] : reading SU(%d) configuration from %s\n",
     param->size[0], param->size[1], param->D-1, SU_N_RANK, filename );
  fflush(stdout);
#endif
  if( !F_N )
    MEMORY_CHECK( F_N = (SU_N *)calloc( param->D*param->ipw[param->D], sizeof(SU_N)), -1);
  for( x[0] = 0; x[0] < param->size[0] ;  x[0]++ )
    for( x[3] = 0; x[3] < lim3;  x[3]++ )
      for( x[2] = 0; x[2] < param->size[2] ;  x[2]++ )
    for( x[1] = 0; x[1] < param->size[1] ;  x[1]++ )
      for( i0 = 0; i0 < param->D; i0++ ){
        if( fread( &buffer[0], sizeof(float), N, f) != N ){
          fclose(f);
          panic("Failed to read chunk");
        };
        make_SU_N( &F_N[link_index( (i0 + 1) % param->D, site_index(x))], buffer);
      };
  fclose(f);
#ifdef VERBOSE
  printf("[geometry %dx%d^%d] : SU(%d) gauge data read, mean plaq. %f\n",
     param->size[0], param->size[1], param->D-1, SU_N_RANK, SU_N_mean_plaq() );
  fflush(stdout);
#endif
  return 0;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int SU_N_write( const char *filename ){
  FILE *f;
  uchar i0, x[4], lim3 = (param->D == 4) ? param->size[3] : 1 ;
  uint N = 2 * SU_N_RANK * SU_N_RANK;
  float buffer[2 * SU_N_RANK * SU_N_RANK];
  if( !(f= fopen( filename, "w")) ) panic("Cannot write open");
  for( x[0] = 0; x[0] < param->size[0] ;  x[0]++ )
    for( x[3] = 0; x[3] < lim3;  x[3]++ )
      for( x[2] = 0; x[2] < param->size[2] ;  x[2]++ )
    for( x[1] = 0; x[1] < param->size[1] ;  x[1]++ )
      for( i0 = 0; i0 < param->D; i0++ ){
        make_buf( &F_N[link_index( (i0 + 1) % param->D, site_index(x))], buffer);
        if( fwrite( &buffer[0], sizeof(float), N, f) != N ){
          fclose(f);
          panic("Failed to dump gauge matrix");
        };
      };
  fclose(f);
  return 0;
};
/* ******************************************************************* */
/* *****  SU(N)/U(N) Monte Carlo routines  *************************** */
/* ******************************************************************* */
/* Calculates staples for link U(i0,m) so that the local action is
   ~ ReTr[ U staples^+ ]    */
static SU_N staples_one_link( uchar i0, uint m ){
  SU_N ret, tmp;
  uchar i1;
  uint m_left, m_up = index_up( m, i0 );
  bzero( &ret, sizeof(SU_N) );
  for( i1 = 0; i1 < param->D ; i1++){
    if( i1 == i0 ) continue;
    m_left = index_down( m, i1 );
    tmp = SU_N_mult( F_N[link_index(i1,m)], F_N[link_index(i0,index_up(m,i1))] );
    ret = SU_N_sum( ret, SU_N_mult( tmp, SU_N_conj( F_N[link_index(i1,m_up)] ) ) );
    tmp = SU_N_mult( SU_N_conj(F_N[link_index(i1,m_left)]), F_N[link_index(i0, m_left)] );
    ret = SU_N_sum( ret, SU_N_mult( tmp, F_N[link_index(i1,index_up(m_left,i0))]) );
  };
  return ret;
};
#ifdef IMPROVEMENT
/*----------------------------------- Improvement ------------------------------------
                     m_up_up 8
                *-------*-------*  
                |       |       |
                |       |9     7|
           -4   |   m_up|  5    |   4   
        *-------*-------*-------*-------*
        |       |               |       |
      -3|       |      i0      6|       |3
        |  -2   |   -1     1    |   2   |
        *-------*-------*-------*-------*
m_left_left   m_left    m      m_right                                     
		|	|	|
		|	|	|
		*-------*-------*
							*/
static SU_N staples_one_link_rect( uchar i0, uint m ){
  SU_N ret, tmp;
  uchar i1;
  uint m_left,m_left_left, m_right, m_up = index_up( m, i0 ), m_up_up = index_up( m_up, i0 ),m_down=index_down(m,i0);
  bzero( &ret, sizeof(SU_N) );
  for( i1 = 0; i1 < param->D ; i1++){
    if( i1 == i0 ) continue;
    m_left = index_down( m, i1 );
    m_left_left = index_down( m_left, i1 );
    m_right = index_up(m, i1);
    tmp= SU_N_mult( SU_N_mult( SU_N_mult( 
               /* 1*2 */                   F_N[link_index(i1,m)], F_N[link_index(i1,m_right)] ), 
               /*  *3 */                   F_N[link_index(i0,index_up(m_right,i1))] ), 
               /*  *4 */                   SU_N_conj(F_N[link_index(i1,index_up(m_right,i0))]) );
    
    ret = SU_N_sum( ret, SU_N_mult( tmp, SU_N_conj( F_N[link_index(i1,m_up)] ) ) );
    
    tmp = SU_N_mult( SU_N_mult( SU_N_mult( 
               /* -1*-2 */                 SU_N_conj(F_N[link_index(i1,m_left)]), SU_N_conj(F_N[link_index(i1, m_left_left)]) ),
               /*   *-3 */                 F_N[link_index(i0,m_left_left)] ),
               /*   *-4 */                 F_N[link_index(i1,index_up(m_left_left,i0))] );
 
    ret = SU_N_sum( ret, SU_N_mult( tmp, F_N[link_index(i1,index_up(m_left,i0))]) );

    tmp= SU_N_mult( SU_N_mult( SU_N_mult( 
               /* 1*6 */                   F_N[link_index(i1,m)], F_N[link_index(i0,m_right)] ), 
               /*  *7 */                   F_N[link_index(i0,index_up(m_up,i1))] ), 
               /*  *8 */                   SU_N_conj(F_N[link_index(i1,m_up_up)]) );
    
    ret = SU_N_sum( ret, SU_N_mult( tmp, SU_N_conj( F_N[link_index(i0,m_up)] ) ) );
    
    tmp = SU_N_mult( SU_N_mult( SU_N_mult( 
               /* -1*-6 */                 SU_N_conj(F_N[link_index(i1,m_left)]), F_N[link_index(i0, m_left)] ),
               /*   *-7 */                 F_N[link_index(i0,index_up(m_left,i0))] ),
               /*   *-8 */                 F_N[link_index(i1,index_down(m_up_up,i1))] );
 
    ret = SU_N_sum( ret, SU_N_mult( tmp, SU_N_conj( F_N[link_index(i0,m_up)] ) ) );

    tmp= SU_N_mult( SU_N_mult( SU_N_mult( 
                                  SU_N_conj(F_N[link_index(i0,m_down)]), F_N[link_index(i1,m_down)] ), 
                                  F_N[link_index(i0,index_up(m_down,i1))] ), 
                                  F_N[link_index(i0,m_right)] );
    
    ret = SU_N_sum( ret, SU_N_mult( tmp, SU_N_conj( F_N[link_index(i1,m_up)] ) ) );
    
    tmp = SU_N_mult( SU_N_mult( SU_N_mult( 
                                SU_N_conj(F_N[link_index(i0,m_down)]), SU_N_conj(F_N[link_index(i1, index_down(m_left,i0))]) ),
                                F_N[link_index(i0,index_down(m_left,i0))] ),
                                F_N[link_index(i0,m_left)] );
 
    ret = SU_N_sum( ret, SU_N_mult( tmp,  F_N[link_index(i1,index_up(m_left,i0))]  ) );
  };
  return ret;
};
/*-------------------------------Parallelograms----------------------------------------------------*/
/*                 
                               4
                           *-------*
                        5 /.      .|
                    m_up / .     . |3     
                . . . . *. . . .*  |
            -4 .       /   .    .  *  
              .       /i0 .     . /    
             *-------*   . 1    ./ 2    
             |  *-------*-------*
          -3 | /m_left  m      m_right      
             |/
             * m_left_back                                                                                              */
static SU_N staples_one_link_parallelograms( uchar i0, uint m ){
  SU_N ret, tmp;
  uchar i1,i2;
  uint m_left,m_left_back, m_right, m_right_back, m_left_front, m_up = index_up( m, i0 );
  bzero( &ret, sizeof(SU_N) );
  for( i1 = 0; i1 < param->D ; i1++){
    if( i1 == i0 ) continue;
    m_left = index_down( m, i1 );
    m_right = index_up(m, i1);
   for( i2 = 0 ; i2 < param->D ; i2++){
    if( i2 == i0 || i2 == i1 ) continue;
    m_left_back = index_down( m_left, i2 );
    m_right_back = index_down( m_right, i2 );
    m_left_front = index_up( m_left, i2 );
    tmp= SU_N_mult( SU_N_mult( SU_N_mult( 
               /* 1*2 */                   F_N[link_index(i1,m)], F_N[link_index(i2,m_right)] ), 
               /*  *3 */                   F_N[link_index(i0,index_up(m_right,i2))] ), 
               /*  *4 */                   SU_N_conj(F_N[link_index(i1,index_up(m_up,i2))]) );
    
    ret = SU_N_sum( ret, SU_N_mult( tmp, SU_N_conj( F_N[link_index(i2,m_up)] ) ) );
   
    tmp = SU_N_mult( SU_N_mult( SU_N_mult( 
               /* -1*-2 */                 SU_N_conj(F_N[link_index(i1,m_left)]), SU_N_conj(F_N[link_index(i2, m_left_back)]) ),
               /*   *-3 */                 F_N[link_index(i0,m_left_back)] ),
               /*   *-4 */                 F_N[link_index(i1,index_up(m_left_back,i0))] );
 
    ret = SU_N_sum( ret, SU_N_mult( tmp, F_N[link_index(i2,index_down(m_up,i2))]) );

    tmp= SU_N_mult( SU_N_mult( SU_N_mult( 
               		                   F_N[link_index(i1,m)], SU_N_conj(F_N[link_index(i2,m_right_back)]) ), 
               		                  F_N[link_index(i0,m_right_back)] ), 
              		                  SU_N_conj(F_N[link_index(i1,index_down(m_up,i2))]) );
    
    ret = SU_N_sum( ret, SU_N_mult( tmp,  F_N[link_index(i2,index_down(m_up,i2))]  ) );
   
    tmp = SU_N_mult( SU_N_mult( SU_N_mult( 
               			                SU_N_conj(F_N[link_index(i1,m_left)]), F_N[link_index(i2, m_left)] ),
              			                F_N[link_index(i0,m_left_front)] ),
               			                F_N[link_index(i1,index_up(m_left_front,i0))] );
 
    ret = SU_N_sum( ret, SU_N_mult( tmp, SU_N_conj( F_N[link_index(i2,m_up)])) );
   };
  };
  return ret;
};
// This computes twisted rectangular staples. Not used in action
static SU_N staples_one_link_twist( uchar i0, uint m ){
  SU_N ret, tmp;
  uchar i1;
  uint m_left,m_left_left, m_right, m_up = index_up( m, i0 ), m_up_up = index_up( m_up, i0 ),m_down=index_down(m,i0);
  bzero( &ret, sizeof(SU_N) );
  for( i1 = 0; i1 < param->D ; i1++){
    if( i1 == i0 ) continue;
    m_left = index_down( m, i1 );
    m_left_left = index_down( m_left, i1 );
    m_right = index_up(m, i1);
    tmp= SU_N_mult( SU_N_mult( SU_N_mult( 
                                  F_N[link_index(i1,m)], F_N[link_index(i1,index_up(m_right,i0))] ), 
                                SU_N_conj(F_N[link_index(i0,index_up(m_right,i1))]) ), 
                                 SU_N_conj(F_N[link_index(i1,m_up)]) );
    
    ret = SU_N_sum( ret, SU_N_mult( tmp, SU_N_conj( F_N[link_index(i1,m_up)] ) ) );
   
    tmp = SU_N_mult( SU_N_mult( SU_N_mult( 
                               SU_N_conj(F_N[link_index(i1,m_left)]), SU_N_conj(F_N[link_index(i1, index_up(m_left_left,i0))]) ),
                               SU_N_conj(F_N[link_index(i0,m_left_left)]) ),
                               F_N[link_index(i1,m_left_left)] );
 
    ret = SU_N_sum( ret, SU_N_mult( tmp, F_N[link_index(i1,index_up(m_left,i0))]) );

    tmp= SU_N_mult( SU_N_mult( SU_N_mult( 
                                 F_N[link_index(i1,m)], F_N[link_index(i0,m_right)] ), 
                                 F_N[link_index(i0,m_up)] ), 
                                 F_N[link_index(i1,m_up_up)] );
    
    ret = SU_N_sum( ret, SU_N_mult( tmp, SU_N_conj( F_N[link_index(i0,index_up(m_up,i1))] ) ) );
    
    tmp = SU_N_mult( SU_N_mult( SU_N_mult( 
                                SU_N_conj(F_N[link_index(i1,m_left)]), F_N[link_index(i0, m_left)] ),
                                F_N[link_index(i0,m_up)] ),
                                SU_N_conj(F_N[link_index(i1,index_down(m_up_up,i1))]) );
 
   ret = SU_N_sum( ret, SU_N_mult( tmp, SU_N_conj( F_N[link_index(i0,index_up(m_left,i0))] ) ) );
   tmp= SU_N_mult( SU_N_mult( SU_N_mult( 
                                  SU_N_conj(F_N[link_index(i0,index_up(m_down,i1))]), SU_N_conj(F_N[link_index(i1,m_down)]) ), 
                                  F_N[link_index(i0,m_down)] ), 
                                  F_N[link_index(i0,m_right)] );
    
    ret = SU_N_sum( ret, SU_N_mult( tmp, SU_N_conj( F_N[link_index(i1,m_up)] ) ) );
    
    tmp = SU_N_mult( SU_N_mult( SU_N_mult( 
                                SU_N_conj(F_N[link_index(i0,index_down(m_left,i0))]), F_N[link_index(i1, index_down(m_left,i0))] ),
                                F_N[link_index(i0,m_down)] ),
                                F_N[link_index(i0,m_left)] );
 
    ret = SU_N_sum( ret, SU_N_mult( tmp,  F_N[link_index(i1,index_up(m_left,i0))]  ) );
  };
  return ret;
};
#endif
/* ******************************************************************* */
/* *****  SU(2) (p,q) subgroup heatbath/cooling  ********************* */
/* ******************************************************************* */
static void SU2_pq_heatbath( SU_N *U, SU_N *staples, SU_N *rectangles, SU_N *parallelograms, uchar p, uchar q, uchar cool){
  SU2 h, tmp;
  double alpha, a0, help, phi, cos_theta, sin_theta;

  uchar i;
  h.alpha = h.beta = 0.0;


  for( i = 0; i < SU_N_RANK; i++ ){
    h.alpha += conj(U->U[i][p]) * staples->U[i][p] + U->U[i][q] * conj(staples->U[i][q]);
    h.beta  += conj(U->U[i][p]) * staples->U[i][q] - U->U[i][q] * conj(staples->U[i][p]);
  };
  alpha = SU2_normalize( &h ) * param->beta/((double) SU_N_RANK);
  if( alpha < FLT_EPSILON ) panic("Degenerate SU2 staples");
  if( cool ){
    h.alpha *= MC_SU_N_COOLING_DELTA;
    h.beta  *= MC_SU_N_COOLING_DELTA;
    h.alpha += 1.0;
    SU2_normalize( &h );
  }else{
    a0 = (alpha < ALPHA_CROSS ) ? Creutz(alpha) : Kennedy(alpha);
    help = sqrt( 1.0 - a0 * a0 );
    phi = M_PI * (2.0 * RND() - 1.0);
    cos_theta = 2.0 * RND() - 1.0;
    sin_theta = sqrt( 1.0 - cos_theta * cos_theta );
    tmp.alpha = a0 + I * help * cos_theta;
    tmp.beta = help * sin_theta * cos(phi) + I * help * sin_theta * sin(phi);
    h = SU2_mult( tmp, h );
  };
  for( i = 0; i < SU_N_RANK; i++ ){
    tmp.alpha = U->U[i][p]; /* store temporaly */
    U->U[i][p] *= h.alpha;
    U->U[i][p] -= conj(h.beta) * U->U[i][q];
    U->U[i][q] *= conj(h.alpha);
    U->U[i][q] += h.beta * tmp.alpha;
  };
};
/* ******************************************************************* */
/* ****  U(1) (p) subgroup heatbath/cooling  ************************* */
/* ******************************************************************* */
static inline double complex U1_heatbath( double a ){
  double help, cs;
  do{
    help = M_PI * ( 2.0 * RND() - 1.0 );
    cs = cos(help);
  }while( cs < 1.0 + log( RND() )/a );
  return cs + I * sin(help);
};
static void U1_p_heatbath( SU_N *U, SU_N *staples, uchar p, uchar cool ){
  double complex z, zeta = 0.0;
  double alpha;
  uchar i;
  for( i = 0; i < SU_N_RANK; i++ ) zeta += conj( U->U[i][p] ) * staples->U[i][p];
  if( (alpha = cabs(zeta)) < FLT_EPSILON ) panic("Degenerate U1 staples");
  zeta /= alpha;
  if( cool ){
    z = 1.0 + MC_SU_N_COOLING_DELTA * zeta;
    alpha = cabs(z);
    z /= alpha;
  }else
    z = zeta * U1_heatbath( param->beta * alpha/((double) SU_N_RANK ) );
  for( i = 0; i < SU_N_RANK; i++ ) U->U[i][p] *= z;
};
/* ******************************************************************* */
/* ******  U(N) overrelaxation is almost trivial ********************* */
/* ******  For SU(N), see, e.g.,  hep-lat/0308033, hep-lat/0503041 *** */
/* ******************************************************************* */
static void overrelax_U_N( SU_N *U, SU_N staples ){
  if( SU_N_projection_U_N( &staples ) ) panic("Degenerate staples");
  *U = SU_N_mult( SU_N_conj(*U), staples );
  *U = SU_N_mult( staples, *U );
};
static void overrelax_SU_N( SU_N *U, SU_N staples, SU_N rectangles, SU_N parallelograms){
  uchar i, k;
  double S_re, S_im, S_new;
  double complex z = 0.0;
  for( i = 0; i < SU_N_RANK; i++ )
    for( k = 0; k < SU_N_RANK; k++ )
     { 
      z += U->U[i][k] * conj(staples.U[i][k]);
     };
  S_re = creal(z);
  S_im = cimag(z);
  if( SU_N_projection_U_N( &staples ) ) panic("Degenerate staples");
  z = SU_N_determinant( staples );
  S_new = 2.0 * atan2( cimag(z), creal(z) )/((double) SU_N_RANK);
  z = cos(S_new) - I * sin(S_new);
  S_new = S_re * creal(z) +  S_im * cimag(z);
  if( S_new - S_re > log(RND()) ){
    *U = SU_N_mult( SU_N_conj(*U), staples );
    *U = SU_N_mult( staples, *U );
    for( i = 0; i < SU_N_RANK; i++ )
      for( k = 0; k < SU_N_RANK; k++ )
    U->U[i][k] *= z;
  };
};
/* ******************************************************************* */
/* *****  Driver for both SU(N)/U(N) Monte-Carlo  ******************** */
/* ******************************************************************* */
static double MC_real( uchar su, uchar cool){
  uint m;
  uchar i0, p, q;
  SU_N *U, staples,rectangles,parallelograms;

  double ret = 0.0;
  if( !param || !F_N ) panic("Geometry/gauge data undefined");
  for( m = 0; m < param->ipw[param->D]; m++)
   for( i0 = 0; i0 < param->D ; i0++ )
   {
    U = &F_N[ link_index(i0, m) ];
    staples = staples_one_link( i0, m );
#ifdef IMPROVEMENT
      rectangles = staples_one_link_rect(i0,m);
      parallelograms = staples_one_link_parallelograms(i0,m);
  uchar i, k;
  for( i = 0; i < SU_N_RANK; i++ )
    for( k = 0; k < SU_N_RANK; k++ )
     { 
      staples.U[i][k] = param->beta1*staples.U[i][k];
      rectangles.U[i][k] = rectangles.U[i][k]*param->beta2;
      parallelograms.U[i][k] = parallelograms.U[i][k]*param->beta3;
     };
    staples = SU_N_sum(staples,SU_N_sum(rectangles,parallelograms));
#endif
    for( p = 0; p < SU_N_RANK; p++ )
    {
     for( q = p+1; q < SU_N_RANK; q++ ){
      SU2_pq_heatbath( U, &staples, &rectangles, &parallelograms, p, q, cool);};
     if( !su ) U1_p_heatbath( U, &staples, p, cool );
    };
    if( !cool )
    {
     if( su ){
      overrelax_SU_N( U, staples, rectangles, parallelograms);}
     else
      overrelax_U_N(  U, staples );
    };
    ret += SU_N_norm( SU_N_mult( *U, SU_N_conj(staples) ) );

   };
  return ret/( 2 * param->D * (param->D-1) * param->ipw[param->D]);
};
/* ------------------------------------------------------------------- */
double MC_U_N(){       return MC_real( 0, 0); };
double MC_SU_N(){      return MC_real( 1, 0);};
double MC_U_N_cool(){  return MC_real( 0, 1); };
double MC_SU_N_cool(){ return MC_real( 1, 1); };
/* ******************************************************************* */
/* ******   Various routines   *************************************** */
/* ******************************************************************* */
double SU_N_mean_plaq(){
  uint m,m1,l;
  uchar i0, i1;
  double mean = 0.0;
  SU_N tmp1, tmp2;
  if( !param || !F_N ) panic("Geometry/fields are not inited");
  for( m = 0; m < param->ipw[ param->D ]; m++)
    for( i0 = 0 ; i0 < param->D-1 ; i0++ ){
      l = link_index(i0,m);
      m1 = index_up( m, i0 );
      for( i1 = i0+1 ; i1 < param->D  ; i1++ ){
        tmp1 = SU_N_mult(F_N[l], F_N[link_index(i1,m1)] );
        tmp2 = SU_N_mult(F_N[link_index(i1,m)], F_N[link_index(i0,index_up(m,i1))] );
        mean += SU_N_norm( SU_N_mult( tmp1, SU_N_conj(tmp2) ) );
      };
    };
  return 2.0 * mean/( param->D * (param->D-1) * param->ipw[param->D]);
};
/* ******************************************************************* */
double U_N_mean_plaq(){
  return SU_N_mean_plaq();
};
/* ******************************************************************* */
/* ****  1-st Chern for U(N) fields  ********************************* */
/* ******************************************************************* */
static double U_N_plaq_angle( uint m, uchar i0, uchar i1 ){
  double complex z;
  SU_N U1 = F_N[ link_index(i0, m) ];
  SU_N U2 = F_N[ link_index(i1, m) ];
  U1 = SU_N_mult( U1, F_N[ link_index(i1, index_up(m, i0)) ] );
  U2 = SU_N_mult( U2, F_N[ link_index(i0, index_up(m, i1)) ] );
  z = SU_N_determinant( SU_N_mult( U1, SU_N_conj(U2) ) );
  return atan2( cimag(z), creal(z) )/(2.0 * M_PI);
};
/* --------------------------------------------------------------------- */
static int U_N_mono_cube( uint m, uchar i0 ){
  double flux = 0.0;
  uchar i1, i2, i3;
  int ret;
  for( i1 = 0; i1 < param->D; i1++ ){
    if( i1 == i0 ) continue;
    D4dual( i0, i1, &i2, &i3 );
    flux += U_N_plaq_angle(          m,     i2, i3 );
    flux += U_N_plaq_angle( index_up(m,i1), i3, i2 );
  };
  ret = (int) rint( flux );
#ifdef DEBUG
  if( fabs( ret - flux ) > FLT_EPSILON ) panic("Not integer");
#endif
  return ret;
};
/* --------------------------------------------------------------------- */
static char *U_N_mono_current = NULL;

double U_N_mono( char **current ){
  uint m;
  uchar i;
  double ret = 0.0;
  char val;
  if( !param || !F_N ) panic("Geometry/fields are not inited");
  if( param->D != 4 ) panic("For D=4 ONLY");
  if( current && !U_N_mono_current
      && !( U_N_mono_current = (char *) calloc( param->D * param->ipw[param->D], sizeof(char))))
    panic("Memory");
  for( m = 0; m < param->ipw[ param->D ]; m++)
    for( i = 0; i < param->D; i++ ){
      val =  U_N_mono_cube( m, i );
      ret += abs( val );
      if( current ) U_N_mono_current[link_index(i, index_down(m,i))] = val;
    };
#ifdef DEBUG
  if( current ){
    int divergence;
    for( m = 0; m < param->ipw[ param->D ]; m++){
      divergence = 0;
      for( i = 0; i < param->D; i++ ){
    divergence += U_N_mono_current[ link_index(i, m) ];
    divergence -= U_N_mono_current[ link_index(i, index_down(m,i)) ];
      };
      if( divergence ) panic("Conservation violation");
    };
  };
#endif
  if( current ) *current = U_N_mono_current;
  return ret/((double) ( param->D * param->ipw[param->D] ));
};
/* ******************************************************************* */
static void SU_N_random_gauge_real( uchar su ){
  SU_N Omega;
  uint m, m1;
  uchar i0;
  if( !param || !F_N ) panic("Geometry/fields are not inited");
  for( m = 0; m < param->ipw[ param->D ]; m++){
    if( su ) SU_N_random( &Omega );
    else     SU_N_random_U_N( &Omega );
    for( i0 = 0; i0 < param->D ; i0++){
      m1 = link_index(i0,m);
      F_N[m1] = SU_N_mult( SU_N_conj(Omega), F_N[m1] );
      m1 = link_index(i0,index_down(m,i0));
      F_N[m1] = SU_N_mult( F_N[m1], Omega );
    };
  };
};
/* ******************************************************************* */
void SU_N_random_gauge(){
  SU_N_random_gauge_real(1);
};
void U_N_random_gauge(){
  SU_N_random_gauge_real(0);
};
/* ******************************************************************* */
static double SU_N_Landau_gauge_real( uchar su ){
  SU_N g_left, g_right;
  double error, help, value;
  uint m, l;
  uchar i0;
#ifdef VERBOSE
  uint count = 0;
#endif
  if( !param || !F_N )  panic("Parameters/fields are not inited");
#ifdef VERBOSE
  printf("Landau gauge fixing "); fflush(stdout);
#endif
  do{
    value = error = 0.0;
    for( m = 0; m < param->ipw[ param->D ]; m++){
      bzero( &g_left, sizeof(SU_N) );
      for( i0 = 0; i0 < param->D ; i0++){
        l = link_index( i0, m );
    g_left = SU_N_sum( g_left, SU_N_conj( F_N[l] ));
        l = link_index( i0, index_down(m,i0) );
    g_left = SU_N_sum( g_left, F_N[l]);
      };
      if( su ) SU_N_normalize( &g_left );
      else     SU_N_projection_U_N( &g_left );
      g_right = SU_N_conj( g_left );
      if( (help=(1.0 - SU_N_norm(g_left))) > error ) error = help;
      for( i0 = 0; i0 < param->D ; i0++){
        l = link_index( i0, m );
        F_N[l] = SU_N_mult( g_left, F_N[l] );
    value += SU_N_norm( F_N[l] );
        l = link_index( i0, index_down(m,i0) );
        F_N[l] = SU_N_mult( F_N[l], g_right );
    value += SU_N_norm( F_N[l] );
      };
    };
    value /= (double)(2 * param->D * param->ipw[ param->D ]);
    value = 1.0 - value;
#ifdef VERBOSE
    if( !((count++) % 20) ){ printf("%.4f ", value ); fflush(stdout); };
#endif
  }while( error > LANDAU_GAUGE_PRECISION );
#ifdef VERBOSE
  printf(" Done, iteration count: %d\n", count); fflush(stdout);
#endif
  return value;
};
double SU_N_Landau_gauge(){
  return SU_N_Landau_gauge_real(1);
};
double  U_N_Landau_gauge(){
  return SU_N_Landau_gauge_real(0);
};
