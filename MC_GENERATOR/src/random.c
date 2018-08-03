#include <random.h>

typedef struct {
    unsigned long int max;
    unsigned long int min;
    size_t size;
    void (*set) (void *state, unsigned long int seed);
    double (*get_double) (void *state);
} s_rnd_type;

typedef struct{
  const s_rnd_type *type;
  void *state;
} s_rnd;

#define N 624   /* Period parameters */
#define M 397

static const unsigned long UPPER_MASK = 0x80000000UL;   
static const unsigned long LOWER_MASK = 0x7fffffffUL;   

typedef struct {
  unsigned long mt[N];
  int mti;
} mt_state_t;

static inline unsigned long int mt_get (void *vstate);
static double mt_get_double (void *vstate);
static void mt_set (void *state, unsigned long int s);

static const s_rnd_type mt_type = {
  0xffffffffUL,                  /* RAND_MAX  */
  0,                             /* RAND_MIN  */
  sizeof (mt_state_t),
  &mt_set,
  &mt_get_double
};


static s_rnd *mt_rnd = NULL;

/* *************************************************************** */
static inline unsigned long mt_get (void *vstate){
  mt_state_t *state = (mt_state_t *) vstate;
  unsigned long k ;
  unsigned long int *const mt = state->mt;

#define MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)

  if ( state->mti >= N ){   /* generate N words at one time */
    int kk;
    for (kk = 0; kk < N - M; kk++) {
      unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + M] ^ (y >> 1) ^ MAGIC(y);
    }
    for (; kk < N - 1; kk++) {
      unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ MAGIC(y);
    }
    {
      unsigned long y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ MAGIC(y);
    }
    state->mti = 0;
  }
  /* Tempering */
  k = mt[state->mti];
  k ^= (k >> 11);
  k ^= (k << 7) & 0x9d2c5680UL;
  k ^= (k << 15) & 0xefc60000UL;
  k ^= (k >> 18);
  state->mti++;
  return k;
}
/* *************************************************************** */
static double mt_get_double (void * vstate){
  return mt_get (vstate) / 4294967296.0 ;
}
/* *************************************************************** */
static void mt_set (void *vstate, unsigned long int s){
  mt_state_t *state = (mt_state_t *) vstate;
  int i;
  if (s == 0) s = 4357;   /* the default seed is 4357 */
  state->mt[0]= s & 0xffffffffUL;
  for (i = 1; i < N; i++) {
    state->mt[i] = (1812433253UL * (state->mt[i-1] ^ (state->mt[i-1] >> 30)) + i);
    state->mt[i] &= 0xffffffffUL;
  }
  state->mti = i;
}
/* *************************************************************** */
static void rnd_free (){
  if( mt_rnd ){
    free( mt_rnd->state );
    free( mt_rnd );
  };
  mt_rnd = NULL;
};
/* *************************************************************** */
void rnd_seed_set(){
  struct timeval tm;
  if( gettimeofday( &tm, NULL) ) tm.tv_usec = (long) time(NULL);
  (mt_rnd->type->set)( mt_rnd->state, (unsigned long int) tm.tv_usec );
};
/* *************************************************************** */
void rnd_init(){
  struct timeval tm;
  if( gettimeofday( &tm, NULL) ) tm.tv_usec = (long) time(NULL);
#if 1
  tm.tv_usec = 0;
  fprintf(stderr, "------------  WARNING - random seed is NOT set ---------------\n");
#endif
  rnd_free();
  if( !(mt_rnd = (s_rnd *) calloc( 1, sizeof(s_rnd)))
      || !(mt_rnd->state = calloc( 1, mt_type.size)) ){
    fprintf(stderr, "Memory allocation failed\n");
    exit(1);
  };
  mt_rnd->type = &mt_type;
  (mt_rnd->type->set)( mt_rnd->state, (unsigned long int) tm.tv_usec );
}
/* *************************************************************** */
double rnd_get(){
  double x;
  do{
    x = (mt_rnd->type->get_double)( mt_rnd->state );
  }while( x == 0 );
  return x;
};
