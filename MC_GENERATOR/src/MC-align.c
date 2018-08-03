#include <MC-align.h>
/* ------------------------------------------------------------------------------------ */
typedef struct {
  SU2 F;
  double complex *z1;
  double complex *z2;
  double *plaq_phase;
  double *triangle_phase;
  double *square_phase;

} MC_align_Bond_Data;

static MC_align_Bond_Data *MC_align_bond_data = NULL;
static double MC_align_delta = 1.0;
/* ************************************************************************************* */
/* **********   Extended Monte Carlo  (ALIGN) ****************************************** */
/* ************************************************************************************* */
typedef struct {
  uint accept, try;
  double mean_plaq, mean_phase;
} MC_align_Monitor;

static MC_align_Monitor MC_align_monitor;
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static inline SU2 MC_align_get_staples( Bond *b ){
  uchar k;
  SU2 ret, tmp;
  Plaq *p;
  ret.alpha = ret.beta = 0;
  for( k = 0; k < 2*(param->D-1); k++ ){
    p = b->plaq[k];
    switch( b->self[k] ){
    case 0:
      tmp = SU2_mult(*(p->bond[1]->F), SU2_conj( SU2_mult( *(p->bond[3]->F),*(p->bond[2]->F) ) ));
      break;
    case 1:
      tmp = SU2_mult(SU2_conj( SU2_mult( *(p->bond[3]->F),*(p->bond[2]->F) )), *(p->bond[0]->F));
      break;
    case 2:
      tmp = SU2_mult(SU2_conj( SU2_mult( *(p->bond[0]->F),*(p->bond[1]->F) )), *(p->bond[3]->F));
      break;
    case 3:
      tmp = SU2_mult(*(p->bond[2]->F), SU2_conj( SU2_mult( *(p->bond[0]->F),*(p->bond[1]->F) ) ));
      break;
    default: panic("What's the hell?!");
    };
    ret.alpha += tmp.alpha;
    ret.beta += tmp.beta;
  };
  return ret;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static inline double MC_align_update_bond_variables( Bond *b, double gamma ){
  uchar i, k;
  Plaq *p;
  double action = 0.0;
  for( k = 0; k < 2 * (param->D-1) ; k++ ){
    p = b->plaq[k];
    for( i=0; i < 4; i++ ){
      MC_align_bond_data->z1[4*k+i] = p->Sp[i]->z1;
      MC_align_bond_data->z2[4*k+i] = p->Sp[i]->z2;
    };
    MC_align_bond_data->plaq_phase[k] = p->phase;
    update_plaq( p );
  };
  for( k = 0; k < N_triangle ; k++ ){
    MC_align_bond_data->triangle_phase[k] = b->triangle[k]->phase;
    update_triangle( b->triangle[k] );
    action += measure_triangle( b->triangle[k] );
  };
  for( k = 0; k < N_square ; k++ ){
    MC_align_bond_data->square_phase[k] = b->square[k]->phase;
    update_square( b->square[k] );
    action += measure_square( b->square[k] );
  };
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  return gamma * action;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static void MC_align_restore_bond_variables( Bond *b ){
  uchar i, k;
  Plaq *p;
  for( k = 0; k < 2 * (param->D-1) ; k++ ){
    p = b->plaq[k];
    for( i=0; i < 4; i++ ){
      p->Sp[i]->z1 = MC_align_bond_data->z1[4*k+i];
      p->Sp[i]->z2 = MC_align_bond_data->z2[4*k+i];
    };
    p->phase = MC_align_bond_data->plaq_phase[k];
  };
  for( k = 0; k < N_triangle ; k++ )
    b->triangle[k]->phase = MC_align_bond_data->triangle_phase[k];
  for( k = 0; k < N_square ; k++ )
    b->square[k]->phase = MC_align_bond_data->square_phase[k];
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* These will enter into the action multiplied by gamma */
static inline double MC_align_bond_gamma( Bond *b ){
  double ret = 0.0;
  uchar k;
  for( k = 0; k < N_triangle ; k++ ) ret += measure_triangle( b->triangle[k] );
  for( k = 0; k < N_square ; k++ )   ret += measure_square( b->square[k] );
  return ret;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static void MC_align_bond_one_link( uchar i0, uint m, double gamma ){
  uchar i, relax;
  SU2 Staples, g, g_new;
  double S_old_gauge, S_old_phase, S_new_gauge, S_new_phase, alpha;
  Bond *b;
#ifdef VERBOSE
  double help;
#endif
  b = &bond[link_index(i0,m)];
  Staples = MC_align_get_staples( b );
  alpha = SU2_normalize( &Staples );
  g = SU2_mult( *(b->F), Staples);
#ifdef VERBOSE
  help = MC_align_bond_gamma(b);
  MC_align_monitor.mean_phase += help;
  S_old_phase = gamma * help;
  //---------------------------------
  help = alpha * creal( g.alpha );
  MC_align_monitor.mean_plaq += help;
  S_old_gauge = param->beta * help;
#else
  S_old_phase  = gamma * MC_align_bond_gamma(b);
  S_old_gauge  = param->beta * alpha * creal( g.alpha );
#endif
  for( i = 0; i < MC_ALIGN_N_TRY ; i++ ){
    relax = i % 2;
    MC_align_bond_data->F = *(b->F);
    if( relax ){ /* Overrelaxation */
      *(b->F) = SU2_conj( SU2_mult( Staples, g ) );
      g_new = SU2_conj( g );
      S_new_gauge = S_old_gauge;
    }else{   /* Random MC drift */
      *(b->F) = SU2_mult( *(b->F), SU2_random_matrix(MC_align_delta) );
      g_new = SU2_mult( *(b->F), Staples);
      S_new_gauge   = param->beta * alpha * creal( g_new.alpha );
      MC_align_monitor.try++ ;
    };
    S_new_phase = MC_align_update_bond_variables( b, gamma );
    if( (S_new_phase + S_new_gauge - S_old_phase - S_old_gauge) > log(RND()) ){
      if( !relax ) MC_align_monitor.accept++ ;
      S_old_phase = S_new_phase;
      S_old_gauge = S_new_gauge;
      g.alpha = g_new.alpha;
      g.beta = g_new.beta;
    }else{
      *(b->F) = MC_align_bond_data->F;
      MC_align_restore_bond_variables( b );
    };
  }; /* i (N_TRY) */
  SU2_normalize( b->F );
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static void MC_align_bond( double gamma ){
  uint m;
  uchar i0;
  double help;
  bzero( &MC_align_monitor, sizeof(MC_align_Monitor) );
  for( m = 0; m < param->ipw[param->D]; m++ )
    for( i0 = 0; i0 < param->D; i0++ )
      MC_align_bond_one_link( i0, m, gamma );
#ifdef VERBOSE
  help = param->D * param->ipw[ param->D ];
  MC_align_monitor.mean_plaq /= 2*(param->D-1)*help;
  if( param->D == 3 ){
    MC_align_monitor.mean_phase /= (3 * 24 + 4 * 26) * help;
  }else{
    MC_align_monitor.mean_phase /= (3 * 232 + 4 * 78) * help;
  };
#endif
  help = ((double)MC_align_monitor.accept) / ((double)MC_align_monitor.try);
  if( help > 0.6 ) MC_align_delta *= 1.25;
  if( help < 0.4 ) MC_align_delta *= 0.75;
  if( MC_align_delta > 2 ) MC_align_delta = 2;
  if( MC_align_delta < 0.01 ) MC_align_delta = 0.01;
#ifdef VERBOSE
  printf("Accept=%d%%, delta=%.3f, beta=%.3f, gamma=%.3f, "
	 "1/2<TrU_p>=%.5f, rho=%.5f\n",
	 ((int)(help*100)), MC_align_delta, param->beta, gamma,
	 MC_align_monitor.mean_plaq, MC_align_monitor.mean_phase );
  fflush( stdout );
#endif
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
int MC_align( uint N_sweeps, double gamma ){
  uint sweep;
  if( param->D != 3 && param->D != 4 ){
    error("Hey, man, this works only in D=3 and D=4!"); return -1;
  };
  if( !spinor || !bond || !plaq || !triangle || !square ){
    error("Hey man, are you in hurry?! (structures undefined)"); return -1;
  };
  if( !MC_align_bond_data ){
    MEMORY_CHECK(MC_align_bond_data=(MC_align_Bond_Data *)calloc(1,sizeof(MC_align_Bond_Data)),-1);
    MEMORY_CHECK(MC_align_bond_data->z1 = (double complex *)calloc(8*(param->D-1), sizeof(double complex)),-1);
    MEMORY_CHECK(MC_align_bond_data->z2 = (double complex *)calloc(8*(param->D-1), sizeof(double complex)),-1);
    MEMORY_CHECK(MC_align_bond_data->plaq_phase = (double *)calloc(2*(param->D-1), sizeof(double)),-1);
    if( N_triangle ) MEMORY_CHECK(MC_align_bond_data->triangle_phase = (double *)calloc( N_triangle, sizeof(double)),-1);
    if( N_square   ) MEMORY_CHECK(MC_align_bond_data->square_phase   = (double *)calloc( N_square,   sizeof(double)),-1);
  };
  for( sweep = 0; sweep < N_sweeps; sweep++ ){
#ifdef VERBOSE
    printf("%3d of %3d:", sweep, N_sweeps );
    fflush( stdout );
#endif
    MC_align_bond( gamma );
  };
  return 0;
};
