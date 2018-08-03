#include <MC-killing.h>
/*
   DEBUG                - include debugging code (slow)

   By default gamma_s kills all diamonds and sites, but diamonds and sites physically
   make two different groups. Each one is responsible for zeros of:
   group 1) (E1 E2 E3), (E1 B2 B3), (B1 E2 B3), (B1 B2 E3) -- diamonds in the second group (16 - 23)
            (B1 B2 B3), (B1 E2 E3), (E1 B2 E3), (E1 E2 B3) -- sites
	    Note that Lorenz components here differ.
   group 2) (B1 E2 B2), (B1 E3 B3), (B2 E1 B1), (B2 E3 B3), (B3 E1 B1), (B3 E2 B2),
            (E1 E2 B2), (E1 E3 B3), (E2 E1 B1), (E2 E3 B3), (E3 E1 B1), (E3 E2 B2)
	    diamonds in the first group (0 - 15)
   Thus losely speaking group 2 have something to do with self-duality - suppression of mag.
   charges in this group causes the fields to be far away from self-duality. Therefore:
   MC_DIAMOND - gamma_s suppresses only 'duality' diamonds (group 2)
   MC_SITE    - gamma_s suppresses only sites and diamonds from group 1
*/
#if (defined(MC_DIAMOND))&&(defined(MC_SITE))
#error "You CANNOT define both MC_DIAMOND and MC_SITE"
#endif
/* ------------------------------------------------------------------------------------ */
typedef struct {
  SU2 F;
  double complex *z1;
  double complex *z2;
  double *plaq_phase;
  double *triangle_phase;
  double *square_phase;
  char   *q_site;
  char   *q_cube;
  char   *q_diamond;
  /* Topology is not implemented, overwise you'll need these
     char   *q_hypersite;
     char   *q_hypercube;
  */
} MC_killing_Bond_Data;
static MC_killing_Bond_Data *MC_killing_bond_data = NULL;
static double MC_killing_delta = 1.0;
/* ************************************************************************************* */
/* **********   Extended Monte Carlo  ************************************************** */
/* ************************************************************************************* */
typedef struct {
  uint accept, try;
  double mean_plaq, mean_s, mean_c;
} MC_killing_Monitor;

static MC_killing_Monitor MC_killing_monitor;
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static inline SU2 MC_killing_get_staples( Bond *b ){
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
static inline double MC_killing_update_bond_variables( Bond *b, double gamma_s, double gamma_c ){
  uchar i, k;
  Plaq *p;
  double sq = 0.0, cq = 0.0;
  for( k = 0; k < 2 * (param->D-1) ; k++ ){
    p = b->plaq[k];
    for( i=0; i < 4; i++ ){
#if defined(SPINOR_SECTION_PLUS)
      MC_killing_bond_data->z2[4*k+i] = p->Sp[i]->z2;
#elif defined(SPINOR_SECTION_MINUS)
      MC_killing_bond_data->z1[4*k+i] = p->Sp[i]->z1;
#else
      MC_killing_bond_data->z1[4*k+i] = p->Sp[i]->z1;
      MC_killing_bond_data->z2[4*k+i] = p->Sp[i]->z2;
#endif
    };
    MC_killing_bond_data->plaq_phase[k] = p->phase;
    update_plaq( p );
  };
  for( k = 0; k < N_triangle ; k++ ){
    MC_killing_bond_data->triangle_phase[k] = b->triangle[k]->phase;
    update_triangle( b->triangle[k] );
  };
  for( k = 0; k < N_square ; k++ ){
    MC_killing_bond_data->square_phase[k] = b->square[k]->phase;
    update_square( b->square[k] );
  };
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for( k = 0; k < N_site ; k++ ){
    MC_killing_bond_data->q_site[k] = b->site[k]->q;
    if( param->D == 3 ){
      sq += fabs( site_charge( b->site[k] ) );
    }else{
#if !defined(MC_DIAMOND)
    sq += fabs( site_charge( b->site[k] ) );
#else
    site_charge( b->site[k] );
#endif
    };
  };
  for( k = 0; k < N_cube ; k++ ){
    MC_killing_bond_data->q_cube[k] = b->cube[k]->q;
    cq += fabs( cube_charge( b->cube[k] ) );
  };
  for( k = 0; k < N_diamond ; k++ ){
    MC_killing_bond_data->q_diamond[k] = b->diamond[k]->q;
    i = fabs( diamond_charge( b->diamond[k] ) );
#if defined(MC_DIAMOND)
    if( !(b->diamond[k]->group & 0x01) ) sq += i;
#elif defined(MC_SITE)
    if( b->diamond[k]->group & 0x01 ) sq += i;
#else
    sq += i;
#endif
  };
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef DEBUG
  for( k = 0; k < N_hypersite ; k++ ) hypersite_charge( b->hypersite[k] );
  for( k = 0; k < N_hypercube ; k++ ) hypercube_charge( b->hypercube[k] );
#endif
  return gamma_s * sq + gamma_c * cq;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static void MC_killing_restore_bond_variables( Bond *b ){
  uchar i, k;
  Plaq *p;
  for( k = 0; k < 2 * (param->D-1) ; k++ ){
    p = b->plaq[k];
    for( i=0; i < 4; i++ ){
#if defined(SPINOR_SECTION_PLUS)
      p->Sp[i]->z2 = MC_killing_bond_data->z2[4*k+i];
#elif defined(SPINOR_SECTION_MINUS)
      p->Sp[i]->z1 = MC_killing_bond_data->z1[4*k+i];
#else
      p->Sp[i]->z1 = MC_killing_bond_data->z1[4*k+i];
      p->Sp[i]->z2 = MC_killing_bond_data->z2[4*k+i];
#endif
    };
    p->phase = MC_killing_bond_data->plaq_phase[k];
  };
  for( k = 0; k < N_triangle ; k++ )
    b->triangle[k]->phase = MC_killing_bond_data->triangle_phase[k];
  for( k = 0; k < N_square ; k++ )
    b->square[k]->phase = MC_killing_bond_data->square_phase[k];
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for( k = 0; k < N_site ; k++ )
    b->site[k]->q = MC_killing_bond_data->q_site[k];
  for( k = 0; k < N_cube ; k++ )
    b->cube[k]->q = MC_killing_bond_data->q_cube[k];
  for( k = 0; k < N_diamond ; k++ )
    b->diamond[k]->q = MC_killing_bond_data->q_diamond[k];
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Normalization factors of gamma_s action for various defines */               
#if defined(MC_DIAMOND)                                                         
#define TOTAL_GAMMA_S   64                                                      
#elif defined(MC_SITE)                                                          
#define TOTAL_GAMMA_S  (38 + 30)                                                
#else                                                                           
#define TOTAL_GAMMA_S  (102 + 30)                                               
#endif                              
/* These will enter into the action multiplied by gamma_s */
static inline uint MC_killing_bond_gamma_s( Bond *b ){
  uint ret = 0;
  uchar k;
  if( param->D == 3 ){
    for( k = 0; k < N_site    ; k++ ) ret += abs( b->site[k]->q );
  }else{
#if defined(MC_DIAMOND)
    for( k = 0; k < N_diamond ; k++ ) if( !(b->diamond[k]->group & 0x01) ) ret += abs( b->diamond[k]->q );
#elif defined(MC_SITE)
    for( k = 0; k < N_site    ; k++ ) ret += abs( b->site[k]->q );
    for( k = 0; k < N_diamond ; k++ ) if(   b->diamond[k]->group & 0x01  ) ret += abs( b->diamond[k]->q );
#else
    for( k = 0; k < N_site    ; k++ ) ret += abs( b->site[k]->q    );
    for( k = 0; k < N_diamond ; k++ ) ret += abs( b->diamond[k]->q );
#endif
  };
  return ret;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* This will enter into the action multiplied by gamma_c */
static inline uint MC_killing_bond_gamma_c( Bond *b ){
  uint ret = 0;
  uchar k;
  for( k = 0; k < N_cube ; k++ ) ret += abs( b->cube[k]->q );
  return ret;
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static void MC_killing_bond_one_link( uchar i0, uint m, double gamma_s, double gamma_c ){
  uchar i, relax;
  SU2 Staples, g, g_new;
  double S_old_gauge, S_old_charge, S_new_gauge, S_new_charge, alpha;
  Bond *b;
#ifdef VERBOSE
  double help;
  uint tmp;
#endif
  b = &bond[link_index(i0,m)];
  Staples = MC_killing_get_staples( b );
  alpha = SU2_normalize( &Staples );
  g = SU2_mult( *(b->F), Staples);
#ifdef VERBOSE
  tmp = MC_killing_bond_gamma_s(b);
  MC_killing_monitor.mean_s += tmp;
  S_old_charge = gamma_s * tmp;
  //---------------------------------
  tmp = MC_killing_bond_gamma_c(b);
  MC_killing_monitor.mean_c += tmp;
  S_old_charge += gamma_c * tmp;
  //---------------------------------
  help = alpha * creal( g.alpha );
  MC_killing_monitor.mean_plaq += help;
  S_old_gauge = param->beta * help;
  //---------------------------------
#else
  S_old_charge  = gamma_s * MC_killing_bond_gamma_s(b);
  S_old_charge += gamma_c * MC_killing_bond_gamma_c(b);
  S_old_gauge   = param->beta * alpha * creal( g.alpha );
#endif
  for( i = 0; i < MC_KILLING_N_TRY ; i++ ){
    relax = i % 2;
    MC_killing_bond_data->F = *(b->F);
    if( relax ){ /* Overrelaxation */
      *(b->F) = SU2_conj( SU2_mult( Staples, g ) );
      g_new = SU2_conj( g );
      S_new_gauge = S_old_gauge;
    }else{   /* Random MC drift */
      *(b->F) = SU2_mult( *(b->F), SU2_random_matrix(MC_killing_delta) );
      g_new = SU2_mult( *(b->F), Staples);
      S_new_gauge   = param->beta * alpha * creal( g_new.alpha );
      MC_killing_monitor.try++ ;
    };
    S_new_charge = MC_killing_update_bond_variables( b, gamma_s, gamma_c );
    if( (S_old_charge - S_new_charge + S_new_gauge - S_old_gauge) > log(RND()) ){
      if( !relax ) MC_killing_monitor.accept++ ;
      S_old_charge = S_new_charge;
      S_old_gauge = S_new_gauge;
      g.alpha = g_new.alpha;
      g.beta = g_new.beta;
    }else{
      *(b->F) = MC_killing_bond_data->F;
      MC_killing_restore_bond_variables( b );
    };
  }; /* i (N_TRY) */
  SU2_normalize( b->F );
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static void MC_killing_bond( double gamma_s, double gamma_c ){
  uint m;
  uchar i0;
  double help;
  bzero( &MC_killing_monitor, sizeof(MC_killing_Monitor) );
  for( m = 0; m < param->ipw[param->D]; m++ )
    for( i0 = 0; i0 < param->D; i0++ )
      MC_killing_bond_one_link( i0, m, gamma_s, gamma_c );
#ifdef VERBOSE
  help = param->D * param->ipw[ param->D ];
  MC_killing_monitor.mean_plaq /= 2*(param->D-1)*help;
  MC_killing_monitor.mean_c /= N_cube * help;
  if( param->D == 3 ){
    MC_killing_monitor.mean_s /= N_site * help;
  }else{
    MC_killing_monitor.mean_s /= TOTAL_GAMMA_S * help;
  };
#endif
  help = ((double)MC_killing_monitor.accept) / ((double)MC_killing_monitor.try);
  if( help > 0.6 ) MC_killing_delta *= 1.25;
  if( help < 0.4 ) MC_killing_delta *= 0.75;
  if( MC_killing_delta > 2 ) MC_killing_delta = 2;
#ifdef VERBOSE
  printf("Accept=%d%%, delta=%.3f, beta=%.3f, gamma_s=%.3f, gamma_c=%.3f, "
	 "1/2<TrU_p>=%.5f, rho_s=%.5f, rho_c=%.5f,   Q_s=%.5f, Q_c=%.5f",
	 ((int)(help*100)), MC_killing_delta, param->beta, gamma_s, gamma_c,
	 MC_killing_monitor.mean_plaq, MC_killing_monitor.mean_s,
	 MC_killing_monitor.mean_c, density_site(), density_cube() );
  if( param->D == 4 )
    printf(", Q_diam_tot=%.5f, Q_diam_0=%.5f, Q_diam_1=%.5f",
	   density_diamond(), density_diamond0(), density_diamond1() );
  printf("\n");
  fflush( stdout );
#endif
};
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
int MC_killing( uint N_sweeps, double gamma_s, double gamma_c ){
  uint sweep;
  if( param->D != 3 && param->D != 4 ){
    error("Hey, man, this works only in D=3 and D=4!"); return -1;
  };
  if( !spinor || !bond || !plaq || !triangle || !square || !site || !cube ){
    error("Hey man, are you in hurry?! (structures undefined)"); return -1;
  };
  if( !MC_killing_bond_data ){
    MEMORY_CHECK(MC_killing_bond_data=(MC_killing_Bond_Data *)calloc(1,sizeof(MC_killing_Bond_Data)),-1);
#if defined(SPINOR_SECTION_PLUS)
    MEMORY_CHECK(MC_killing_bond_data->z2 = (double complex *)calloc(8*(param->D-1), sizeof(double complex)),-1);
#elif defined(SPINOR_SECTION_MINUS)
    MEMORY_CHECK(MC_killing_bond_data->z1 = (double complex *)calloc(8*(param->D-1), sizeof(double complex)),-1);
#else
    MEMORY_CHECK(MC_killing_bond_data->z1 = (double complex *)calloc(8*(param->D-1), sizeof(double complex)),-1);
    MEMORY_CHECK(MC_killing_bond_data->z2 = (double complex *)calloc(8*(param->D-1), sizeof(double complex)),-1);
#endif
    MEMORY_CHECK(MC_killing_bond_data->plaq_phase = (double *)calloc(2*(param->D-1), sizeof(double)),-1);
    if( N_triangle ) MEMORY_CHECK(MC_killing_bond_data->triangle_phase = (double *)calloc( N_triangle, sizeof(double)),-1);
    if( N_square   ) MEMORY_CHECK(MC_killing_bond_data->square_phase   = (double *)calloc( N_square,   sizeof(double)),-1);
    if( N_site     ) MEMORY_CHECK(MC_killing_bond_data->q_site         = (char *)  calloc( N_site,     sizeof(char)),-1);
    if( N_cube     ) MEMORY_CHECK(MC_killing_bond_data->q_cube         = (char *)  calloc( N_cube,     sizeof(char)),-1);
    if( N_diamond  ) MEMORY_CHECK(MC_killing_bond_data->q_diamond      = (char *)  calloc( N_diamond,  sizeof(char)),-1);
    /* Topology is not implemented, overwise you'll need these
       if( N_hypersite) MEMORY_CHECK(MC_killing_bond_data->q_hypersite    = (char *)  calloc( N_hypersite,sizeof(char)),-1);
       if( N_hypercube) MEMORY_CHECK(MC_killing_bond_data->q_hypercube    = (char *)  calloc( N_hypercube,sizeof(char)),-1);
    */
  };
  if( param->D == 4 ){
#if defined(MC_DIAMOND)
    printf("<<===  MC_killing [compiled with MC_DIAMOND], gamma_s kills 'duality' diamonds in D=4 ===>>\n");
#elif defined(MC_SITE)
    printf("<<===  MC_killing [compiled with MC_SITE], gamma_s kills group 1 diamonds and sites in D=4 ===>>\n");
#else
    printf("<<===  MC_killing [compiled with no special options], gamma_s kills all diamonds and sites in D=4 ===>>\n");
#endif
  };
  for( sweep = 0; sweep < N_sweeps; sweep++ ){
#ifdef VERBOSE
    printf("%3d of %3d:", sweep, N_sweeps );
    fflush( stdout );
#endif
    MC_killing_bond( gamma_s, gamma_c );
  };
  return 0;
};
