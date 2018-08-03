/*!
 * \file lattice.c
 *
 * \brief Lattice routines
 *
 *
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \author Pavel Buividovich, email: gbuividovich@gmail.com (implemented background magnetic field, chemical potential, SU(3) gauge group in 2008 - 2009)
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <memory.h>
#include <ctype.h>

#include "defs.h"
#include "types.h"
#include "lattice.h"
#include "panic.h"

#ifdef SU2
#include "su2math.h"
#endif

#ifdef SU3
#include "su3math.h"
#endif

#if (DIM == 4)
static int lat_size[DIM] = {LS, LS, LS, LT};
#elif (DIM == 3)
static int lat_size[DIM] = {LS, LS, LT};
#elif (DIM == 2)
static int lat_size[DIM] = {LS, LT};
#endif

t_mv * lat_mov;

/*!
 * \brief  Initialize movement array
 *
 * \param void
 *
 * \return void
 */
void lat_mov_init(void)
{
    int x[DIM], tmp[DIM];
    int mu;
    int idx_x, idx_moved;

    if (lat_mov != NULL)
        warn("lat_mov_init: lat_mov already initialized");
    else
        if ((lat_mov = malloc(sizeof(t_mv))) == NULL)
            panic("lat_mov_init: can't allocate memory for lat_mov");
    FORALL(x) {
        for (mu = 0; mu < DIM; mu++) {
            idx_x = INDEX_X(x);
            if (idx_x >= VOL) {
                fprintf(stderr, "x: (%i %i %i %i)\n",
                    x[0], x[1], x[2], x[3]);
                panic("lat_mov_init: idx_x(=%i) >= VOL(=%i)", idx_x, VOL);
            }
            memcpy(tmp, x, sizeof(x));
            X_FWD(tmp, mu);
            idx_moved = INDEX_X(tmp);
#ifdef DEBUG_LAT_MOV
            fprintf(stderr, "F: (%i %i %i %i) -> (%i) -> (%i %i %i %i)\n",
                    x[0], x[1], x[2], x[3], mu, tmp[0], tmp[1], tmp[2], tmp[3]);
#endif
            (*lat_mov)[idx_x][mu][FWD] = idx_moved;
            memcpy(tmp, x, sizeof(x));
            X_BWD(tmp, mu);
            idx_moved = INDEX_X(tmp);
#ifdef DEBUG_LAT_MOV
            fprintf(stderr, "B: (%i %i %i %i) -> (%i) -> (%i %i %i %i)\n",
                    x[0], x[1], x[2], x[3], mu, tmp[0], tmp[1], tmp[2], tmp[3]);
#endif
            (*lat_mov)[idx_x][mu][BWD] = idx_moved;
        }
    }
}

/*!
 * \brief create gauge configuration
 *
 * - allocate memory for gauge field configuration
 * - set all gauge matrices to ZERO
 *
 * \param void
 *
 * \return t_gauge_vector *
 */
t_gauge_vector * lat_gauge_create(void)
{
    int i, mu;
    t_gauge_vector *tmp;

    if ((tmp = malloc(sizeof(t_complex)*N_GMAT_COMP*VOL*DIM)) == NULL)
        panic("lat_gauge_create: can't allocate memory for gauge_vector");
    for(i = 0; i<VOL; i++)
     for(mu = 0; mu < DIM; mu++)
     {
      MATCALL(LOADZERO)((*tmp)[i][mu]);
     };
    return tmp;
}

/*!
 * \brief delete gauge configuration
 *
 * \param gv
 *
 * \return void
 */
void lat_gauge_destroy(t_gauge_vector * gv)
{
    if (gv != NULL)
        free(gv);
}

/*!
 * \brief info file structure
 *
 * Format of info file:
 * \verbatim
 * format #format
 * precision #precision
 * lattice #L_s #L_s #L_s #L_t
 * datafile #datafile
 * sun #sun
 * endian #endian
 * \endverbatim
 *
 * \#format is a format of datafile. Available formats:
 * <ul>
 *  <li> 1 - "Shadow" format.
 *  \f{eqnarray*}
 *  &&\dots <U_{j-1}><U_j><U_{j+1}> \dots \\
 *  &&j = (x_0 + L_s*(x_1 + L_s*(x_2 + L_s*(x_3))))*DIM + \mu \\
 *  &&x_0 \in [0, L_s-1] \\
 *  &&x_1 \in [0, L_s-1] \\
 *  &&x_2 \in [0, L_s-1] \\
 *  &&x_3 \in [0, L_t-1] \\
 *  &&\mu \in [0, DIM] \\
 *  &&j \in [0, 1, \dots, L_s^3*L_t*DIM] \\
 *  \f}
 *  \f$ U_j \f$ is a \f$ SU(2) \f$ matrix given by
 *  \f[ U_j = u_0 + i \sigma_i u_i, \f]
 *  where \f$ u_i \f$ is a float-point
 *  number of precision #precision.
 *  <li> 2 - "Bornyakov" format. Almost the same as "Shadow" format
 *  \f{eqnarray*}
 *  &&<SIZE>...<U_{j-1}><U_j><U_{j+1}>...<SIZE>,
 *  \f}
 *  where \f$<SIZE> = L_s^3*L_t*DIM*NMAT*sizeof(u_i)\f$.
 *  For this format only single (\#precision=4) precision is allowed.
 * </ul>
 *
 * For SU(3) (gi.format=1) (21.05.2009) is the format of Kanazawa's SU(3) mc generator:
 * 4 bytes of fortran header + Fortran array (DIM, NCOLORS, NCOLORS, NSITES) which stores all (NC x NC) components
 * of gauge matrix for all DIM links going from each of NSITES sites
 * (gi.format=2) is the format of Fortran array U1(Lt,L,L,L,4,3,2) saved to unformatted file where only the first two columns of the (3x3) gauge matrix are stored
 * This format is currently (21.05.2009) not tested - I had no configurations with such format
 * Note that before the calculation, the configuration is tested to check that all link matrices are in SU(3)
 * This excludes to a large extent reading the "wrong" configurations
 */
typedef struct t_gauge_info {
    int format;         //!< datafile format
    int precision;      //!< 4 - float, 8 - double
    int lattice[DIM];   //!< lattice size, should match LS, LT
    char datafile[256]; //!< filename of gauge configuration
    int sun_n;          //!< 2- su2, 3 -su3.
    int endian;         //!< l - little, b - big. !!!Unsupported for now!!!
} t_gauge_info;

/*!
 * \brief Read gauge configuration parameters
 *
 * \param gi
 * \param info_fname
 *
 * \return void
 */
void lat_info_read(t_gauge_info * gi, char * info_fname)
{
    FILE * info = NULL;
    int i;
    char ch;

    if ((info = fopen(info_fname, "rb")) == NULL)
        panic("lat_gauge_load: can't open info file(=%s)", info_fname);
    /* format */
    fscanf(info, "format %i\n", &(gi->format));
#ifdef SU2
    if (gi->format != 1 && gi->format != 2 && gi->format != 3)
        panic("lat_info_read: wrong info version(=%i)", gi->format);
    /* precision */
    fscanf(info, "precision %i\n", &(gi->precision));
    if (gi->format == 1)
    {
        if (!(gi->precision == 4 || gi->precision == 8))
            panic("lat_info_read: unsupported precision(=%i) for this format(=%i)",
                    gi->precision, gi->format);
    } else if (gi->format == 2)
    {
        if (!(gi->precision == 4))
            panic("lat_info_read: unsupported precision(=%i) for this format(=%i)",
                    gi->precision, gi->format);
    }
#endif
#ifdef SU3
    if(gi->format!=1  && gi->format != 2)
     panic("lat_info_read: wrong info version(=%i) for SU(%i) group", gi->format, NCOLORS);
    /* precision */
    fscanf(info, "precision %i\n", &(gi->precision));
    if (gi->format == 1 || gi->format == 2)
    {
     if (!(gi->precision == 4 || gi->precision == 8))
      panic("lat_info_read: unsupported precision(=%i) for this format(=%i)", gi->precision, gi->format);
    };
#endif
    /* lattice size */
    fscanf(info, "lattice %i", &(gi->lattice[0]));
    if (gi->lattice[0] != LS)
        panic("lat_info_read: lattice Ssize(=%i) != LS(=%i)",
                gi->lattice[0], LS);
    for (i = 1; i < DIM-1; i++) {
        fscanf(info, "%i", &(gi->lattice[i]));
        if (gi->lattice[i] != LS)
            panic("lat_info_read: lattice Ssize(=%i) != LS(=%i)",
                    gi->lattice[i], LS);
    }
    fscanf(info, "%i\n", &(gi->lattice[DIM-1]));
    if (gi->lattice[DIM-1] != LT)
        panic("lat_info_read: lattice Tsize(=%i) != LT(=%i)",
                gi->lattice[DIM-1], LT);
    /* datafile */
    while(!isspace(ch = (char) fgetc(info)));
    memset(gi->datafile, 0, 256);
    while(isspace(ch = (char) fgetc(info)));
    ungetc(ch, info);
    i = 0;
    do {
        ch = (char) fgetc(info);
        gi->datafile[i++] = ch;
    } while (!isspace(ch));
    gi->datafile[i-1] = 0;
    /* sun */

    fscanf(info, "\nsun %i\n", &(gi->sun_n));
    if(gi->sun_n != NCOLORS)
     panic("lat_info_read: unsupported gauge group sun(=%i), NCOLORS = %i", gi->sun_n, NCOLORS);
    /* endian */
    fscanf(info, "endian %i\n", &(gi->endian));
    fclose(info);
}

/*!
 * \brief load gauge configuration
 *
 * - read info file (description) of configuration
 * - read gauge configuration
 *
 * \param gv t_gauge_vector
 * \param info_fname ifno filename
 *
 * \return void
 */
void lat_gauge_load(t_gauge_vector * gv, char * info_fname)
{
    t_gauge_info gi;
    FILE * gfile = NULL;
    int i, mu, osize;
    void *u;
//    float* u1,u2;

    if (gv == NULL)
        panic("lat_gauge_load: gv == NULL");
    lat_info_read(&gi, info_fname);
    if ((gfile = fopen(gi.datafile, "rb")) == NULL)
        panic("lat_gauge_load: can't open gauge file(=%s) for reading", gi.datafile);
#ifdef SU2
    osize = gi.precision*N_GMAT_COMP*2; /* 2 - because of complex components */
    if ((u = malloc(osize)) == NULL)
        panic("lat_gauge_load: can't allocate memory for u");
    if (gi.format == 1) {
    /* shadow format */
        if (gi.precision == 4)
            for (i = 0; i < VOL; i++)
                for (mu = 0; mu < DIM; mu++) {
                    fread(u, osize, 1, gfile);
                    if (gi.sun_n == NCOLORS) {
                        (*gv)[i][mu][0] = *((float *)u)   + I*(*((float *)u+3));
                        (*gv)[i][mu][1] = *((float *)u+2) + I*(*((float *)u+1));
                    }
                }
        else if (gi.precision == 8)
            for (i = 0; i < VOL; i++)
                for (mu = 0; mu < DIM; mu++) {
                    fread(u, osize, 1, gfile);
                    if (gi.sun_n == NCOLORS) {
                        (*gv)[i][mu][0] = *((double *)u)+I*(*((double *)u+3));
                        (*gv)[i][mu][1] = *((double *)u+2)+I*(*((double *)u+1));
                    }
                }
    } else if (gi.format == 2) {
    /* bornyakov format */
        if (gi.precision != 4)
            panic("%s: unsupported precision(=%i) for format(=%i)", __func__,
                    gi.precision, gi.format);
        if (gi.sun_n != NCOLORS)
            panic("lat_gauge_load: unsupported gauge group(=%i) for format(=%i)",
                    gi.sun_n, gi.format);
        fread(u, 4, 1, gfile); /* read 4 bytes of fortran header */
        for (i = 0; i < VOL; i++)
            for (mu = 0; mu < DIM; mu++) {
                fread(u, osize, 1, gfile);
                (*gv)[i][mu][0] = *((float *)u)+I*(*((float *)u+3));
                (*gv)[i][mu][1] = *((float *)u+2)+I*(*((float *)u+1));
            }
    } else if (gi.format == 3) {
        /* bornyakov U(1) (phase format) */
        if (gi.precision != 4)
            panic("%s: unsupported precision(=%i) for format(=%i)",
                    __func__, gi.precision, gi.format);
        if (gi.sun_n != NCOLORS)
            panic("%s: unsupported gauge group(=%i) for format(=%i)",
                    __func__, gi.sun_n, gi.format);
        fread(u, 4, 1, gfile); /* read 4 bytes of fortran header */
        for (i = 0; i < VOL; i++)
            for (mu = 0; mu < DIM; mu++) {
                fread(u, sizeof(float), 1, gfile);
                /* cos(phi) + I sin(phi) */
                (*gv)[i][mu][0] = cos(*((float *)u)) + I*sin(*((float *)u));
                /* 0.0 */
                (*gv)[i][mu][1] = 0.0;
            }
    }
#endif
#ifdef SU3
 u = NULL;
 if(gi.format == 1)
 {
  int id1, id2;
  osize = gi.precision*2; /* 2 - because of complex components */
  if((u = malloc(osize)) == NULL)
   panic("lat_gauge_load: can't allocate memory for u");
//  fread(u, 4, 1, gfile); /* read 4 bytes of fortran header */
//  for(mu = 0; mu < DIM; mu++)
//   for(id1 = 0; id1 < NCOLORS; id1++)
//    for(id2 = 0; id2 < NCOLORS; id2++)
//     for(i=0; i<VOL; i++)
    for(i=0;i<VOL;i++)
     for(mu=0;mu<DIM;mu++)
      for(id1=0;id1<NCOLORS;id1++)
       for(id2=0;id2<NCOLORS;id2++)
     {
      fread(u, osize, 1, gfile);
      if(gi.precision==4)
       (*gv)[i][mu][id2*NCOLORS + id1] = *((float *)u)+I*(*((float *)u+1));
      if(gi.precision==8)
       (*gv)[i][mu][id2*NCOLORS + id1] = *((double *)u) + I*(*((double *)u+1));
     };
 };
 if(gi.format == 2)
 {
  int id1, id2, xt;
  osize = gi.precision*2; /* 2 - because of complex components */
  if((u = malloc(osize)) == NULL)
   panic("lat_gauge_load: can't allocate memory for u");
//  fread(u, 4, 1, gfile); /* read 4 bytes of fortran header */
  /* U1(Lt,L,L,L,4,3,2) */
//  for(id1 = 0; id1 < 2; id1++)
//   for(id2 = 0; id2 < 3; id2++)
  for(xt=0; xt<LT; xt++)
   for(i=0;i<LS3;i++)
    for(mu = 0; mu < DIM; mu++)
//     for(i=0; i<LS3; i++)
//      for(xt=0; xt<LT; xt++)
      for(id1=0;id1<2;id1++)
       for(id2=0;id2<3;id2++)
      {
       int nmu = (mu==0)? 3: mu-1;
       fread(u, osize, 1, gfile);
       if(gi.precision==4)
        (*gv)[LS3*xt + i][nmu][id2*NCOLORS + id1] = *((float *)u)+I*(*((float *)u+1));
       if(gi.precision==8)
        (*gv)[LS3*xt + i][nmu][id2*NCOLORS + id1] = *((double *)u) + I*(*((double *)u+1));
      };
   /* Restoring the last column of link matrices */
  for(i=0; i<VOL; i++)
   for(mu = 0; mu < DIM; mu++)
   {
    (*gv)[i][mu][2] = CCALL(conj)((*gv)[i][mu][3]*(*gv)[i][mu][7] - (*gv)[i][mu][6]*(*gv)[i][mu][4]);
    (*gv)[i][mu][5] = CCALL(conj)((*gv)[i][mu][6]*(*gv)[i][mu][1] - (*gv)[i][mu][0]*(*gv)[i][mu][7]);
    (*gv)[i][mu][8] = CCALL(conj)((*gv)[i][mu][0]*(*gv)[i][mu][4] - (*gv)[i][mu][3]*(*gv)[i][mu][1]);
   };
 };

#endif
    free(u);
    fclose(gfile);
}

#define EPS_UNITARITY       10e-3   //!< allowed nonunitarity: (det(A)-1.0) < EPS_UNITARITY

/*!
 * \brief Check unitarity of each gauge matrix in configuration
 *
 * \param gv t_gauge_vector
 *
 * \return void
 */
void lat_gauge_check_unit(t_gauge_vector * gv)
{
 int j, mu, k, l;
 t_real sum;
 t_gmat t;

 for(j = 0; j<VOL; j++)
 {
  for(mu = 0; mu < DIM; mu++)
  {
   memcpy(&t, (*gv)[j][mu], sizeof(t_gmat));
#ifdef SU2
#ifdef DEBUG_CHECK_UNIT
  MATCALL(print_mat)(t);
#endif
   sum = MATCALL(det_A)(t);
   if((sum - 1.0) > EPS_UNITARITY)
   {
    warn("%s: matrix not unitary (j = %i, mu = %i, sum = %f)\n", __func__, j, mu, sum);
  MATCALL(print_mat)(t);
   };
#endif
#ifdef SU3
   int tmp_i, tmp_j, tmp_k;
   t_gmat td, tdt;
   memcpy(&td, t, sizeof(t_gmat));
   MATCALL(hconj_A)(td);
   MATCALL(C_eq_A_mul_B)(tdt, td, t);
   sum = 0.0;
   for(k=0; k<NCOLORS; k++)
    for(l=0; l<NCOLORS; l++)
     sum += CCALL(cabs)(tdt[k*NCOLORS + l] - ((k==l)? 1.0 : 0.0) );
   if(sum > EPS_UNITARITY)
   {
    warn("%s: matrix not unitary (j = %i, mu = %i, sum = %f)\n", __func__, j, mu, sum);
    MATCALL(print_mat)(t);
   };
   sum = CCALL(cabs)(MATCALL(det_A)(t));
   if(fabs(sum - 1.0) > EPS_UNITARITY)
   {
    warn("%s: Absolute value of the link determinant not equal to one (j = %i, mu = %i, det = %f)\n", __func__, j, mu, sum);
    MATCALL(print_mat)(t);
   };
#endif
  };
 };
}

/*!
 * \brief calculate gauge action
 *
 * \f$ A = 1.0 - \frac{1}{N_c}\sum Re(Tr(U_{plaq})) \f$
 *
 * \param gv
 *
 * \return t_real
 */
t_real lat_gauge_action(t_gauge_vector * gv)
{
    int idx1, idx2, idx3;
#ifdef SU3
    int tmp_i, tmp_j, tmp_k;
#endif
    int mu, nu;
    t_gmat m3, m4;
    t_real sum = 0.0,  A = 0.0;

    for (idx1 = 0; idx1 < VOL; idx1++)
    {
     for (mu = 0; mu < DIM; mu++)
     {
      for (nu = mu+1; nu < DIM; nu++)
      {
       /*
        *idx3---->+
        *   ^     ^
        *   |     |
            *idx1---->idx2
        *
            * U_p = U[idx2][nu] * U[idx3][mu]^+ * U[idx1][nu]^+ * U[idx1][mu]
        */
       idx2 = (*lat_mov)[idx1][mu][FWD];
       idx3 = (*lat_mov)[idx1][nu][FWD];
       MATCALL(C_eq_A_mul_Bd)(m3, (*gv)[idx2][nu], (*gv)[idx3][mu]);
       MATCALL(C_eq_A_mul_Bd)(m4, m3, (*gv)[idx1][nu]);
       sum += MATCALL(re_tr_A_mul_B)(m4, (*gv)[idx1][mu]);
      };
     };
    };
    A = 1.0 - sum/(N_PLAQS*NCOLORS);
    return A;
}

/*!
 * \brief Load identity configuration
 *
 * \param gv t_gauge_vector
 *
 * \return void
 */
void lat_gauge_identity(t_gauge_vector * gv)
{
    int i, mu;
    for(i = 0; i < VOL; i++)
     for(mu = 0; mu < DIM; mu++)
     {
      MATCALL(LOADIDENTITY)((*gv)[i][mu]);
     };
}

void   xcoord(int x, int* xc)
{
 int q = x;
 xc[3] = q/LS3;
 q = q%LS3;
 xc[2] = q/LS2;
 q = q%LS2;
 xc[1] = q/LS;
 q = q%LS;
 xc[0] = q;
}

/*{{{*/
/*!
 * \brief Calculate "taxi driver distance" between two points
 *
 * \param x
 * \param y
 *
 * \return int \f$ = \sum_{\mu} |x_\mu - y_\mu| \f$
 *//*}}}*/
#define xmax(_x, _y, _t)                                            \
    if ( (_x) < (_y) ) { (_t) = (_y); (_y) = (_x); (_x) = (_t);}
#define min(_a, _b) ((_a) < (_b) ? (_a) : (_b))
int taxi_distance(int* x, int* y)  /*{{{*/
{
    int i, tmp;
    int tx, ty;
    int sum = 0;

    for (i = 0; i < DIM-1; i++)  {
        tx = x[i];
        ty = y[i];
        xmax(tx, ty, tmp);
        sum += abs(min(tx - ty , LS - tx + ty));
    }
    tx = x[DIM-1];
    ty = y[DIM-1];
    xmax(tx, ty, tmp);
    sum += abs(min(tx - ty, LT - tx + ty));

    return sum;
}/*}}}*/
#undef min
#undef xmax

t_real distance4D(int x1, int x2)
{
 int xc1[DIM], xc2[DIM], mu, mdist;
 xcoord(x1, xc1);
 xcoord(x2, xc2);
 t_real dist = 0.0;
 for(mu=0; mu<DIM; mu++)
 {
  mdist = (xc1[mu] - xc2[mu]);
  if(abs(mdist - lat_size[mu])<abs(mdist))
   mdist = mdist - lat_size[mu];
  if(abs(mdist + lat_size[mu])<abs(mdist))
   mdist = mdist + lat_size[mu];
  dist += (t_real)(mdist*mdist);
 };
 return sqrt(dist);
}

void lat_gauge_aperiodic(t_gauge_vector* gv)
{
 int s[DIM];
 int idx, mu, ic;
 t_gmat m;
 s[3] = LT-1;
 for(s[0] = 0; s[0] < LS; s[0]++)
  for(s[1] = 0; s[1] < LS; s[1]++)
   for(s[2] = 0; s[2] < LS; s[2]++)
   {
    idx = INDEX_X(s);
    mu = 3;
    memcpy(&m, (*gv)[idx][mu], sizeof(t_gmat));
    for(ic=0; ic<N_GMAT_COMP; ic++)
     m[ic] *= -1.0;
    memcpy((*gv)[idx][mu], &m, sizeof(t_gmat));
   };
}

#ifdef EFIELD
void init_efield(int E, int H)
{
 int s[DIM];
 int idx;
 t_real FE = 2*M_PI*(t_real)E/(t_real)(LS*LT);
 t_real FH = 2*M_PI*(t_real)H/(t_real)(LS*LS);
 FORALL_4D(s)
 {
  idx = INDEX_X_4D(s);
  ev[idx][3] = cos(0.5*FE*(t_real)s[0]) - I*sin(0.5*FE*(t_real)s[0]);
  ev[idx][2] = 1.0;
  ev[idx][1] = cos(0.5*FH*(t_real)s[0]) - I*sin(0.5*FH*(t_real)s[0]);
  ev[idx][0] = cos(0.5*FH*(t_real)s[1] +  0.5*FE*(t_real)s[3]) + I*sin(0.5*FH*(t_real)s[1] +  0.5*FE*(t_real)s[3]);
 };
 s[3] = LT-1;
 for(s[0] = 0; s[0] < LS; s[0]++)
  for(s[1] = 0; s[1] < LS; s[1]++)
   for(s[2] = 0; s[2] < LS; s[2]++)
   {
    idx = INDEX_X_4D(s);
    ev[idx][3] *= cos(0.5*FE*(t_real)LT*(t_real)s[0]) - I*sin(0.5*FE*(t_real)LT*(t_real)s[0]);
   };
 s[1] = LS-1;
 for(s[0] = 0; s[0] < LS; s[0]++)
  for(s[2] = 0; s[2] < LS; s[2]++)
   for(s[3] = 0; s[3] < LT; s[3]++)
   {
    idx = INDEX_X_4D(s);
    ev[idx][1] *= cos(0.5*FH*(t_real)LS*(t_real)s[0]) - I*sin(0.5*FH*(t_real)LS*(t_real)s[0]);
   };
 s[0] = LS-1;
 for(s[1] = 0; s[1] < LS; s[1]++)
  for(s[2] = 0; s[2] < LS; s[2]++)
   for(s[3] = 0; s[3] < LT; s[3]++)
   {
    idx = INDEX_X_4D(s);
    ev[idx][0] *= cos(0.5*FE*(t_real)LS*(t_real)s[3] + 0.5*FH*(t_real)LS*(t_real)s[1]) + I*sin(0.5*FE*(t_real)LS*(t_real)s[3] + 0.5*FH*(t_real)LS*(t_real)s[1]);
   };
}
#endif

/*!
 * \brief Create one instanton
 *
 * \param gv
 *
 * \return
 */
void lat_gauge_instanton(t_gauge_vector * gv, t_real rho, t_real* x0)
{
    int x[DIM];
    t_real phi[DIM], phi0[DIM];
    t_real delta[DIM];
    t_real a[DIM][3];
    t_real eta[3];
    int mu;
    t_real num, den, norm, den0, num0;
    int j, k;
    t_gmat m;
    int idx;

    FORALL(x) {
        idx = INDEX_X(x);
        for (mu = 0; mu < DIM; mu++) {
            memcpy(&m , (*gv)[idx][mu], sizeof(t_gmat));
            /* define Phi_\mu(x) */
            for (j = 0; j < DIM; j++) {
                num = 0.0;
                den = 0.0;
                for (k = 0; k < DIM; k++) {
                    if (k == j)
                        den += (x[k]-x0[k])*(x[k]-x0[k]);
                    else {
                        num += (x[k]-x0[k])*(x[k]-x0[k]);
                        den += (x[k]-x0[k])*(x[k]-x0[k]);
                    }
                }
                num0 = sqrt(num);
                den0 = den;
                den0 += (x[j]-x0[j]);
                phi0[j] = atan2(num0, den0)/num0;

                num += rho*rho;
                num = sqrt(num);
                den += rho*rho + (x[j]-x0[j]);
                phi[j] = atan2(num, den)/num;
            }
            /* define Delta_\mu */
            for (j = 0; j < DIM; j++)
                delta[j]  = x[j] - x0[j];
            a[0][0] = -delta[3];
            a[0][1] = delta[2];
            a[0][2] = -delta[1];

            a[1][0] = -delta[2];
            a[1][1] = -delta[3];
            a[1][2] = delta[0];

            a[2][0] = delta[1];
            a[2][1] = -delta[0];
            a[2][2] = -delta[3];

            a[3][0] = delta[0];
            a[3][1] = delta[1];
            a[3][2] = delta[2];
            for (k = 0; k < 3; k++) {
                eta[k] = 0.0;
                for (j = 0; j < DIM; j++)
                    eta[k] += a[j][k]*(phi0[j]-phi[j]);
            }
            norm = 0.0;
            for (k = 0; k < 3; k++)
                norm += eta[k]*eta[k];
            norm = sqrt(norm);
            if(norm != 0.0)
            {
             for (k = 0; k < 3; k++)
              eta[k] /= norm;
#ifdef SU2
             m[0] = cos(norm) + I*eta[2]*sin(norm);
             m[1] = sin(norm)*eta[1] + I*eta[0]*sin(norm);
#endif
#ifdef SU3
             m[0] = cos(norm) + I*eta[2]*sin(norm);
             m[1] = sin(norm)*eta[1] + I*eta[0]*sin(norm);
             m[3] = - sin(norm)*eta[1] + I*eta[0]*sin(norm);
             m[4] = cos(norm) - I*eta[2]*sin(norm);
             m[2] = 0.0 + 0.0*I;
             m[5] = 0.0 + 0.0*I;
             m[8] = 1.0 + 0.0*I;
             m[6] = 0.0 + 0.0*I;
             m[7] = 0.0 + 0.0*I;
#endif
            }
            else
            {
             MATCALL(LOADIDENTITY)(m);
            }
            memcpy((*gv)[idx][mu], &m, sizeof(t_gmat));
        }
    }
}
