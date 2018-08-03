/*{{{*/
/*!
 * \file propagators.c
 *
 * \brief Invert Overlap Dirac operator
 *
 *
 * $Id$
 *
 * \author Pavel Buividovich, email: gbuividovich@gmail.com (implemented background magnetic field, chemical potential, SU(3) gauge group in 2008 - 2009)
 */
/*}}}*/
#define _XOPEN_SOURCE   600
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <time.h>
#include "defs.h"
#include "types.h"
#include "mt19937ar.h"
#include "lattice.h"
#include "arpack.h"
#include "dirac.h"
#include "panic.h"
#include "timing.h"
#include "linal.h"
#include "minmax.h"
#include "polacc.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#undef DEFLATE_PROP //TODO: in Makefile.var
#define GAUSSIAN_SOURCE

#ifdef SU2
#include "su2math.h"
#endif

#ifdef SU3
#include "su3math.h"
#endif

#define OP_TYPE_OV                      (0)
#define OP_TYPE_OV_CC                   (1)
#define OP_TYPE_MAX                     (2)

#define DEFAULT_DELTA                   (1e-8)

#define DEFAULT_TOL                     (1e-6)
#define IMAX                            (1000)

#ifdef MU
#define DEFAULT_DEG                     (50)
#else
#define DEFAULT_DEG                     (20)
#endif

#define DEFAULT_XB                      (1.0)
#define DEFAULT_NCV_FACTOR_NOM          (5)

#define MAX_ZERO_EV  (0.000001)

void printhelp(void) __attribute__ ((noreturn));

/*{{{*/
/*!
 * \brief Print help/usage in case of wrong option
 *
 * \param void
 *
 * \return void
 */
/*}}}*/
void printhelp(void) /*{{{*/
{

    fprintf(stdout, "Overlap Dirac Operator Spectrum\n");
    fprintf(stdout, "Usage:\n");
    fprintf(stdout, "-o fname\t -- outputfilename\n");
    fprintf(stdout, "[-f]\t\t -- free operator\n");
    fprintf(stdout, "[-I]\t\t -- instanton\n");
    fprintf(stdout, "[-i]\t\t -- info of gauge configuration file\n");
#ifndef MU
    fprintf(stdout, "[-t]\t\t -- 0 - D_ov, 1 - D_ov_cc\n");
#endif
    fprintf(stdout, "[-a]\t\t -- aperiodic boundary conditions\n");
    fprintf(stdout, "-O fname\t -- output file.\n");
    fprintf(stdout, "[-R]\t -- calculate left and right eigenvectors\n");
    fprintf(stdout, "-r rho\t\t -- mass parameter rho. rho > 0\n");
    fprintf(stdout, "[-m num]\t -- overlap mass term\n");
    fprintf(stdout, "[-k num]\t -- anisotropy factor, the ratio of kappa_t to kappa_s\n");
    fprintf(stdout, "[-d num]\t -- delta. Default is %12.10f.\n", DEFAULT_DELTA);
    fprintf(stdout, "[-T num]\t -- tolerance. Default is %12.10f.\n", DEFAULT_TOL);
#ifdef EFIELD
    fprintf(stdout, "-E num\t -- external electric field. Default is zero.\n");
    fprintf(stdout, "-H num\t -- external magnetic field. Default is zero.\n");
#endif
#ifdef MU
    fprintf(stdout, "-M mu\t --  chemical potential mu.\n");
    fprintf(stdout, "-[L] num\t -- tolerance for the calculation of the sign of the inner matrix.\n");
#endif
    fprintf(stdout, "[-p]\t\t -- DO polynomial acceleration\n");
    fprintf(stdout, "[-X num]\t -- polynomial acceleration xb. [DEFAULT:  %f]\n", DEFAULT_XB);
    fprintf(stdout, "[-P num]\t -- degree of polynomial. [DEFAULT: %i]\n", DEFAULT_DEG);
    exit(EXIT_SUCCESS);
}/*}}}*/

int sort_ev(t_complex* evecs, t_complex* evals, int nev) //TODO: declare in linal.h or dirac.h
{
 t_cds_vector* tmp;
 if(posix_memalign((void*)&tmp, 16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for tmp (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 int i, k, ki;
 for(i=0; i<nev; i++)
 {
  t_real minlambda = 100000000;
  ki = 0;
  for(k=i; k<nev; k++)
   if(CCALL(cabs)(evals[k]) < minlambda)
   {
    minlambda = CCALL(cabs)(evals[k]);
    ki = k;
   };
  k = ki;
  memcpy(                           tmp, &(evecs[i*VOL*NCOLORS*NDIRAC]), VOL*NCOLORS*NDIRAC*sizeof(t_complex));
  memcpy(&(evecs[i*VOL*NCOLORS*NDIRAC]), &(evecs[k*VOL*NCOLORS*NDIRAC]), VOL*NCOLORS*NDIRAC*sizeof(t_complex));
  memcpy(&(evecs[k*VOL*NCOLORS*NDIRAC]),                            tmp, VOL*NCOLORS*NDIRAC*sizeof(t_complex));
  t_complex tmpval = evals[i];
  evals[i] = evals[k];
  evals[k] = tmpval;
 };
 free(tmp);

 int nzero = 0;
 while(CCALL(cabs)(evals[nzero])<MAX_ZERO_EV)
  nzero++;

 return nzero;
}

/*{{{*/
/*!
 * \brief  main routine
 *
 * \param argc
 * \param argv[]
 *
 * \return int
 */
/*}}}*/
int main (int argc, char * argv[]) /*{{{*/
{
    srand(time(NULL));
    t_gauge_vector * gv;
    int c;
    int i, j, a, id, x, x0, ret;
    char *iname = NULL;
    char *oname = NULL;
    char *OVname = NULL;
    t_real kappa = -1.0, rho = -1.0, dist;
    t_real r0 = 1.0; // Radius of a Gaussian source, TODO: input by user
    int freewd = 0;
    t_real action;
    int optype = -1;
    t_ov_data_mm ov_data;
    int snev = 30; /* No of smalles ev to deflate TODO: control by users */
    int polacc = 0;
    t_real delta = DEFAULT_DELTA, mass = 0.0, tol = DEFAULT_TOL;
#ifndef MU
    t_real eps;
    t_complex* tmpevals;
    t_complex  zb;
    double lambda1 = 0.0;
    int polacc_deg = DEFAULT_DEG;
    t_real polacc_xb = DEFAULT_XB;
    t_polacc_data polacc_data;
    t_op_ev evd;
#endif
    int res;
    int instanton = 0;
#ifdef EFIELD
    int E = 0, H = 0;
#endif
    anisotropy_factor = 1.0;
#ifdef MU
    Mu = 0.0;
    t_double_real sign_tol = 0.0000001;
    int deg = 200;
#endif
    t_real instR = 2.0;
    t_real ic[DIM] = {0.0, 0.0, 0.0, 0.0};
    int aperiodic = 0;
    int ncv_factor_nom = DEFAULT_NCV_FACTOR_NOM, ncv_factor_denom = 2;

    while ((c = getopt(argc, argv, "pr:E:H:M:L:o:ft:i:O:d:m:T:X:P:Iak:")) != -1) /*{{{*/
    {
        switch (c) {
        case 'a':
            aperiodic = 1;
            break;
        case 'k':
            anisotropy_factor = atof(optarg);
            break;
        case 'I':
            instanton = 1;
            break;
#ifndef MU
        case 'P':
            polacc_deg = atoi(optarg);
            if (polacc_deg <= 0) {
                fprintf(stdout, "Polynomial acceleration polynomial degree < 0\n");
                printhelp();
            }
            break;
        case 'X':
            polacc_xb = atof(optarg);
            /* TODO: not only massless check */
            if (polacc_xb >= 2.0 || polacc_xb <= 0.0 ) {
                fprintf(stdout, "Polynomial acceleration xb should be in (0;2.0) range\n");
                printhelp();
            }
            break;
#endif
        case 'p':
            polacc = 1;
            break;
        case 'm':
            mass = atof(optarg);
            break;
        case 'd':
            delta = atof(optarg);
            break;
        case 'o':
            oname = malloc(sizeof(char)*(strlen(optarg)+1));
            memset(oname, 0, sizeof(char)*(strlen(optarg)+1));
            strcpy(oname, optarg);
            break;
        case 'r':
            rho = atof(optarg);
            break;
#ifdef EFIELD
        case 'E':
            E = atoi(optarg);
            break;
        case 'H':
            H = atoi(optarg);
            break;
#endif
#ifdef MU
        case 'M':
            Mu = atof(optarg);
            break;
        case 'L':
            sign_tol = atof(optarg);
            break;
        case 'P':
            deg = atoi(optarg);
            break;
#endif
        case 'f':
            freewd = 1;
#ifdef OV_DIRAC_ORTHO_PROJ /*{{{*/
            fprintf(stdout, "Projection doesn't work with degenerate spectrum\n");
            fprintf(stdout, "Turn Projection off in Makefile.var\n");
            exit(EXIT_FAILURE);
#endif
            break;
        case 'T':
            tol = atof(optarg);
            break;
        case 'i':
            iname = malloc(sizeof(char)*strlen(optarg));
            strcpy(iname, optarg);
            break;
        case 't':
            optype = atoi(optarg);
            if (optype >= OP_TYPE_MAX || optype < 0) {
                fprintf(stdout, "Unsupported operator type(=%i)\n", optype);
                printhelp();
            }
            switch (optype) {
                case OP_TYPE_OV:
                    evd.op = &ov_dirac_mm;
                    evd.data = (t_ov_data_mm*)&ov_data;
                    break;
                case OP_TYPE_OV_CC:
                    evd.op = &ov_dirac_mm_cc;
                    evd.data = (t_ov_data_mm*)&ov_data;
                    break;
            }
            break;
        case 'O':
            OVname = malloc(sizeof(char)*(strlen(optarg)+1));
            memset(OVname, 0, sizeof(char)*(strlen(optarg)+1));
            strcpy(OVname, optarg);
            break;
        case '?':
            if (isprint(optopt))
                fprintf(stderr, "Unknown option '-%c'.\n", optopt);
            else
                fprintf(stderr,
                        "Unknown option character '\\x%x'.\n", optopt);
            printhelp();
        default:
            printhelp();
        }
    }/*}}}*/
#ifdef MU
    if(polacc==1)
    {
     fprintf(stdout, "\n Polinomial acceleration does not work at nonzero mu, turned OFF automatically \n");
     polacc = 0;
    };
#endif
    if (!instanton && !freewd && iname == NULL)
        printhelp();
    else if ((instanton || freewd) && iname != NULL)
        printhelp();
    if (!oname || optype < 0)
        printhelp();
    if (OVname == NULL) {
        fprintf(stdout, "You should specify ovdata file\n");
        printhelp();
    }
    /* evd initialization */
    evd.evals = NULL;
    evd.evecs = NULL;
    res = read_ev(OVname, &evd, SM_SECTION);
    if ((res & E_READ_CANTREAD) || (res & E_READ_CANTOPEN) || (res & E_READ_SECNOTHERE)) {
        fprintf(stdout, "Some errors! res = %i\n", res);
        printhelp();
    };
    ov_data.n_proj  = evd.nev;
    kappa = rho2kappa(rho);
    ov_data.kappa = kappa;
    ov_data.mass = mass;
#ifdef MU
    ov_data.sign_tol = sign_tol;
    ov_data.deg      = deg;
    ov_data.evals  = evd.evals;
    ov_data.revecs = evd.evecs;
    evd.evals = NULL;
    evd.evecs = NULL;
#else
    if ((tmpevals = malloc(sizeof(t_complex)*evd.nev)) == NULL)
        panic("%s: can't allocate memory for tmpevals (%i Kb)", __func__, sizeof(t_complex)*evd.nev/1024);
    memcpy(tmpevals, evd.evals, sizeof(t_complex)*evd.nev);
    if (evd.nev >= 2)
        qsort(tmpevals, evd.nev, sizeof(t_complex), &cmplx_cmp_mod);
#ifdef OV_DIRAC_ORTHO_PROJ /*{{{*/
    /* init projectors */
    ov_data.evecs = evd.evecs;
    ov_data.evals = evd.evals;
    ov_data.l_sq_min = CCALL(creal)(tmpevals[ov_data.n_proj-1]*CCALL(conj)(tmpevals[ov_data.n_proj-1]));
    /* enlarge spectral window to avoid nonconvergance of Cheb. series */
    ov_data.l_sq_min -= ov_data.l_sq_min*SPEC_MULT;
    /* */
    evd.evals = NULL;
    evd.evecs = NULL; /*}}}*/
#else /*{{{*/
    ov_data.l_sq_min = CCALL(creal)(tmpevals[0])*CCALL(creal)(tmpevals[0]);
    /* enlarge spectral window to avoid nonconvergance of Cheb. series */
    ov_data.l_sq_min -= SPEC_MULT*ov_data.l_sq_min;
    if (evd.evals != NULL) {
        free(evd.evals);
        evd.evals = NULL;
    }
    if (evd.evecs != NULL) {
        free(evd.evecs);
        evd.evecs = NULL;
    }
#endif/*}}}*/
    free(tmpevals);
#endif
    res = read_ev(OVname, &evd, LM_SECTION);
    if ((res & E_READ_CANTREAD) || (res & E_READ_CANTOPEN) || (res & E_READ_SECNOTHERE)) {
        fprintf(stdout, "Some errors! res = %i\n", res);
        printhelp();
    }
#ifdef MU
    /* Sorting eigenvalues so that left/right agree */
    int k, ki = 0, exchg = 0;
    t_cds_vector* tmp;
    if(posix_memalign((void*)&tmp, 16, sizeof(t_cds_vector)))
     panic("%s: can't allocate memory for tmp (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
    fprintf(stdout, "\n Sorting left evals/evecs ... \n");
    fflush(stdout);
    for(i=0; i<ov_data.n_proj-1; i++)
    {
     exchg = 0;
     t_real mindist = 100000000;
     for(k=i; k<evd.nev; k++)
     {
      t_real dist = CCALL(cabs)(evd.evals[k] - CCALL(conj)(ov_data.evals[i]));
      if(dist < mindist)
      {
       mindist = dist;
       ki = k;
      };
     };
     k = ki;
     memcpy(                               tmp, &(evd.evecs[i*VOL*NCOLORS*NDIRAC]), VOL*NCOLORS*NDIRAC*sizeof(t_complex));
     memcpy(&(evd.evecs[i*VOL*NCOLORS*NDIRAC]), &(evd.evecs[k*VOL*NCOLORS*NDIRAC]), VOL*NCOLORS*NDIRAC*sizeof(t_complex));
     memcpy(&(evd.evecs[k*VOL*NCOLORS*NDIRAC]),                                tmp, VOL*NCOLORS*NDIRAC*sizeof(t_complex));
     t_complex tmpval = evd.evals[i];
     evd.evals[i] = evd.evals[k];
     evd.evals[k] = tmpval;
#ifdef NOISY_OUTPUT
     fprintf(stdout, "Exchanging %i and %i ...\n", i, k);
     fflush(stdout);
#endif
    };
    free(tmp);
    fprintf(stdout, "\n DONE: sorting left evals/evecs! \n\n");
    fflush(stdout);
    ov_data.levecs = evd.evecs;
    evd.evecs = NULL;
    evd.evals = NULL;
    /* Additionally rotating revecs in span(revecs), so that they are completely orthogonal */

    fprintf(stdout, "\n Additionally rotating revecs in span(revecs), so that they are completely orthogonal ... \n");
    fflush(stdout);
    reorthogonalize(ov_data.levecs, ov_data.revecs, ov_data.n_proj);
    fprintf(stdout, "\n DONE: Additionally rotating revecs in span(revecs) ... \n");
    fflush(stdout);

#else
    if (evd.nev >= 2)
        qsort(evd.evals, evd.nev, sizeof(t_complex), &cmplx_cmp_mod);
    ov_data.l_sq_max = CCALL(creal)(evd.evals[evd.nev-1])*CCALL(creal)(evd.evals[evd.nev-1]);
    /* enlarge spectral window to avoid nonconvergance of Cheb. series */
    ov_data.l_sq_max += SPEC_MULT*ov_data.l_sq_max;
    eps = ov_data.l_sq_min/ov_data.l_sq_max;
    ov_data.coef = NULL;
    minmax_pol(eps, delta, &(ov_data.deg), &(ov_data.coef));
    eps = sqrt(eps);
    if (ov_data.coef == NULL)
        panic("%s: ov_data.coef == NULL", __func__);
    if (evd.evals != NULL) {
        free(evd.evals);
        evd.evals = NULL;
    }
    if (evd.evecs != NULL) {
        free(evd.evecs);
        evd.evecs = NULL;
    }
#endif
    /* stats */

#ifdef _OPENMP
    fprintf(stdout,"Running in parallel mode with %i threads \n", omp_get_max_threads());
#endif
    fprintf(stdout,"Output file:\t%s\n", oname);
    fprintf(stdout,"Size       :\t%ix%ix%ix%i\n", LS, LS, LS, LT);
    fprintf(stdout,"Group      :\tSU(%i)\n", NCOLORS);
    fprintf(stdout,"Tolerance  :\t%12.10f\n", tol);
    fprintf(stdout,"rho        :\t%12.10f\n", rho);
    fprintf(stdout,"kappa      :\t%12.10f\n", kappa);
    fprintf(stdout,"mass       :\t%12.10f\n", mass);

    fprintf(stdout,"anisotropy :\t%12.10f\n", anisotropy_factor);
    fprintf(stdout,"type       :\t%i\n", optype);
#ifdef GAUSSIAN_SOURCE
    fprintf(stdout,"\t\t GAUSSIAN SOURCE with radius %2.4lf\n", (double)r0);
#endif
#ifndef MU
    fprintf(stdout,"l_min      :\t%12.10f\n", sqrt(ov_data.l_sq_min));
    fprintf(stdout,"l_max      :\t%12.10f\n", sqrt(ov_data.l_sq_max));
    fprintf(stdout,"eps        :\t%12.10f\n", eps);
    fprintf(stdout,"delta      :\t%12.10f\n", delta);
#endif
    fprintf(stdout,"deg        :\t%i\n", ov_data.deg);
#ifdef MU
    fprintf(stdout, "\n\n !!!!!!!!!!!!!Running with nonzero Mu!!!!!!!!!!!!!!!!! \n\n");
    fprintf(stdout, "Mu = %4.4f\n", Mu);
    evd.which[0] = 'S';
    evd.which[1] = 'M';
    evd.data = &(ov_data);
#else
#ifdef OV_DIRAC_ORTHO_PROJ
    fprintf(stdout,"projection :\tON\n");
#else
    fprintf(stdout,"projection :\tOFF\n");
#endif
    /* */
    if (polacc && freewd) {
        fprintf(stdout, "Polynomial accelaration doesn't work with degenerate spectrum\n");
        fprintf(stdout, "Do not turn Polynomial accelaration.\n");
        exit(EXIT_FAILURE);
    }
    if (polacc) {
        fprintf(stdout, "PolAcc     :\tON\n");
        evd.which[0] = 'L';
        evd.which[1] = 'M';
        fprintf(stdout, "lambda1    :\t%f\n", lambda1);
        fprintf(stdout, "P_deg      :\t%i\n", polacc_deg);
        build_mapping_ov(lambda1, polacc_xb, &polacc_data);
        fprintf(stdout, "XB         :\t%f\n", polacc_xb);
        fprintf(stdout, "C          :\t%f\n", polacc_data.c);
        fprintf(stdout, "E          :\t%f\n", polacc_data.e);
        polacc_data.op = evd.op;
        polacc_data.deg = polacc_deg;
        polacc_data.data = &(ov_data);
        evd.op = &pol_acc_ov;
        evd.data = &polacc_data;
    }  else {
        fprintf(stdout, "PolAcc     :\tOFF\n");
        evd.which[0] = 'S';
        evd.which[1] = 'M';
        evd.data = &(ov_data);
    }
#endif
    evd.val_vec = 1;
    evd.tol = tol;
    evd.nev = snev;
    gv = lat_gauge_create();
    lat_mov_init();



    if (freewd)
    {
     lat_gauge_identity(gv);
    }
    else
     if(instanton)
     {
      ic[0] = LS/2.0+0.5;
      ic[1] = LS/2.0+0.5;
      ic[2] = LS/2.0+0.5;
      ic[3] = LT/2.0+0.5;
      lat_gauge_instanton(gv, instR, ic);
     }
     else
     {
      lat_gauge_load(gv, iname);
      action = lat_gauge_action(gv);
      fprintf(stdout,"action     :\t%12.10f\n", action);
     };

    lat_gauge_check_unit(gv);
    if(aperiodic)
     lat_gauge_aperiodic(gv);

#ifdef EFIELD
    init_efield(E,H);
    t_real FE = 2*M_PI*(double)E/(double)(LS*LT);
    t_real FH = 2*M_PI*(double)H/(double)(LS*LS);
    fprintf(stdout,"Constant electric and magnetic fields of magnitude %4.4lf and %4.4lf in lattice units are ON!!!\n",FE,FH);
    fprintf(stdout,"E is directed along X, H is directed along Z.\n");
    fflush(stdout);
#endif
#ifdef MU
    /* Restoring the square roots of eigenvalues with a correct sign!!! */
    fprintf(stdout, "\n Restoring the square roots of eigenvalues with a correct sign... \n");
    fflush(stdout);
    t_complex* artmp = (t_complex *)malloc(VOL*NDIRAC*NCOLORS*sizeof(t_complex));
    for(i=0; i<ov_data.n_proj; i++)
    {
     t_cds_vector* Ain  = (t_cds_vector *)&(ov_data.revecs[i*VOL*NDIRAC*NCOLORS]);
     t_cds_vector* Aout = (t_cds_vector *)&(artmp[0]);
     h_wdirac(&(ov_data), gv, Aout, Ain);
     ov_data.evals[i] = innprod(&(ov_data.levecs[i*VOL*NDIRAC*NCOLORS]), &(artmp[0]));
    };
    free(artmp);
    fprintf(stdout, "DONE: Restoring the square roots of eigenvalues with a correct sign.\n");
    fflush(stdout);
#endif

#ifdef DEFLATE_PROP
 fprintf(stdout, "\n\n\t Calculating %i smallest eigenvalues for deflation ... \n\n", snev);
#ifdef TIMING
 init_time();
#endif
 ov_data.mass = 0.0;
 arnoldi(gv, &evd, ncv_factor_nom, ncv_factor_denom);
#ifdef TIMING
 print_time();
#endif
 int nzero = sort_ev(evd.evecs, evd.evals, snev);
 t_complex* evecs = (t_complex* )malloc(nzero*VOL*NCOLORS*NDIRAC*sizeof(t_complex));
 ov_data.mass = mass;
 memcpy(evecs, evd.evecs, nzero*VOL*NCOLORS*NDIRAC*sizeof(t_complex));
#endif

 t_cds_vector *chck, *evec;
 if(posix_memalign((void*)&chck,        16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for chck      (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&evec,        16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for chck      (%i Kb)", __func__, sizeof(t_cds_vector)/1024);

#ifdef DEFLATE_PROP
 ov_data.mass = 0.0;
 fprintf(stdout, "\n\n\t Overlap: eigenvalues and ev.norms (%i smallest eigenvalues, %i zero modes):\n\n", snev, nzero);
 for(i=0; i<nzero; i++)
 {
  memcpy((t_complex* )evec, &(evecs[i*VOL*NCOLORS*NDIRAC]), VOL*NCOLORS*NDIRAC*sizeof(t_complex));
  ov_dirac_mm(&ov_data, gv, chck, evec);
  fprintf(stdout, "\n\n\t %i: |f_i|: %4.6E, |D f_i - l_i f_i|: %4.6E \n", i, vnorm(&(evecs[i*VOL*NCOLORS*NDIRAC])), vnorm((t_complex *)chck));
 };

 fprintf(stdout, "\n\n\t Inn. products:\n\n");
 for(i=0; i<nzero; i++)
 {
  for(j=0; j<nzero; j++)
   fprintf(stdout, "%2.2E ", CCALL(cabs)( innprod( &(evecs[i*VOL*NCOLORS*NDIRAC]), &(evecs[j*VOL*NCOLORS*NDIRAC]) ) ) );
  fprintf(stdout, "\n");
 };
 ov_data.mass = mass;
#endif

// !!!!!!!!!!!!!!!
 t_cds_vector *source, *solution, *tsol;
 if(posix_memalign((void*)&source,      16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for source (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&solution,    16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for solution      (%i Kb)", __func__, sizeof(t_cds_vector)/1024);
 if(posix_memalign((void*)&tsol,        16, sizeof(t_cds_vector)))
  panic("%s: can't allocate memory for tsol      (%i Kb)", __func__, sizeof(t_cds_vector)/1024);

 if(posix_memalign((void *)&evd.evecs,        16, NCOLORS*NDIRAC*sizeof(t_cds_vector)))
   panic("%s: can't allocate memory for evd.evecs        (%i Kb)", __func__, NCOLORS*NDIRAC*sizeof(t_cds_vector)/1024);

 for(a=0; a<NCOLORS; a++)
  for(id=0; id<NDIRAC; id++)
  {
   i = NDIRAC*a + id;
   fprintf(stdout, "\n\n\t Inverting the Overlap Dirac operator for source %i ... \n\n", i+1);
   set_zero((t_complex *)source);
   // preparing a gaussian source
   x0 = 0;
#ifdef GAUSSIAN_SOURCE
   for(x = 0; x<VOL; x++)
   {
    dist = distance4D(x, x0);
    *((t_complex *)source + NDIRAC*NCOLORS*x + i) = exp(-0.5*dist*dist/(r0*r0));
   };
#else
   *((t_complex *)source + NDIRAC*NCOLORS*x0 + i) = 1.0 + I*0.0;
#endif
#ifdef DEFLATE_PROP
   // projecting out zero modes
   for(i=0; i<nzero; i++)
   {
    zb = innprod(&(evecs[i*VOL*NCOLORS*NDIRAC]), (t_complex *)source);
    xpcby((t_complex *)source, &(evecs[i*VOL*NCOLORS*NDIRAC]), -1.0*zb);
   };
#endif
   dist = vnorm((t_complex *)source);
   ax((t_complex *)source, 1.0/dist);

#ifdef TIMING
   init_time();
#endif

   ret = shumr(&ov_data, gv, source, solution, tol, IMAX);
   if(ret!=0)
   {
    fprintf(stdout, "\n\n\t SHUMR has not reached the precision %4.4E in allowed number %i of iterations at source %i. Program terminates.\n\n", tol, IMAX, NDIRAC*a + id);
    fprintf(stderr, "\n\n\t SHUMR has not reached the precision %4.4E in allowed number %i of iterations at source %i. Program terminates.\n\n", tol, IMAX, NDIRAC*a + id);
    return EXIT_FAILURE;
   };

   ov_dirac_mm(&ov_data, gv, chck, solution);
   xpcby((t_complex* )chck, (t_complex* )source, -1.0);
   fprintf(stdout, "\n\n\t Reconstructed source differs from the original by %6.6E \n\n", (double)vnorm((t_complex*)chck));
   ov_data.mass = 0.0;
   ov_dirac_mm(&ov_data, gv, tsol, solution);
   ov_data.mass = mass;
   xpby((t_complex *)solution, (t_complex *)tsol, -0.5);
#ifdef TIMING
   print_time();
#endif
   memcpy(&(evd.evecs[i*VOL*NCOLORS*NDIRAC]), (t_complex *)solution, sizeof(t_cds_vector));
  };

  evd.nev   = NCOLORS*NDIRAC;
  evd.evals = (t_complex *)malloc(NCOLORS*NDIRAC*sizeof(t_complex));
  for(i=0; i<NCOLORS*NDIRAC; i++)
   evd.evals[i] = (double)i;
  write_ev(oname, "shadow", &evd, SM_SECTION);

  if(evd.evals != NULL)
  {
   free(evd.evals);
   evd.evals = NULL;
  };
  if(evd.evecs != NULL)
  {
   free(evd.evecs);
   evd.evecs = NULL;
  };

  lat_gauge_destroy(gv);
  free(source);
  free(solution);
  free(tsol);
  free(chck);
#ifdef DEFLATE_PROP
  free(evecs);
#endif
  free(evec);

  fprintf(stdout,"\n\n OVERLAP run succesfully completed!!! \n\n");
  return EXIT_SUCCESS;
}/*}}}*/
