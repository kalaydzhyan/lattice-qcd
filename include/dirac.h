#ifndef _DIRAC_H_
#define _DIRAC_H_

#include "types.h"

#define SPEC_MULT			(1e-2)

#define OV_DIRAC_CONV_FNAME	("ov_dirac_conv.log")	//!< convergence history

#define	SM_SECTION			(0)
#define	LM_SECTION			(1)
#define L_SECTION           (2) // for storing left and right eigenvectors
#define R_SECTION           (3)


#define E_READ_NOERROR		(0x0)	//!< no errors during reading
#define E_READ_UNSUPFORM	(0x1)
#define E_READ_CANTREAD		(0x2)	//!< critical error: can't read file
#define E_READ_CANTOPEN		(0x4)	//!< critical error: can't open file
#define E_READ_SECNOTHERE	(0x8)	//!< critical error:
#define E_READ_BADDIM		(0x10)
#define E_READ_BADLS		(0x20)
#define E_READ_BADLT		(0x40)
#define E_READ_BADNCOLORS	(0x80)
#define E_READ_BADNDIRAC	(0x100)

t_real kappa2rho(const t_real kappa);
t_real rho2kappa(const t_real kappa);

#ifdef MU
 t_real Mu;
#endif
 t_real anisotropy_factor;

void wdirac(void* data, t_gauge_vector * gv, t_cds_vector * out, t_cds_vector * in);

void h_wdirac(void* data, t_gauge_vector * gv, t_cds_vector * out, t_cds_vector * in);
void h_wdirac_sq(void* data, t_gauge_vector * gv, t_cds_vector * out, t_cds_vector * in);

void ov_dirac_mm(void * data, t_gauge_vector* gv, t_cds_vector* out, t_cds_vector* in); 
void ov_dirac_mmc(void * data, t_gauge_vector* gv, t_cds_vector* out, t_cds_vector* in);
void ov_dirac_mm_cc(void * data, t_gauge_vector* gv, t_cds_vector* out, t_cds_vector* in); 

void arnoldi(t_gauge_vector * gv, t_op_ev* evd, int ncv_factor_nom, int ncv_factor_denom);
int shumr(void *data, t_gauge_vector *gv, t_cds_vector *source, t_cds_vector *solution, t_real tol, int imax);

int read_ev(char * fname, t_op_ev * evd, int sm_lm);
void write_ev(char * fname, char * cname, t_op_ev * evd, int sm_lm);
int cmplx_cmp_mod(const void * a, const void * b);

void proj_ortho_comp(t_complex* in, t_complex* evecs, int n_proj, t_complex* dp);


void scalar(void* data, t_gauge_vector * gv, t_cds_vector * out, t_cds_vector * in);

#endif
