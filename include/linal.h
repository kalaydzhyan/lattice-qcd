#ifndef _LINAL_H_
#define _LINAL_H_

#include "types.h"

void ag5xpby(t_complex* x, const t_real a, const t_complex* y, const t_real b);
void ag5xpcby(t_complex* x, const t_complex a, const t_complex* y, const t_complex b);
void xbyg5(t_complex* x);
void xpby(t_complex* x, const t_complex* y, const t_real b);
void xpcby(t_complex* x, const t_complex* y, const t_complex b);
void xcpy(t_complex* x, const t_complex* y);
void axpby(t_complex *x, const t_real a, const t_complex* y, const t_real b);
void axpbyn(t_complex *x, const t_real a, const t_complex* y, const t_real b, int n);
void axpbynd(t_double_complex *x, const t_double_real a, const t_double_complex* y, const t_double_real b, int n);
void  ax(t_complex *x, const t_real a);
void axc(t_complex *x, const t_complex a);
t_complex innprod(t_complex* x, t_complex* y);
t_real    vnorm(t_complex* x);
t_real    adiffn(t_complex* x, t_complex* y, int n); // ||x - y||
t_double_real    adiffnd(t_double_complex* x, t_double_complex* y, int n); // ||x - y||
void set_zero(t_complex* x); // sets all elements of x to zero
void g5x(t_complex* x, t_complex* y);

void inverse(t_double_complex *A, t_double_complex *R, int size);
void reorthogonalize(t_complex *L, t_complex *R, int n_vecs); // Reorthogonalizes L and R, so that L_i*R_k = \delta_ik, i,k = 0..n_vecs-1
void reorthogonalize_gs(t_complex *vecs, int n_vecs); //reorthogonalizes vecs

void check_ortho(int nev, int vol, t_complex* evals, t_complex* evecs);
void init_CDSV(t_cds_vector*, int*, int, int);
#endif
