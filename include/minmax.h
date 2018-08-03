#ifndef _MINMAX_H_
#define _MINMAX_H_

double eval_pol_cheb(double eps, int n, double* c, double x, double* C);
double eval_pol_cheb_cc(double eps, int n, double* c, double x);
double minmax_pol_deg(double eps, int deg, double* c);
void minmax_pol(double eps, double delta, int* deg, double** c);

#endif
