#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

typedef double t_real;
typedef double complex t_complex;

#define CCALL(callname) (callname)

/* tests the SHUMR method, A. Borici, A. Allkoci, hep-lat/0602015 */

#define NC 100

t_complex V[NC][NC];
t_real    c1 = 0.3; // D = c1*I + c2*V, V is unitary
t_real    c2 = 0.5;

void init_V()
{
 int i, j;
 double rev, imv;
 FILE *vfile = fopen("/home/pools/1/buividovich/overlap/rand_unit100x100.dat", "r");
 for(i=0; i<NC; i++)
  for(j=0; j<NC; j++)
  {
   fscanf(vfile, "%lf %lf", &rev, &imv);
   V[i][j] = (t_real)rev + I*(t_real)imv;
  }; 
 fclose(vfile);  
}

void print_V()
{
 int i, j;
 for(i=0; i<NC; i++)
  for(j=0; j<NC; j++)
   printf("%i %i: %2.4lf + I*%2.4lf\n", i, j, CCALL(creal)(V[i][j]), CCALL(cimag)(V[i][j])); 
}

void print_carr(t_complex* x, int n)
{
 int i;
 for(i=0; i<n; i++)
  printf("%i: %6.6lf + I*%6.6lf\n", i, CCALL(creal)(x[i]), CCALL(cimag)(x[i]));    
}

void check_unit()
{
 t_complex  sp;
 t_real    nsp;
 int i, j, k;
 for(i=0; i<NC; i++)
  for(j=0; j<NC; j++)
  {
   sp = 0.0 + I*0.0;
   for(k=0; k<NC; k++)
    sp += V[i][k]*CCALL(conj)(V[j][k]);
   sp  = (i==j)? sp - 1.0 : sp;
   nsp = CCALL(cabs)(sp);
   if(nsp>1E-7)
    printf("Matrix not unitary!!! nsp = %2.4E\n", nsp);       
  }; 
}

void D_op(t_complex *in, t_complex *out)
{
 int i, j;
 for(i=0; i<NC; i++)
 {
  out[i] = c1*in[i];
  for(j=0; j<NC; j++)
   out[i] += c2*V[i][j]*in[j];
 }; 
}

void set_zero(t_complex *x)
{
 int i;
 for(i=0; i<NC; i++)
  x[i] = 0.0 + I*0.0;
}  
  
t_real vnorm(t_complex *x)
{
 int i;
 t_real q = 0.0;
 for(i=0; i<NC; i++)
  q += CCALL(creal)(x[i]*CCALL(conj)(x[i]));
 return sqrt(q); 
}

void  axc(t_complex *x, const t_complex a)
{
 int i;
 for(i=0; i<NC; i++)
  x[i] *= a;      
}

t_complex innprod(t_complex* x, t_complex* y)
{
	int i;
	double complex sum = 0.0 + I*0.0;
	for (i = 0; i<NC; i++)
     sum += CCALL(conj)(x[i])*y[i];

	return (t_complex) sum;
}

void xpcby(t_complex* x, const t_complex* y, const t_complex b)
{
 int i;
 for(i = 0; i < NC; i++)
  x[i] += b*y[i];
}

void shumr(t_complex *source, t_complex *x, t_real tol, int imax) //x is the solution
{
 t_complex c_k, s_k, a, t11, mu, nu, mu_k, nu_k, omega, theta, gamma, L11, iqv, eps;
 t_real L21, rnorm_p;
 
 t_complex *q_tilde     = (t_complex* )malloc(NC*sizeof(t_complex));
 t_complex *Dp          = (t_complex* )malloc(NC*sizeof(t_complex));
 t_complex *v    = (t_complex* )malloc(NC*sizeof(t_complex));
 t_complex *w    = (t_complex* )malloc(NC*sizeof(t_complex));
 t_complex *s    = (t_complex* )malloc(NC*sizeof(t_complex)); 
 t_complex *p    = (t_complex* )malloc(NC*sizeof(t_complex));
 t_complex *rp     = (t_complex* )malloc(NC*sizeof(t_complex));
 t_complex *xp     = (t_complex* )malloc(NC*sizeof(t_complex));
 
 set_zero(x);
 
 t_complex *r     = (t_complex* )malloc(NC*sizeof(t_complex));
 memcpy(r, source, NC*sizeof(t_complex));
 
 t_real    rho   = vnorm(r);
 t_real    rnorm = rho;
 t_complex alpha = rho; 
 
 t_complex *q     = (t_complex* )malloc(NC*sizeof(t_complex));
 memcpy(q, source, NC*sizeof(t_complex));
 axc(q, 1.0/rho);
 
 t_complex u12       = 0.0 + I*0.0;
 t_complex beta      = 1.0 + I*0.0; 
 t_complex L11_tilde = 1.0 + I*0.0; 
 
 t_complex *q_old     = (t_complex* )malloc(NC*sizeof(t_complex));
 set_zero(q_old); 
 t_complex *v_old     = (t_complex* )malloc(NC*sizeof(t_complex));
 set_zero(v_old); 
 t_complex *w_old     = (t_complex* )malloc(NC*sizeof(t_complex));
 set_zero(w_old);
 t_complex *s_old     = (t_complex* )malloc(NC*sizeof(t_complex));
 set_zero(s_old);
 
 t_complex c_km1 = 1.0 + I*0.0;
 t_complex s_km1 = 0.0 + I*0.0;
 t_complex c_km2 = 0.0 + I*0.0;
 t_complex s_km2 = 0.0 + I*0.0;
 
 t_complex *p1    = (t_complex* )malloc(NC*sizeof(t_complex));
 set_zero(p1); 
 t_complex *p2    = (t_complex* )malloc(NC*sizeof(t_complex));
 set_zero(p2);
 t_complex *Dp1    = (t_complex* )malloc(NC*sizeof(t_complex));
 set_zero(Dp1); 
 t_complex *Dp2    = (t_complex* )malloc(NC*sizeof(t_complex));
 set_zero(Dp2);
 
 int counter = 0;
 
 while((rnorm > tol) && (counter<=imax))
 {
  D_op(q, v);
  if(counter > 0)
  {
   iqv = innprod(q_old, v_old);
   //if(CCALL(cabs)(iqv)>1E-10) //TODO:  compare with tol?
   u12 = - 1.0*innprod(q_old, v)/iqv;
   /*else
   {
    printf("iqv = 0!!!\n");
    return;   
   };*/      
  }; 
  gamma = -1.0*c1*u12;
  L11 = innprod(q, v) + u12*innprod(q, v_old);
  memcpy(q_tilde, v, NC*sizeof(t_complex));
  xpcby(q_tilde,     q, -1.0*L11);
  xpcby(q_tilde, v_old, u12);
//L21 is real
  L21 = vnorm(q_tilde);
  if(L21 <= tol)
  {
   printf("Broken at %i: L21 = %4.6E <= tol\n", counter, L21);
   break;
  };
  
  memcpy(w, q, NC*sizeof(t_complex));
  xpcby(q, q_old, u12);
  xpcby(q, w_old, gamma/L11_tilde);

//s=c1*(q+q_old*u12)+c2*(v+v_old*u12)+s_old*gamma/L11_tilde;
  memcpy(s, q, NC*sizeof(t_complex));
  axc(s, c1);
  xpcby(s, q_old, c1*u12);
  xpcby(s,     v,     c2);
  xpcby(s, v_old, c2*u12);
  xpcby(s, s_old, gamma/L11_tilde);
  
  L11_tilde = c1 + c2*L11 - beta*gamma/L11_tilde;
  alpha     = alpha*beta/L11_tilde;

//x=x+w*alpha;
//r=r-s*alpha;
  xpcby(x, w,      alpha);
  xpcby(r, s, -1.0*alpha);
//q_old=q; v_old=v; w_old=w; s_old=s;
  memcpy(q_old, q, NC*sizeof(t_complex));
  memcpy(v_old, v, NC*sizeof(t_complex));
  memcpy(w_old, w, NC*sizeof(t_complex));
  memcpy(s_old, s, NC*sizeof(t_complex));
  
//q=q_tilde/L21;
  memcpy(q, q_tilde, NC*sizeof(t_complex)); 
  axc(q_tilde, 1.0/L21);
  
  beta = -1.0*c2*L21;
  rnorm = vnorm(r);

  t11 = c1 + c2*L11;
  mu  = t11*c_km1 + gamma*CCALL(conj)(s_km1)*c_km2;
  nu  = c2*L21;
  if(CCALL(cabs)(mu) > 1E-10) //TODO: tolerance??? 
  {
    c_k = CCALL(cabs)(mu)/sqrt(CCALL(cabs)(mu)*CCALL(cabs)(mu) + CCALL(cabs)(nu)*CCALL(cabs)(nu));
    s_k = CCALL(conj)(c_k*nu/mu);
  }
  else
  {
    printf("Exceptional mu = %6.6E encountered!!!\n", mu);  
    c_k = 0.0 + I*0.0;
    s_k = 1.0 + I*0.0;
  };
  omega = nu*alpha*s_k;
  mu_k  = c_k*mu + s_k*nu;
  eps   = t11*s_km1 - gamma*c_km1*c_km2;
  theta = -1.0*gamma*s_km2;
  //p = q/mu_k + q_old*u12/mu_k - p1*eps/mu_k - p2*theta/mu_k;
  memcpy(p, q, NC*sizeof(t_complex));
  axc(p, 1.0/mu_k);
  xpcby(p, q_old, u12/mu_k);
  xpcby(p, p1, -1.0*eps/mu_k);
  xpcby(p, p2, -1.0*theta/mu_k);
  //Dp = c1/mu_k*q + c1*u12/mu_k*q_old + c2/mu_k*v + c2*u12/mu_k*v_old - Dp1/mu_k*eps - Dp2*theta/mu_k;
  memcpy(Dp, q, NC*sizeof(t_complex));
  axc(Dp, c1/mu_k);
  xpcby(Dp, q_old,     c1*u12/mu_k);
  xpcby(Dp,     v,         c2/mu_k);
  xpcby(Dp, v_old,     c2*u12/mu_k);
  xpcby(Dp,   Dp1,   -1.0*eps/mu_k);
  xpcby(Dp,   Dp2, -1.0*theta/mu_k);

  //rnorm_p = vnorm(r + omega*Dp);
  memcpy(rp, r, NC*sizeof(t_complex));
  xpcby(rp, Dp, omega);
  rnorm_p = vnorm(rp);
  
  //xp=x-omega*p;
  memcpy(xp, x, NC*sizeof(t_complex));
  xpcby(xp, p, -1.0*omega);

  c_km2 = c_km1;
  s_km2 = s_km1;
  c_km1 = c_k;
  s_km1 = s_k;
  memcpy( p2,  p1, NC*sizeof(t_complex));
  memcpy(Dp2, Dp1, NC*sizeof(t_complex));
  memcpy( p1,   p, NC*sizeof(t_complex));
  memcpy(Dp1,  Dp, NC*sizeof(t_complex));
  
  counter++;
  printf("Hello from SHUMR's %i's iteration, rnorm = %6.6E...\n", counter, rnorm);
 };
}

int main ()
{
 printf("Hello world!!!\n");   
 init_V();
// printf("The matrix V:\n");
// print_V();
 check_unit();
 
 t_complex *source     = (t_complex* )malloc(NC*sizeof(t_complex));
 set_zero(source);
 source[5] = 1.0 + I*0.0;
 t_complex *solution     = (t_complex* )malloc(NC*sizeof(t_complex));
 set_zero(solution);
 
 shumr(source, solution, 0.0000001, 200);
 
 printf("Solution:\n");
// print_carr(solution, NC);
 
 t_complex *chck     = (t_complex* )malloc(NC*sizeof(t_complex));
 D_op(solution, chck);
 
 printf("Reconstructed source:\n\n");
// print_carr(chck, NC);
 
 xpcby(chck, source, -1.0);
 printf("The reconstructed source differs from the original by %6.6lf\n", vnorm(chck));
 
/* solution[0] =  0.777866 - I*0.117242;
 solution[1] = -1.1654   - I*0.104288;
 solution[2] =  0.467216 - I*0.435803;
 solution[3] = -0.702001 + I*0.094379; 
 solution[4] =  0.592778 + I*0.28956;
 
 D_op(solution, chck);
 
 printf("Reconstructed source (exact solution):\n\n");
 print_carr(chck, NC); */
 
 system("PAUSE");
 return EXIT_SUCCESS;
}
