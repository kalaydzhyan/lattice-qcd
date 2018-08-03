#include<stdlib.h>
#include<stdio.h>
#include<unistd.h>
#include<stdbool.h>
#include<time.h>
#include<complex.h>

#include"SU2-utils.h"
#include"MC-SU2.h"
#include"hist/hist.h"

#define N_RUNS 200
#define MEASUREMENT_INTERVAL 50

int get_plaq_sign(int x, int mu, int nu)
{
 SU2 tmp1, tmp2;
 tmp1 = SU2_mult(F[link_index(mu,x)], F[link_index(nu,index_up(x,mu))] );
 tmp2 = SU2_mult(F[link_index(nu,x)], F[link_index(mu,index_up(x,nu))] );
 tmp1 = SU2_mult(tmp1,SU2_conj(tmp2));
 int q = 1;
 if(crealf(tmp1.alpha)<0)
  q=-1;
 return q;
}

void cube_dual_to_link(int mu, int *a, int *b, int *c)
{
 if(mu==0) {*a=1; *b=2; *c=3;};
 if(mu==1) {*a=0; *b=2; *c=3;};
 if(mu==2) {*a=0; *b=1; *c=3;};
 if(mu==3) {*a=0; *b=1; *c=2;};
}

int z2monopole(int x, int mu) //Finds whether a Z2 monopole sits on (x,mu) link of the dual lattice. returns 1 for links with monopoles
{
 int q=1;
 int a,b,c;
 cube_dual_to_link(mu, &a, &b, &c);
 q = q*get_plaq_sign(x,a,b)*get_plaq_sign(index_up(x,c),a,b);
 q = q*get_plaq_sign(x,a,c)*get_plaq_sign(index_up(x,b),a,c);
 q = q*get_plaq_sign(x,b,c)*get_plaq_sign(index_up(x,a),b,c);
 q = (1 - q)/2;
 return q;
}

float z2m_dens()
{
 int nz2m = 0;
 for(int x=0; x<param->ipw[param->D]; x++)
  for(int mu=0; mu<param->D; mu++)
   nz2m+=z2monopole(x,mu);
 return (float)nz2m/(float)(param->ipw[param->D]*param->D);
}

int main()
{
 fields_init(NULL, 12, 12, 4);

 FILE *pfile = fopen("z2m_dens.dat","w");

 for(int nb=0; nb<8; nb++)
 {
  float az2m = 0;
  int nz2m = 0;
  param->beta = 2.25 + 0.05*nb;
  for(int r=0; r<N_RUNS; r++)
  {
   double mp = MC_SU2();
   printf("beta = %f, %i runs of %i, mean plaquette: %f \n",param->beta,r,N_RUNS,mp);
   if(r>MEASUREMENT_INTERVAL)
   {
    az2m+=z2m_dens();
    nz2m++;
   };
  };
  az2m=az2m/(float)nz2m;
  float ls = lattice_spacing(param->beta);
  fprintf(pfile,"%4.6f %4.6f %4.6f %4.6f %4.6f\n",param->beta, ls, az2m, az2m*param->D/(ls*ls*ls), az2m*param->D/ls);
 };

 fields_free();
 fclose(pfile);
 
 return 0;
};
 
