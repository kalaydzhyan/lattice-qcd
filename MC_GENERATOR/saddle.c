#include<stdlib.h>
#include<stdio.h>
#include<unistd.h>
#include<stdbool.h>
#include<time.h>
#include<complex.h>

#include"SU2-utils.h"
#include"MC-SU2.h"

#define LS 16

float DMCP_func()
{
 float q=0;     
 for(int x=0; x<param->ipw[param->D]; x++)
  for(int mu=0; mu<param->D; mu++)
   q+=fabs(creal(F[link_index(mu,x)].alpha));
 return q/(param->D*param->ipw[param->D]);
}

float get_plaq(int x, int mu, int nu)
{
 SU2 tmp1, tmp2;
 tmp1 = SU2_mult(F[link_index(mu,x)], F[link_index(nu,index_up(x,mu))] );
 tmp2 = SU2_mult(F[link_index(nu,x)], F[link_index(mu,index_up(x,nu))] );
 tmp1 = SU2_mult(tmp1,SU2_conj(tmp2));
 return crealf(tmp1.alpha);
}

void empty_latt()
{
 for(int x=0; x<param->ipw[param->D]; x++)
  for(int mu=0; mu<param->D; mu++)
  {    
   F[link_index(mu,x)].alpha = 1;
   F[link_index(mu,x)].beta  = 0;
  }; 
}

void print_latt()
{
 for(int x=0; x<param->ipw[param->D]; x++)
  for(int mu=0; mu<param->D; mu++)
   printf(" %f \n",fabs(creal(F[link_index(mu,x)].alpha)));    
}

int main()
{
    srand ( time(NULL) );
    
    fields_init(NULL, LS, LS, 4);
    param->beta = 1.;
    print_latt();
    MC_SU2();
   // print_latt();
    MaA_center_gauge_SU2();
   // print_latt();
        
    fields_free();
    return 0;
}

