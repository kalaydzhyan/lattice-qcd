#include<stdlib.h>
#include<stdio.h>
#include<unistd.h>
#include<stdbool.h>
#include<time.h>
#include<complex.h>

#include"SU2-utils.h"
#include"MC-SU2.h"
#include"hist/hist.h"

#define N_RUNS 100000
#define MEASUREMENT_INTERVAL 25

histogram *hist;

double polyakov_line()
{
 uchar ix[4] =  {0, 0, 0, 0};
 SU2 pl;
 double tpl = 0.0;
 int ntpl = 0;
 for(uint m=0; m<param->ipw[param->D]; m++)
 {
  if(m%param->size[0]==0)
  {
   if(m!=0)
   {
    tpl+=creal(pl.alpha);
    add_point(hist,creal(pl.alpha),1);
    ntpl++;
   };
   pl.alpha = 1.0;
   pl.beta  = 0.0; 
  };
  pl = SU2_mult(pl,F[link_index(0,m)]);
 };
 return tpl/ntpl; 
}

int main()
{
 srand ( time(NULL) );
    
 fields_init(NULL, 8, 16, 4);
 param->beta = 2.65;
 
 double avpl = 0.0;
 int npl=0;
 
 double avplaq = 0.0;
 int nplaq = 0; 
 
 hist = init_hist(-1, 1, 20);
  
 for(int r=0; r<N_RUNS; r++)
 {
  double mp = MC_SU2();
  printf("mean plaquette: %f \n",mp);
  printf("Polyakov line: %f \n",polyakov_line());
  if(r>50)
  {
   avplaq+=mp;
   nplaq++;
  };
  if(r>50 && (r+1)%MEASUREMENT_INTERVAL==0)
  {
   avpl+=polyakov_line();
   npl++;
   printf("\n ************************* \n \n");
   printf(" Average polyakov line after %i runs: %5.3f\n",r,avpl/npl);
   printf(" Average plaquette after %i runs: %5.3f\n",r,avplaq/nplaq);
   printf("\n ************************* \n \n");
   plot_hist_cons(hist);
  };
 }; 
  
 free_hist(hist); 
 fields_free();
 return 0;
};
