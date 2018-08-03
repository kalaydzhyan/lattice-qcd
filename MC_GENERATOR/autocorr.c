#include <stdio.h>
#include <unistd.h>
#include<complex.h>
#include<stdlib.h>

/*#include"SU_N-utils.h"
#include"MC-SU_N.h"

#define LS 8*/
////#define SYSTEM_RANDOM
/*float DMCP_func()
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
void print_latt_SU_N()
{
 for(int x=0; x<param->ipw[param->D]; x++)
  for(int mu=0; mu<param->D; mu++){
	for( int i = 0; i < SU_N_RANK; i++ ){
	printf("\n");
      for(int  j = 0; j < SU_N_RANK; j++ )
   printf(" %f \t",fabs(creal(F_N[link_index(mu,x)].U[i][j])));  
};  printf(" \n");
};
}
void cold_start_SU_N(uchar T, uchar L, uchar D )
{
  if( D != 3 && D != 4 ) panic("I'm too lazy, only D=3,4 is implemented");
  if( !param ){
    if( param_init(T,L,D) ) panic("Parameters initialization failed");
  }else{
    if( param->D != D || param->size[0] != T || param->size[1] != L )
      panic("Floating geomety?!");
  };
  if( !F_N )
    MEMORY_CHECK( F_N = (SU_N *)calloc( param->D*param->ipw[param->D], sizeof(SU_N)), -1);
  for(int x=0; x<param->ipw[param->D]; x++)
  for(int mu=0; mu<param->D; mu++)
    for( int i = 0; i < SU_N_RANK; i++ )
    {
	for(int  j = 0; j < SU_N_RANK; j++ )
	  F_N[link_index(mu,x)].U[i][j]=0;
        F_N[link_index(mu,x)].U[i][i]=1;
    };
  
}*/
double  autocorr(double *a,double *m,int dim)
{
int i,l;
	double mean=0.,sqmean=0.,mean1=0.,sqmean1=0.;
int dim1=dim-dim/4;
	for(i=0;i<dim1;i++)
	{
		mean=mean + m[i];
		sqmean=sqmean + m[i]*m[i];
	};
	if(mean==sqmean) fprintf(stderr,"division by zero/n");
	mean1=mean;
	sqmean1=sqmean;
	for( l=0;l<dim/4;l++)
	{
		a[l]=0.;
		for(i=0;i<dim1;i++)
			a[l] =a[l] + m[i]*m[i+l];
		a[l]=(a[l]-mean*mean1/(double)(dim1))/sqrt(sqmean-mean*mean/(double)(dim1))/sqrt(sqmean1-mean1*mean1/(double)(dim1));
		mean1-=m[l];
		mean1+=m[dim1+l];
		sqmean1-=m[l]*m[l];
		sqmean1+=m[dim1+l]*m[dim1+l];
	};
return mean;
}
int main()
{

double* m;
double mean;
double* a;
a=malloc(25*sizeof(double));
m=malloc(100*sizeof(double));

int i=0;
 
  //  SU_N_init(NULL, LS, LS, 4);
   //cold_start_SU_N( LS, LS, 4);
 //   param->beta = 4.08;
 //   print_latt_SU_N();
printf("---------------------");


for(i=0;i<100;i++)
m[i]=i;
mean=autocorr(a,m,100);
for( i=0;i<25;i++){

    fprintf(stderr,"%f\n",a[i]);

};

free(a);
free(m);
    return 0;
}
