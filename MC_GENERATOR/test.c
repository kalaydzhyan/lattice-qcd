#include <stdio.h>
#include <unistd.h>
#include<complex.h>

#include"SU_N-utils.h"
#include"MC-SU_N.h"

#define LS 6
////#define SYSTEM_RANDOM


float get_plaq(int x, int mu, int nu)
{
 SU2 tmp1, tmp2;
 tmp1 = SU2_mult(F[link_index(mu,x)], F[link_index(nu,index_up(x,mu))] );
 tmp2 = SU2_mult(F[link_index(nu,x)], F[link_index(mu,index_up(x,nu))] );
 tmp1 = SU2_mult(tmp1,SU2_conj(tmp2));
 return crealf(tmp1.alpha);
}
double get_plaq_SU_N(int x, int mu, int nu)
{
 SU_N tmp1, tmp2;
 tmp1 = SU_N_mult(F_N[link_index(mu,x)], F_N[link_index(nu,index_up(x,mu))] );
 tmp2 = SU_N_mult(F_N[link_index(nu,x)], F_N[link_index(mu,index_up(x,nu))] );
 tmp1 = SU_N_mult(tmp1,SU_N_conj(tmp2));
 return SU_N_norm(tmp1);
}

void print_latt_SU_N()
{
 for(int x=0; x<param->ipw[param->D]; x++)
  for(int mu=0; mu<param->D; mu++){
	for( int i = 0; i < SU_N_RANK; i++ ){
		printf("\n");
     		for(int  j = 0; j < SU_N_RANK; j++ )
  			 printf(" %f \t",fabs(creal(F_N[link_index(mu,x)].U[i][j])));  
	}; 
	printf(" \n");
  }
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
  
}
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
return mean/(double)dim1;
}
int main()
{
FILE* f;
f=fopen("1.dat","w");

	int i=0,j=0;
	    SU_N_init(NULL, LS, LS, 4);
            
	  cold_start_SU_N( LS, LS, 4);
	rnd_seed_set();
param->beta = 2.;

//param->uzero = .1*param->beta;
  //param_improvement();
for(j=0;j<20;j++){
	    param->beta += 1.;
#ifdef IMPROVEMENT
	    param->uzero = .791627;
    param_improvement();
#endif
for(i=0;i<20;i++){
	MC_SU_N();
	fprintf(stderr,"%f\n",SU_N_mean_plaq());
};
fprintf(f,"%f\n",SU_N_mean_plaq());
};
/*
		sprintf(filename,"DATA/%ix%ib%fu%f.0.dat",LS,LS,param->beta,param->uzero);
	
	for(i=0;i<100;i++){
		sprintf(filename,"DATA/%ix%ib%fu%f.%i.dat",LS,LS,param->beta,param->uzero,i+1);
		SU_N_init(filename, LS, LS, 4);
		fprintf(stderr,"%s\t%f\n",filename,SU_N_mean_plaq());

		};*/

	SU_N_free();
fclose(f);
    return 0;
}
