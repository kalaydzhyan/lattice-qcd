#include <stdio.h>
#include <unistd.h>
#include<complex.h>

#include"SU_N-utils.h"
#include"MC-SU_N.h"
//#undef IMPROVEMENT
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
};  printf(" \n");
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
#ifdef IMPROVEMENT
FILE* pot=fopen("pot6withoutor.dat","w");
#else
FILE* pot=fopen("potwil.dat","w");
#endif
	char filename[36];
	int x=0, i=0,j=0,k=0,l=0,m0=0,m1=0,m2=0,N_RUNS=50;
	uchar  nu=0,mu=0,sigma=0;
	double W[LS-2][LS-2],W_diag2[LS-2],W_diag3[LS-2],W_disp[LS-2][LS-2],W_disp_diag2[LS-2],W_disp_diag3[LS-2];
	for(i=0;i<LS-2;i++){
		for(j=0;j<LS-2;j++){
			W[i][j]=0.;
			W_disp[i][j]=0.;
		};
		W_diag2[i]=0.;
		W_diag3[i]=0.;	
		W_disp_diag2[i]=0.;
		W_disp_diag3[i]=0.;
		};
	SU_N tmp1, tmp2,tmp3,tmp4;

            
	cold_start_SU_N( LS, LS, 4);
	rnd_seed_set();
	param->beta = 6.8;
#ifdef IMPROVEMENT
	param->uzero = .4637;
	param_improvement();
#endif
	for(l=0;l<N_RUNS;l++){
#ifdef IMPROVEMENT
		sprintf(filename,"DATA/%ix%ib%fu%f.%i.dat",LS,LS,param->beta,param->uzero,l+1);
#else
		sprintf(filename,"DATA/WIL/%ix%ib%f.%i.dat",LS,LS,param->beta,l+1);
#endif
		SU_N_init(filename, LS, LS, 4);
		  for( x = 0; x < param->ipw[ param->D ]; x++)
			for(nu=1;nu<param->D;nu++)
				for(i=0;i<LS-2;i++){
					tmp1 = F_N[link_index(0,x)];
					m0=index_up(x,0);
					for(k=0;k<i;k++){
						tmp1=SU_N_mult(tmp1, F_N[link_index(0,m0)]);
						m0=index_up(m0,0);
					};
					for(j=0;j<LS-2;j++){
						m2=m0;
						tmp3 = F_N[link_index(nu,x)];
						tmp2 = F_N[link_index(nu,m2)];
						for(k=0;k<j;k++){
							m2=index_up(m2,nu);
							tmp2=SU_N_mult(tmp2, F_N[link_index(nu,m2)]);
						};
						m1=index_up(x,nu);
						for(k=0;k<j;k++){
							tmp3=SU_N_mult(tmp3, F_N[link_index(nu,m1)]);
							m1=index_up(m1,nu);
						};
						tmp4 = F_N[link_index(0,m1)];
						for(k=0;k<i;k++){
							m1=index_up(m1,0);
							tmp4=SU_N_mult(tmp4, F_N[link_index(0,m1)]);
						};
						if(index_up(m2,nu)!=index_up(m1,0)) fprintf(stderr,"!!!!!!!!!!!!!!!");
						tmp2 = SU_N_mult(tmp1, tmp2);
						tmp3 = SU_N_mult(tmp3, tmp4 );
						tmp2 = SU_N_mult(tmp2,SU_N_conj(tmp3));	
						W[i][j] += SU_N_norm(tmp2);
						
					};
/*------------------------diag2---------------------------------*/
					for(mu=1;mu<param->D;mu++){
						if(mu==nu) continue;
						m1=index_up(x,mu);
						tmp2=SU_N_mult(F_N[link_index(nu,m0)],F_N[link_index(mu,index_up(m0,nu))]);
						tmp3=SU_N_mult(F_N[link_index(mu,x)],F_N[link_index(nu,m1)]);
						m1=index_up(m1,nu);
						tmp4 = F_N[link_index(0,m1)];
						for(k=0;k<i;k++){
							m1=index_up(m1,0);
							tmp4=SU_N_mult(tmp4, F_N[link_index(0,m1)]);
						};
						if(index_up(index_up(m0,nu),mu)!=index_up(m1,0)) fprintf(stderr,"!!!!!!!!!!!!!!!");
						tmp2 = SU_N_mult(tmp1, tmp2);
						tmp3 = SU_N_mult(tmp3, tmp4 );
						tmp2 = SU_N_mult(tmp2,SU_N_conj(tmp3));	
						W_diag2[i] += SU_N_norm(tmp2);

						m1=index_up(x,nu);
						tmp2=SU_N_mult(F_N[link_index(mu,m0)],F_N[link_index(nu,index_up(m0,mu))]);
						tmp3=SU_N_mult(F_N[link_index(nu,x)],F_N[link_index(mu,m1)]);
						m1=index_up(m1,mu);
						
						tmp2 = SU_N_mult(tmp1, tmp2);
						tmp3 = SU_N_mult(tmp3, tmp4 );
						tmp2 = SU_N_mult(tmp2,SU_N_conj(tmp3));	
						W_diag2[i] += SU_N_norm(tmp2);

						m1=index_up(x,mu);
						tmp2=SU_N_mult(F_N[link_index(mu,m0)],F_N[link_index(nu,index_up(m0,mu))]);
						tmp3=SU_N_mult(F_N[link_index(mu,x)],F_N[link_index(nu,m1)]);
						m1=index_up(m1,mu);
						
						tmp2 = SU_N_mult(tmp1, tmp2);
						tmp3 = SU_N_mult(tmp3, tmp4 );
						tmp2 = SU_N_mult(tmp2,SU_N_conj(tmp3));	
						W_diag2[i] += SU_N_norm(tmp2);
						m1=index_up(x,nu);
						tmp2=SU_N_mult(F_N[link_index(nu,m0)],F_N[link_index(mu,index_up(m0,nu))]);
						tmp3=SU_N_mult(F_N[link_index(nu,x)],F_N[link_index(mu,m1)]);
						m1=index_up(m1,mu);
						
						tmp2 = SU_N_mult(tmp1, tmp2);
						tmp3 = SU_N_mult(tmp3, tmp4 );
						tmp2 = SU_N_mult(tmp2,SU_N_conj(tmp3));	
						W_diag2[i] += SU_N_norm(tmp2);
/*--------------------diag3--------------------------------------*/
						sigma=6-mu-nu;
						m1=index_up(x,sigma);
						tmp2=SU_N_mult(SU_N_mult(F_N[link_index(nu,m0)],F_N[link_index(mu,index_up(m0,nu))]),F_N[link_index(sigma,index_up(index_up(m0,nu),mu))]);
						tmp3=SU_N_mult(SU_N_mult(F_N[link_index(sigma,x)],F_N[link_index(mu,m1)]),F_N[link_index(nu,index_up(m1,mu))]);
						m1=index_up(index_up(m1,mu),nu);
						tmp4 = F_N[link_index(0,m1)];
						for(k=0;k<i;k++){
							m1=index_up(m1,0);
							tmp4=SU_N_mult(tmp4, F_N[link_index(0,m1)]);
						};
					//	if(index_up(index_up(m0,nu),mu)!= index_up(m1,0)) fprintf(stderr,"!!!!!!!!!!!!!!!");
						tmp2 = SU_N_mult(tmp1, tmp2);
						tmp3 = SU_N_mult(tmp3, tmp4 );
						tmp2 = SU_N_mult(tmp2,SU_N_conj(tmp3));	
						W_diag3[i] += SU_N_norm(tmp2);
					};
				};
					
};
		for(i=0;i<LS-2;i++){fprintf(stderr,"\n");
fprintf(pot,"\n");
			for(j=0;j<LS-2;j++){
				W[i][j]/=(double)(param->ipw[ param->D ]*(param->D-1)*N_RUNS);
				fprintf(stderr,"%f\t",W[i][j]);
				fprintf(pot,"%f\t",W[i][j]);
				if(j==0)
				{
					W_diag2[i]/=(double)(param->ipw[ param->D ]*(param->D-1)*(param->D-2)*N_RUNS*4.);
					fprintf(stderr,"%f\t",W_diag2[i]);
					fprintf(pot,"%f\t",W_diag2[i]);	
					W_diag3[i]/=(double)(param->ipw[ param->D ]*(param->D-1)*(param->D-2)*N_RUNS);
					fprintf(stderr,"%f\t",W_diag3[i]);
					fprintf(pot,"%f\t",W_diag3[i]);	
				};
		};};

//		};

	SU_N_free();
fclose(pot);
    return 0;
}
