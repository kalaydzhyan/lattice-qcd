#include<stdlib.h>
#include<stdio.h>
#include<unistd.h>
#include"MC-SU2.h"
#include"geometry.h"
#include"plotter.h"

int vol[3];

int main()
{
    vol[0] = 2; vol[1] = 1; vol[2] = 1;
    param_init(2, 4, 2, vol, 3);
    uchar cc1[3],cc2[3];
    uchar c;
        
    printf("%i \n",param->T1);
    
    for(uint x=0; x<6; x++)
     for(uint y = 0; y<4; y++)
     {
      cc1[0]=x; cc1[1]=y;
      uint m1 = site_index(cc1);
      uint m2r = index_up(m1,0,1);
      uint m2l = index_down(m1,0,1);
      uint m2u = index_up(m1,1,1);
      uint m2d = index_down(m1,1,1);
      site_coordinates(cc2, m2r);
      printf("from (%i,%i) right: (%i,%i) \n",x,y,cc2[0],cc2[1]);
      site_coordinates(cc2, m2l);
      printf("from (%i,%i) left: (%i,%i) \n",x,y,cc2[0],cc2[1]);
      site_coordinates(cc2, m2u);
      printf("from (%i,%i) up: (%i,%i) \n",x,y,cc2[0],cc2[1]);
      site_coordinates(cc2, m2d);
      printf("from (%i,%i) down: (%i,%i) \n",x,y,cc2[0],cc2[1]);
      scanf("%c",&c);
     }
    
    param_free();
	return 0;
}

/*  vol[0] = 1; vol[1] = 1; vol[2] = 1;
	fields_init(4,8,4,vol,2);
    param->beta = 2.35;
    double mp = MC_SU2();
    printf("mean plaquette: %f \n",mp);
    fields_free(); */

/*  printf("%i \n",param->size[0]);
    printf("%i \n",param->T1);
    double *mp, *bt;
    mp = (double *)calloc(5,sizeof(double));
    bt = (double *)calloc(5,sizeof(double));
    for(int nb=0; nb<5; nb++)
    {
	 param->beta = 2.5 - 0.1*(float)nb;
     for(int k = 0; k<10; k++)
      mp[nb] = MC_SU2();
     bt[nb] = param->beta;
     printf("%f %f \n",bt[nb],mp[nb]);
    }; 
    plot_array(bt,mp,5,"mp.eps");*/
