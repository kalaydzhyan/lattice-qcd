#include<stdlib.h>
#include<stdio.h>
#include<strings.h>
#include<math.h>
#include"hist.h"

histogram* init_hist(double x_min, double x_max, int n_bins)
{
 histogram *hist;
 hist = (histogram* )malloc(sizeof(histogram));
 hist->bins = (double* )malloc(n_bins*sizeof(double));
 for(int i=0; i<n_bins; i++)
  hist->bins[i]=0;
 hist->sbins = (double* )malloc(n_bins*sizeof(double));
 for(int i=0; i<n_bins; i++)
  hist->sbins[i]=0;
 hist->nbins = (int* )malloc(n_bins*sizeof(int));
 for(int i=0; i<n_bins; i++)
  hist->nbins[i]=0;
 hist->n_bins = n_bins;
 hist->x_min = x_min;
 hist->x_max = x_max;
 hist->npoints=0;
 return hist; 
}

void free_hist(histogram *hist)
{
 free(hist->bins);
 free(hist->sbins);
 free(hist->nbins);
 free(hist);
}

void add_point(histogram *hist, double x, double s)
{
 if(x>=hist->x_min && x<=hist->x_max)
 {
  double q = (x - hist->x_min)/(hist->x_max - hist->x_min);
  int i=(int)((double)hist->n_bins*q);
  hist->nbins[i]++;
  hist->bins[i]+=s;
  hist->sbins[i]+=s*s;
  hist->npoints++; 
 };
}

void print_hist(histogram *hist)
{
 double x1,x2,avs;
 printf("Histogram: \n");
 printf("Number of points: %i \n",hist->npoints);
 for(int i=0; i<hist->n_bins; i++)
 {
  avs = (hist->nbins[i]!=0)? hist->bins[i]/hist->nbins[i] : 0;       
  x1 = hist->x_min + i/((double)hist->n_bins)*(hist->x_max - hist->x_min);
  x2 = hist->x_min + (i+1)/((double)hist->n_bins)*(hist->x_max - hist->x_min);       
  printf("x in [%2.6f , %2.6f]: %i points, <s(x)> = %2.6f \n",x1,x2, hist->nbins[i],avs);       
 };
 printf("\n \n");
}

void save_hist (histogram *hist, char *filename) // Saves the histogram in a form which is used by gnuplot
{
 FILE *data;
 double avs,avsq,x,err;
 double avn; // These vars are used to find the distribution of x's
 int tvn = 0;
 data = fopen(filename,"w");
 for(int i=0; i<hist->n_bins; i++)
  tvn+=hist->nbins[i];
 for(int i=0; i<hist->n_bins; i++)
 {
  avs  = (hist->nbins[i]!=0)? hist->bins[i]/hist->nbins[i] : 0;
  avn  = (tvn!=0)? ((double)hist->nbins[i]/(double)tvn) : 0;
  avsq = (hist->nbins[i]!=0)? hist->sbins[i]/hist->nbins[i] : 0;
  err  = (hist->nbins[i]!=0)? sqrt((avsq - avs*avs)/hist->nbins[i]) : 0;
  x    = hist->x_min + ((double)i+0.5)/((double)hist->n_bins)*(hist->x_max - hist->x_min);       
  fprintf(data,"%2.6f %2.6f %2.6f %2.6f\n",x,avs,err,avn);
 };
 fclose(data); 
}

void plot_hist (histogram *hist, char *filename, int options)
{
 FILE *gnup, *data;
 double avs,avsq,x,err;
 double avn; // These vars are used to find the distribution of x's
 int tvn = 0;
 gnup = popen("gnuplot - ","w");
 data = fopen("data123tmp.txt","w");
 for(int i=0; i<hist->n_bins; i++)
  tvn+=hist->nbins[i];
 for(int i=0; i<hist->n_bins; i++)
 {
  avs  = (hist->nbins[i]!=0)? hist->bins[i]/hist->nbins[i] : 0;
  avn  = (tvn!=0)? ((double)hist->nbins[i]/(double)tvn) : 0;
  avsq = (hist->nbins[i]!=0)? hist->sbins[i]/hist->nbins[i] : 0;
  err  = (hist->nbins[i]!=0)? sqrt((avsq - avs*avs)/hist->nbins[i]) : 0;
  x    = hist->x_min + ((double)i+0.5)/((double)hist->n_bins)*(hist->x_max - hist->x_min);
  if(options==2)       
   fprintf(data,"%2.6f %2.6f %2.6f \n",x,avs,err);
  if(options==0) 
   fprintf(data,"%2.6f %2.6f \n",x,avs); 
  if(options==1) 
   fprintf(data,"%2.6f %2.6f \n",x,avn);  
 };
 fclose(data); 
 printf("gnuplot pipe opened \n");
 fprintf(gnup,"set term postscript color enhanced solid landscape \"Helvetica\" 12\n");
 fprintf(gnup,"set out '%s'\n",filename);
 
 if(options==2)
  fprintf(gnup,"plot [%2.2f : %2.2f ] 'data123tmp.txt' with yerrorbars\n", hist->x_min - 0.1*(hist->x_max - hist->x_min), hist->x_max + 0.1*(hist->x_max - hist->x_min));
 if(options==0)
  fprintf(gnup,"plot [%2.2f : %2.2f ] 'data123tmp.txt' using 1:2:(1.0) smooth acsplines\n", hist->x_min - 0.1*(hist->x_max - hist->x_min), hist->x_max + 0.1*(hist->x_max - hist->x_min));
 if(options==1)
  fprintf(gnup,"plot [%2.2f : %2.2f ] 'data123tmp.txt' using 1:2 notitle smooth csplines\n", hist->x_min - 0.1*(hist->x_max - hist->x_min), hist->x_max + 0.1*(hist->x_max - hist->x_min));
  
 fprintf(gnup,"q\n");
 pclose(gnup);  
 system("rm data123tmp.txt");     
}

void plot_2hist (histogram *hist1, histogram *hist2, char *filename, int options)
{
 FILE *gnup, *data;
 double avs1,avsq1,x1,err1;
 double avs2,avsq2,x2,err2;
 gnup = popen("gnuplot - ","w");
 data = fopen("data123tmp.txt","w");
 for(int i=0; i<hist1->n_bins; i++)
 {
  avs1  = (hist1->nbins[i]!=0)? hist1->bins[i]/hist1->nbins[i] : 0;
  avsq1 = (hist1->nbins[i]!=0)? hist1->sbins[i]/hist1->nbins[i] : 0;
  err1  = (hist1->nbins[i]!=0)? sqrt((avsq1 - avs1*avs1)/hist1->nbins[i]) : 0;
  x1    = hist1->x_min + i/((double)hist1->n_bins)*(hist1->x_max - hist1->x_min);
  
  avs2  = (hist2->nbins[i]!=0)? hist2->bins[i]/hist2->nbins[i] : 0;
  avsq2 = (hist2->nbins[i]!=0)? hist2->sbins[i]/hist2->nbins[i] : 0;
  err2  = (hist2->nbins[i]!=0)? sqrt((avsq2 - avs2*avs2)/hist2->nbins[i]) : 0;
  x2    =  hist2->x_min + i/((double)hist2->n_bins)*(hist2->x_max - hist2->x_min);
  
  if( (hist1->nbins[i]!=0) && (hist2->nbins[i]!=0) )
  {
   if(options!=0)       
    fprintf(data,"%2.6f %2.6f %2.6f %2.6f %2.6f %2.6f \n",x1,avs1,err1,x2,avs2,err2);
   else
    fprintf(data,"%2.6f %2.6f %2.6f %2.6f \n",x1,avs1,x2,avs2); 
  };  
 };
 fclose(data); 
 printf("gnuplot pipe opened \n");
 fprintf(gnup,"set term postscript color enhanced solid landscape \"Helvetica\" 12\n");
 fprintf(gnup,"set out '%s'\n",filename);
 double lrx=hist1->x_min - 0.1*(hist1->x_max - hist1->x_min); // Upper and lower plotting limits
 double urx=hist1->x_max + 0.1*(hist1->x_max - hist1->x_min);
 if(options!=0)
  fprintf(gnup,"plot [%2.2f : %2.2f ] 'data123tmp.txt' using 1:2:3 with yerrorbars title \"hist1\",'data123tmp.txt' using 4:5:6 with yerrorbars title \"hist2\", 'data123tmp.txt' using 1:2 smooth csplines notitle,  'data123tmp.txt' using 4:5 smooth csplines notitle \n", lrx, urx);
 else 
  fprintf(gnup,"plot [%2.2f : %2.2f ] 'data123tmp.txt' using 1:2:(1.0) smooth acsplines,  'data123tmp.txt' using 3:4:(1.0) smooth acsplines \n", lrx, urx);

 fprintf(gnup,"q\n");
 pclose(gnup);  
 system("rm data123tmp.txt");        
}

histogram2D* init_hist2D(double x_min, double x_max, double y_min, double y_max, int n_bins_x, int n_bins_y)
{
 histogram2D *hist;
 hist = (histogram2D* )malloc(sizeof(histogram2D));
 hist->bins = (double* )malloc(n_bins_x*n_bins_y*sizeof(double));
 for(int i=0; i<n_bins_x; i++)
  for(int j=0; j<n_bins_y; j++)
   hist->bins[n_bins_y*i + j]=0;
 hist->sbins = (double* )malloc(n_bins_x*n_bins_y*sizeof(double));
 for(int i=0; i<n_bins_x; i++)
  for(int j=0; j<n_bins_y; j++)
   hist->sbins[n_bins_y*i + j]=0;
 hist->nbins = (int* )malloc(n_bins_x*n_bins_y*sizeof(int));
 for(int i=0; i<n_bins_x; i++)
  for(int j=0; j<n_bins_y; j++)
   hist->nbins[n_bins_y*i + j]=0;
 hist->n_bins_x = n_bins_x;
 hist->n_bins_y = n_bins_y;
 hist->x_min = x_min;
 hist->x_max = x_max;
 hist->y_min = y_min;
 hist->y_max = y_max;

 hist->npoints=0;
 return hist; 
}

void free_hist2D(histogram2D *hist)
{
 free( hist->bins);    
 free(hist->sbins);
 free(hist->nbins);
 free(hist);
}

void add_point2D(histogram2D *hist, double x, double y, double s)
{
 if((x>=hist->x_min) && (x<=hist->x_max) && (y>=hist->y_min) && (y<=hist->y_max))
 {
  double qx = (x - hist->x_min)/(hist->x_max - hist->x_min);
  double qy = (y - hist->y_min)/(hist->y_max - hist->y_min);           
  int i=(int)((double)hist->n_bins_x*qx);
  int j=(int)((double)hist->n_bins_y*qy);
  hist->nbins[hist->n_bins_y*i + j]++;
  hist->bins[hist->n_bins_y*i + j]+=s;
  hist->sbins[hist->n_bins_y*i + j]+=s*s;
  hist->npoints++;           
 };               
}

void plot_hist2D(histogram2D *hist, char *filename, int options)
{
 FILE *gnup, *data;
 double avs,avsq,x,y,err;
 gnup = popen("gnuplot - ","w");
 data = fopen("data123tmp.txt","w");
 for(int i=0; i<hist->n_bins_x; i++)
 {
  for(int j=0; j<hist->n_bins_y; j++)
  {
   avs  = (hist->nbins[i*hist->n_bins_y + j]!=0)? hist->bins[i*hist->n_bins_y + j]/hist->nbins[i*hist->n_bins_y + j] : 0;
   avsq = (hist->nbins[i*hist->n_bins_y + j]!=0)? hist->sbins[i*hist->n_bins_y + j]/hist->nbins[i*hist->n_bins_y + j] : 0;
   err  = (hist->nbins[i*hist->n_bins_y + j]!=0)? sqrt((avsq - avs*avs)/hist->nbins[i*hist->n_bins_y + j]) : 0;
   x    = hist->x_min + i/((double)hist->n_bins_x)*(hist->x_max - hist->x_min);
   y    = hist->y_min + j/((double)hist->n_bins_y)*(hist->y_max - hist->y_min);
   fprintf(data,"%2.6f %2.6f %2.6f \n",x,y,avs);
  };
  fprintf(data,"\n");
 }; 
 fclose(data); 
 printf("gnuplot pipe opened \n");
 fprintf(gnup,"set term postscript color enhanced solid landscape \"Helvetica\" 12\n");
 fprintf(gnup,"set out '%s'\n",filename);
 fprintf(gnup,"splot 'data123tmp.txt' with lines\n");
 fprintf(gnup,"q\n");
 pclose(gnup);  
 system("rm data123tmp.txt");     
}

void save_hist2D(histogram2D *hist, char *filename)
{
 FILE *data;
 double avs,avsq,x,y,err;
 data = fopen(filename,"w");
 for(int i=0; i<hist->n_bins_x; i++)
 {
  for(int j=0; j<hist->n_bins_y; j++)
  {
   avs  = (hist->nbins[i*hist->n_bins_y + j]!=0)? hist->bins[i*hist->n_bins_y + j]/hist->nbins[i*hist->n_bins_y + j] : 0;
   avsq = (hist->nbins[i*hist->n_bins_y + j]!=0)? hist->sbins[i*hist->n_bins_y + j]/hist->nbins[i*hist->n_bins_y + j] : 0;
   err  = (hist->nbins[i*hist->n_bins_y + j]!=0)? sqrt((avsq - avs*avs)/hist->nbins[i*hist->n_bins_y + j]) : 0;
   x    = hist->x_min + i/((double)hist->n_bins_x)*(hist->x_max - hist->x_min);
   y    = hist->y_min + j/((double)hist->n_bins_y)*(hist->y_max - hist->y_min);
   fprintf(data,"%2.6f %2.6f %2.6f %2.6f \n",x,y,avs,err);
  };
  fprintf(data,"\n");
 }; 
 fclose(data);
}

double lattice_spacing(double beta)
{
 static const double rbeta[16]    = {  2.25,   2.30,   2.35,   2.40,   2.45,  2.475,   2.50, 2.5115,   2.55,   2.60,   2.635,   2.65,   2.70,    2.74,    2.85,   2.96};  
 static const double rspacing[16] = {0.1898, 0.1655, 0.1394, 0.1193, 0.0996, 0.0904, 0.0854, 0.0823, 0.0713, 0.0601,  0.0541, 0.0532, 0.0455, 0.04086,  0.0283, 0.0222};      
 double q = 0;
 if(beta< rbeta[0])
  q = rspacing[0];
 if(beta>=rbeta[15])
  q = rspacing[15]; 
 for(int l=0; l<15; l++)
  if(beta<rbeta[l+1] && beta>=rbeta[l])
   q = rspacing[l] +  (rspacing[l+1] - rspacing[l])*(beta - rbeta[l])/(rbeta[l+1] - rbeta[l]);
 return q;    
}

double max_in_hist(histogram *hist) // Searches for a maximal value in bins
{
 double mx = -1E20;
 for(int i=0; i<hist->n_bins; i++)
  if(hist->bins[i]>=mx)
   mx=hist->bins[i];
 return mx;   
};

void plot_hist_cons(histogram *hist)
{
 double mx = max_in_hist(hist);
 printf(" \n Histogram: \n \n");
 for(int i=0; i<hist->n_bins; i++)
 {
  int k = roundf(20.0*hist->bins[i]/mx);
  for(int l=0; l<k; l++)
   printf("#");
  printf("\n"); 
 };
};

void plot_2Darray_as_contour(double *z, int nx, int ny, double minx, double maxx, double miny, double maxy, char *filename)
{
 FILE *gnup, *data;
 gnup = popen("gnuplot - ","w");
 data = fopen("data123tmp.txt","w");

 for(int i=0; i<nx; i++)
 {
  for(int j=0; j<ny; j++)
  {
   double x = minx + (maxx - minx)/((double)nx - 1.0)*i; 
   double y = miny + (maxy - miny)/((double)ny - 1.0)*j;
   fprintf(data,"%2.6f %2.6f %2.6f \n",x,y,z[i + nx*j]);
  };
  fprintf(data,"\n");
 };
 fclose(data); 

 printf("gnuplot pipe opened \n");
 fprintf(gnup,"set term postscript color enhanced solid landscape \"Helvetica\" 12\n");
 fprintf(gnup,"set out '%s'\n",filename);

 fprintf(gnup,"set pm3d map\n");

 fprintf(gnup,"splot 'data123tmp.txt' with lines\n");

 fprintf(gnup,"q\n");
 pclose(gnup);  
 system("rm data123tmp.txt");     
};

void save_2Darray(double *z, int nx, int ny, char *filename)
{
 FILE *data;
 data = fopen(filename,"w");
 for(int i=0; i<nx; i++)
 {
  for(int j=0; j<ny; j++)
   fprintf(data,"%2.6f ",z[i + nx*j]);
  fprintf(data,"\n");
 };
 fclose(data); 
};

