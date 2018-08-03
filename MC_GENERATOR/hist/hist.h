#include<stdlib.h>
#include<stdio.h>

typedef struct{
 double x_min, x_max;
 int   n_bins;        
 double *bins;   // sum of all points which are in a given range of coords
 double *sbins;  // sum of squares
 int   *nbins;  // number of points in each bin
 int   npoints; // total number of points in a histogram       
} histogram;

typedef struct{
 double x_min, x_max, y_min, y_max;
 int   n_bins_x, n_bins_y;        
 double *bins;   // sum of all points which are in a given range of coords
 double *sbins;  // sum of squares
 int   *nbins;  // number of points in each bin
 int   npoints; // total number of points in a histogram       
} histogram2D;


histogram* init_hist(double x_min, double x_max, int n_bins);
void free_hist(histogram *hist);
void add_point(histogram *hist, double x, double s);

histogram2D* init_hist2D(double x_min, double x_max, double y_min, double y_max, int n_bins_x, int n_bins_y);
void free_hist2D(histogram2D *hist);
void add_point2D(histogram2D *hist, double x, double y, double s);
void plot_hist2D(histogram2D *hist, char *filename, int options);
void save_hist2D(histogram2D *hist, char *filename);

double max_in_hist(histogram *hist); // Searches for a maximal value in bins

void print_hist(histogram *hist);
void plot_hist_cons(histogram *hist);
void plot_hist  (histogram *hist, char *filename, int options); // filename is the name of the EPS file with output, option 0 - smooth without errbars, 1-errbars
void plot_2hist (histogram *hist1, histogram *hist2, char *filename, int options); // filename is the name of the EPS file with output, option 0 - smooth without errbars, 1-errbars
void save_hist  (histogram *hist, char *filename); // Saves hist such that it can be plotted with GNUplot

double lattice_spacing(double beta); // Two-loop corrections to RG equation, to use for scaling tests
void plot_2Darray_as_contour(double *z, int nx, int ny, double minx, double maxx, double miny, double maxy, char *filename); // Plots a contour plot using Gnuplot, z is assumed to be a [nx,ny] array
void save_2Darray(double *z, int nx, int ny, char *filename); // Saves array z as a table in  a specified text file

