#include<stdlib.h>
#include<stdio.h>

typedef struct{
 int nX;         // Total number of X values
 double* aX;     // Average X * n_points
 double* aX2;    // Average X^2 * n_points
 int     n_points;
 char    name[100];
 char    units[100];
 double  scale; // Rescale final averages by this factor
} value;

typedef struct{
 int nX;         // Total number of X values
 int nY;         // Total number of Y values
 double* aX;     // Average X * n_points
 double* aY;     // Average Y * n_points
 double* aX2;    // Average X^2 * n_points
 double* aY2;    // Average Y^2 * n_points
 double* aXY;    // Average XY
 int    n_points;
 char   name[100];
} correlator;

value*       init_value(int nX, char* name, char* units, double scale);
correlator*  init_correlator(int nX, int nY, char* name);

void free_value(value *v);
void free_correlator(correlator *c);

void add_point(value *v, double *x);

void add_pointc(correlator *c, double *x, double *y);

void print_value(value *v, FILE* f, char* token);
void print_correlator(correlator *c, FILE* f, char* token);
void get_value(value *v, double* mX, double* dX);
void get_correlator(correlator *c, double *cXY, double *dcXY);
