#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<statistics.h>

value*       init_value(int nX, char* name, char* units, double scale)
{
 int i;
 value *v = (value *)malloc(sizeof(value));

 v->aX     = (double *)malloc(nX*sizeof(double));
 for(i=0; i<nX; i++)
  v->aX[i] = 0.0;

 v->aX2    = (double *)malloc(nX*sizeof(double));
 for(i=0; i<nX; i++)
  v->aX2[i] = 0.0;

 v->n_points = 0;
 v->nX = nX;
 v->scale = scale;

 strcpy(v->name, name);
 strcpy(v->units, units);

 return v;
}

correlator*  init_correlator(int nX, int nY, char* name)
{
 int i;
 correlator *c = (correlator *)malloc(sizeof(correlator));

 c->aX  = (double *)malloc(nX*sizeof(double));
 for(i=0; i<nX; i++)
  c->aX[i] = 0.0;
 c->aX2 = (double *)malloc(nX*sizeof(double));
 for(i=0; i<nX; i++)
  c->aX2[i] = 0.0;

 c->aY  = (double *)malloc(nY*sizeof(double));
 for(i=0; i<nY; i++)
  c->aY[i] = 0.0;
 c->aY2 = (double *)malloc(nY*sizeof(double));
 for(i=0; i<nY; i++)
  c->aY2[i] = 0.0;

 c->aXY = (double *)malloc(nX*nY*sizeof(double));
 for(i=0; i<nX*nY; i++)
  c->aXY[i] = 0.0;


 c->nX = nX;
 c->nY = nY;
 c->n_points = 0;
 strcpy(c->name, name);

 return c;
}


void free_value(value *v)
{
 free(v->aX);
 free(v->aX2);
 free(v);
}

void free_correlator(correlator *c)
{
 free(c->aX);
 free(c->aY);
 free(c->aX2);
 free(c->aY2);
 free(c->aXY);
 free(c);
}

void add_point(value *v, double *x)
{
 int i;
 for(i=0; i<v->nX; i++)
 {
  v->aX[i]  += x[i];
  v->aX2[i] += x[i]*x[i];
 };
 v->n_points++;
}

void add_pointc(correlator *c, double *x, double *y)
{
 int i, j;

 for(i=0; i<c->nX; i++)
 {
  c->aX[i]  += x[i];
  c->aX2[i] += x[i]*x[i];
 };

 for(j=0; j<c->nY; j++)
 {
  c->aY[j]  += y[j];
  c->aY2[j] += y[j]*y[j];
 };

 for(i=0; i<c->nX; i++)
  for(j=0; j<c->nY; j++)
   c->aXY[i*c->nY + j] += x[i]*y[j];

 c->n_points++;
}

void print_value(value *v, FILE* f, char* token)
{
 int i;
 fprintf(f, "%s: %i points, ", token, v->n_points);
 for(i=0; i<v->nX; i++)
 {
  double mX  = v->aX[i]/(double)(v->n_points);
  double mX2 = v->aX2[i]/(double)(v->n_points);
  double dX  = sqrt(fabs( (mX2 - mX*mX)/((double)(v->n_points) - 1.0) ));
  fprintf(f, "%2.4E +/- %2.4E, ", v->scale*mX, v->scale*dX);
 };
 fprintf(f, "%s", v->units);
 fprintf(f, "\n\n");
 fflush(f);
}

void get_value(value *v, double* mX, double* dX)
{
 int i;
 for(i=0; i<v->nX; i++)
 {
  mX[i]  = v->aX[i]/(double)(v->n_points);
  double mX2 = v->aX2[i]/(double)(v->n_points);
  dX[i]  = sqrt(fabs( (mX2 - mX[i]*mX[i])/((double)(v->n_points) - 1.0) ));
  mX[i]  = mX[i]*v->scale;
  dX[i]  = dX[i]*v->scale;
 };
}

void print_correlator(correlator *c, FILE* f, char* token)
{
 int i,j;
 fprintf(f, "%s: %i points\n Correlation matrix: \n\n", token, c->n_points);

 double max_err = 0;

 for(i=0; i<c->nX; i++)
 {
  double mX  = c->aX[i]/(double)(c->n_points);
  double mX2 = c->aX2[i]/(double)(c->n_points);
  double sdX  = sqrt(mX2 - mX*mX);

  for(j=0; j<c->nY; j++)
  {
   double mY  = c->aY[j]/(double)(c->n_points);
   double mY2 = c->aY2[j]/(double)(c->n_points);
   double sdY  = sqrt(mY2 - mY*mY);

   double mXY = c->aXY[i*c->nY + j]/(double)(c->n_points);

   double  rXY = (mXY - mX*mY)/(sdX*sdY);
   double drXY = (1.0 - rXY*rXY)/sqrt((double)(c->n_points));
   max_err = (drXY>max_err)? drXY : max_err;

   fprintf(f, "%+02.4lf ", rXY);
  };
  fprintf(f, "\n");
 };
 fprintf(f, "\n Max. error: %2.4E \n", max_err);
 fflush(f);
}

void get_correlator(correlator *c, double *cXY, double *dcXY)
{
 int i,j;

 for(i=0; i<c->nX; i++)
 {
  double mX  = c->aX[i]/(double)(c->n_points);
  double mX2 = c->aX2[i]/(double)(c->n_points);
  double sdX  = sqrt(mX2 - mX*mX);

  for(j=0; j<c->nY; j++)
  {
   double mY  = c->aY[j]/(double)(c->n_points);
   double mY2 = c->aY2[j]/(double)(c->n_points);
   double sdY  = sqrt(mY2 - mY*mY);

   double mXY = c->aXY[i*c->nY + j]/(double)(c->n_points);

    cXY[i*c->nY + j] = (mXY - mX*mY)/(sdX*sdY);
   dcXY[i*c->nY + j] = (1.0 - cXY[i*c->nY + j]*cXY[i*c->nY + j])/sqrt((double)(c->n_points));
  };
 };
}

