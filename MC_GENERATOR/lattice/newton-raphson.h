#ifndef _NEWTON_RAPHSON_H_
#define _NEWTON_RAPHSON_H_

#include <geometry.h>
#include <lapack.h>

/* Convergence criterion on \delta x. */
#define NR_TOLX   FLT_EPSILON

/*  Ensures sufficient decrease in function value. */
#define NR_ALF     1.0e-4

/* Maximum number of iterations */
#define NR_MAXITS 200

/* Sets the convergence criterion on function values */
//#define NR_TOLF 1.0e-8
#define NR_TOLF FLT_EPSILON

/* Sets the criterion for deciding whether spurious convergence to a minimum has occurred */
//#define NR_TOLMIN 1.0e-6
#define NR_TOLMIN   FLT_EPSILON



/* Given an initial guess x for a root in n dimensions, find the root by a globally convergent
   Newton-Raphson method. The vector of functions \vec{F} to be zeroed is given by the user-supplied
   callback. Return status is \pm 1 (sign of Jacobian at root) on a normal termination, otherwise
   return status is 0 indicating that routine has converged to a local minimum of the function
   \vec{F} \vec{F}. In this case try restarting from a different initial guess. */
extern int newton_raphson( uint n, double *x,
			   double *low, double *high,
			   void *data, int (*callback)( double *, void *, double *),
			   int *status );

/* ---------------------------------------------------------------------------- */
/* The preferred way to solve the linear set of equations A x = b :
       double *a,*b, det;
       int    n,*indx;
       ...
       ludcmp(a,n,indx, &det);
       lubksb(a,n,indx,b);
   The answer x will be given back in b. Your original matrix A will have been
   destroyed. If you subsequently want to solve a set of equations with the same A
   but a different right-hand side b, you repeat only
       lubksb(a,n,indx,b);
   not, of course, with the original matrix A, but with a and indx as were already
   set by ludcmp. */

/* Zero return status indicates that input matrix is singular - solution is NOT available.
   Optional det is matrix determinant on output (if you don't need it pass NULL to ludcmp).
   Normal return status is \pm 1 - sign of det(A) */
extern int  ludcmp( double *A, uint n, uint *indx, double *det  );
extern void lubksb( double  *A, uint n, uint *indx, double *b );

/* ---------------------------------------------------------------------------- */

#define NR_R             0.61803399
#define NR_C             (1.0 - NR_R)
#define NR_PRES          FLT_EPSILON

/* Golden section search of maximum (sign=1) / minimum (sign=-1) of given function
   in the interval [a,c]. Triple a,b,c MUST bracket extremum. Returns function's
   extremal value, 'double *t' points to extremum location on exit. */
extern double find_extremum( double a, double b, double c, double *t, int sign,
			     void *data, double (*callback)( double , void * ) );

/* Finds root of the function 'callback' on the interval [start ; end] (with known
   function's values F_start, F_end). Root MUST be located in the above interval.
   Returns 'almost zero' value of function at the root and sets pointer *t accordingly. */
extern double find_root( double start, double F_start, double end,   double F_end, double *t,
			 void *data, double (*callback)( double, void * ) );
#endif
