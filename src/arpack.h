/*!
 * \file arpack.h
 *
 * \brief C-interface/wrapper to ARPACK library
 *
 *
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \date   Чтв Сен 23 01:47:40 MSD 2004
 */
#ifndef _ARPACK_H_
#define _ARPACK_H_

#include "defs.h"
#include "types.h"


extern void UNDERSC(znaupd)
  (int*, char*, int*, char*, int*, double*, double complex*, 
   int*, double complex*, int*, int*, int*, double complex*, double complex*,
   int*, double* ,int*, int, int);
extern void UNDERSC(zneupd)
  (int*, char*, int*, double complex*, double complex*, int*, double complex*,
   double complex*, 
   char*, int*, char*, int*, double*, double complex*, 
   int*, double complex*, int*, int*, int*, double complex*, double complex*,
   int*, double* ,int*,
   int,int,int);
extern void UNDERSC(cnaupd)
  (int*, char*, int*, char*, int*, float*, float complex*, 
   int*, float complex*, int*, int*, int*, float complex*, float complex*,
   int*, float* ,int*, int, int);
extern void UNDERSC(cneupd)
  (int*, char*, int*, float complex*, float complex*, int*, float complex*,
   float complex*, 
   char*, int*, char*, int*, float*, float complex*, 
   int*, float complex*, int*, int*, int*, float complex*, float complex*,
   int*, float* ,int*,
   int,int,int);


#ifdef SINGLE_PREC
void naupd(int* a, char* b, int* c, char* d, int* e, float* f, t_complex * g,
		int* h, t_complex* i, int* j, int* k, int* l, t_complex* m,
		t_complex* n,	int* o, float* p,int* q, int r, int s);
void neupd(int* a, char* b, int* c, t_complex* d, t_complex* e, int* f,
		t_complex* g, t_complex* h,char* i, int* j, char* k, int* l,
		float* m, t_complex* n, int* o, t_complex* p, int* q, int* r,
		int* s, t_complex* t, t_complex* u, int* v, float* w, int* x,
		int y,int z ,int a1);
#else
void naupd(int* a, char* b, int* c, char* d, int* e, double* f, t_complex* g, 
		  int* h, t_complex* i, int* j, int* k, int* l, t_complex* m,
		  t_complex* n,	int* o, double* p,int* q, int r, int s);
void neupd(int* a, char* b, int* c, t_complex* d, t_complex* e, int* f,
		  t_complex* g, t_complex* h,char* i, int* j, char* k, int* l,
		  double* m, t_complex* n, int* o, t_complex* p, int* q, int* r,
		  int* s, t_complex* t, t_complex* u, int* v, double* w,int* x,
		  int y,int z ,int a1);
#endif

#endif
