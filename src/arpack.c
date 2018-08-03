#include "types.h"
#include "arpack.h"

#ifdef SINGLE_PREC
void naupd(int* a, char* b, int* c, char* d, int* e, float* f, t_complex * g,
		int* h, t_complex* i, int* j, int* k, int* l, t_complex* m,
		t_complex* n,	int* o, float* p,int* q, int r, int s)
{
  UNDERSC(cnaupd)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s);
}

void neupd(int* a, char* b, int* c, t_complex* d, t_complex* e, int* f,
		t_complex* g, t_complex* h,char* i, int* j, char* k, int* l,
		float* m, t_complex* n, int* o, t_complex* p, int* q, int* r,
		int* s, t_complex* t, t_complex* u, int* v, float* w, int* x,
		int y,int z ,int a1)
{
  UNDERSC(cneupd)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a1);
}
#else
void naupd(int* a, char* b, int* c, char* d, int* e, double* f, t_complex* g, 
		  int* h, t_complex* i, int* j, int* k, int* l, t_complex* m,
		  t_complex* n,	int* o, double* p,int* q, int r, int s)
{
  UNDERSC(znaupd)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s);
}

void neupd(int* a, char* b, int* c, t_complex* d, t_complex* e, int* f,
		  t_complex* g, t_complex* h,char* i, int* j, char* k, int* l,
		  double* m, t_complex* n, int* o, t_complex* p, int* q, int* r,
		  int* s, t_complex* t, t_complex* u, int* v, double* w,int* x,
		  int y,int z ,int a1)
{
  UNDERSC(zneupd)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a1);
}
#endif

