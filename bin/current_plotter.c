#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<complex.h>

#define LS  	14
#define LT  	14
#define VOL 	(LS*LS*LS*LT)
#define sps     (LS*LS*LS)
#define ND  	4		 // Number of Dirac components
#define NC  	2		 // Number of colour components
#define BLOCK	(VOL*ND*NC + 1)
//#define csp     0.135225         // beta = 3.9410, improved action
#define csp       0.226044	   // beta = 3.2810, improved action
//#define csp   1.0
#define jsc     1000000.0

#define LS2     7
#define sps2    LS2*LS2*LS2

typedef double complex t_complex;

int dxc[4]   = {LS*LS*LS, LS*LS, LS,  1};
int size[4]  = {      LT,    LS, LS, LS};

int up(int x, int mu)
{
 int q = x;
 if((x/dxc[mu])%size[mu]==(size[mu]-1))
  q-=dxc[mu]*(size[mu]-1);
 else
  q+=dxc[mu];
 return q; 
};

int down(int x, int mu)
{
 int q = x;
 if((x/dxc[mu])%size[mu]==0)
  q+=dxc[mu]*(size[mu]-1);
 else
  q-=dxc[mu];
 return q; 
};

int xc2x(int xc[4])
{
 return dxc[0]*xc[0] + dxc[1]*xc[1] + dxc[2]*xc[2] + dxc[3]*xc[3];
};

void x2xc(int x, int *xc)
{
 for(int mu=0; mu<4; mu++)
  xc[mu] = (x/dxc[mu])%size[mu];
};

int main(int argc, char *argv[])
{
 if(argc<4) 
 {
  printf("Correct usage: current_plotter evn infile outfile\n");
  return EXIT_FAILURE;
 }; 
  
 //char* evreader = "/home/pools/2/lena/overlap_2/bin/read_ev_248s";
 //char* evreader = "/home/pools/1/buividovich/overlap/bin/read_ev_66s";
 char* evreader = "/home/pools/2/lena/overlap_14x14_ext/bin/read_ev_1414s";
 
 int iev;
 sscanf(argv[1], "%i", &iev);
 
 // Testing the up-down procedures
 
 /*int x, xc[4];
 
 for(int mu=0; mu<4; mu++)
 {
  xc[0]=7; xc[1]=0; xc[2]=0; xc[3]=23;
  x = xc2x(xc);
  x = up(x,mu);
  x2xc(x,xc);
  printf("up(x,%i): %i %i %i %i\n", mu, xc[0], xc[1], xc[2], xc[3]);
 
  xc[0]=7; xc[1]=0; xc[2]=0; xc[3]=23;
  x = xc2x(xc);
  x = down(x,mu);
  x2xc(x,xc);
  printf("down(x,%i): %i %i %i %i\n", mu, xc[0], xc[1], xc[2], xc[3]);
 };*/
  
 char cmd[500];
 strcpy(cmd, evreader);
 strcat(cmd, " -i ");
 strcat(cmd, argv[2]);
 strcat(cmd, " -o evtmp1.dat");
 
 printf("Started processing the file %s, eigenvalue number is %i ...\n", argv[2], iev);
 
 system(cmd);
 system("mkdir -p ./jslices");
 
 double Re, Im;
 t_complex X;
 
 t_complex *psi = (t_complex *)malloc(VOL*NC*ND*sizeof(t_complex));
 if(!psi)
 {
  printf("Cannot allocate memory for psi!!!");
  return EXIT_FAILURE;
 };
 
 double    *j   = (double *)malloc(VOL*4*sizeof(double));
 if(!j)
 {
  printf("Cannot allocate memory for j!!!");
  return EXIT_FAILURE;
 };
 for(int i=0; i<VOL*4; i++)
  j[i] = 0.0;
 
 int counter   = 0;
 int evcounter = 0;
 
 FILE *infile = fopen("evtmp1.dat", "r");
 
 while(fscanf(infile, "%lf %lf\n", &Re, &Im)!=EOF)
 {
  X = Re + I*Im;
  if((counter % BLOCK) == 0)
  {
   double dnorm  = fabs(cabs(X - 1.0) - 1);
   double lambda = fabs(2.0*440.0*Im/(csp*(2.0 - Re)));
   printf("counter: %i, evcounter: %i, BLOCK: %i, lambda: %4.4lf\n", counter, evcounter, BLOCK, lambda);
   if(dnorm>0.0001)
   {
    printf("Bad eigenvalue Nr. %i, dnorm = %2.4lf: lies off the circle\n", evcounter + 1, dnorm);
    printf("Re: %2.4lf Im: %2.4lf delta: %2.4lf\n", Re, Im, dnorm);
   };
   evcounter++;
  }
  else
   if(evcounter == iev)
   {
    int ixc = (counter % BLOCK) - 1;
    psi[ixc] = Re + I*Im;
   }; 
  counter++;
 };
 
 fclose(infile);
 
 // Here we have all the components stored in an array psi

 t_complex migkg0[4][ND][ND] = { // minus I times gamma_0 times gamma_k, k = 1,2,3 - spatial indices
  {{   1.0,   0.0,    0.0,    0.0},
   {   0.0,   1.0,    0.0,    0.0},
   {   0.0,   0.0,    1.0,    0.0},
   {   0.0,   0.0,    0.0,    1.0}}, 
    
  {{   0.0,   1.0,    0.0,    0.0},
   {   1.0,   0.0,    0.0,    0.0},
   {   0.0,   0.0,    0.0,   -1.0},
   {   0.0,   0.0,   -1.0,    0.0}},
   
  {{   0.0,-1.0*I,    0.0,    0.0},
   { 1.0*I,   0.0,    0.0,    0.0},
   {   0.0,   0.0,    0.0,  1.0*I},
   {   0.0,   0.0, -1.0*I,    0.0}},
   
  {{   1.0,   0.0,    0.0,    0.0},
   {   0.0,  -1.0,    0.0,    0.0},
   {   0.0,   0.0,   -1.0,    0.0},
   {   0.0,   0.0,    0.0,    1.0}}
   
 };
 
 t_complex gamma[4][ND][ND] = { // Euclidean Gamma-matrices in the chiral represnetation
  {{   0.0,    0.0,    1.0,    0.0},
   {   0.0,    0.0,    0.0,    1.0},
   {   1.0,    0.0,    0.0,    0.0},
   {   0.0,    1.0,    0.0,    0.0}}, 
    
  {{   0.0,    0.0,    0.0, -1.0*I},
   {   0.0,    0.0, -1.0*I,    0.0},
   {   0.0,  1.0*I,    0.0,    0.0},
   { 1.0*I,    0.0,    0.0,    0.0}},
   
  {{   0.0,    0.0,    0.0,   -1.0},
   {   0.0,    0.0,    1.0,    0.0},
   {   0.0,    1.0,    0.0,    0.0},
   {  -1.0,    0.0,    0.0,    0.0}},
   
  {{   0.0,    0.0, -1.0*I,    0.0},
   {   0.0,    0.0,    0.0,  1.0*I},
   { 1.0*I,    0.0,    0.0,    0.0},
   {   0.0, -1.0*I,    0.0,    0.0}}
   
 };
 
 char outfname[300], genfname[300];
 
 double density_norm = 0.0;
 
 for(int sidx = 0; sidx < LT; sidx++)
  for(int x = 0; x<sps; x++)
  {
   int idx = sps*sidx + x;
   
   for(int ic = 0; ic<NC; ic++)
    for(int id1 = 0; id1<ND; id1++)
    {
     density_norm += creal(conj(psi[NC*ND*idx + ND*ic + id1])*psi[NC*ND*idx + ND*ic + id1]);
     for(int id2 = 0; id2<ND; id2++)
      for(int k = 0; k<4; k++)
       j[4*sps*sidx + 4*x + k] += creal(conj(psi[NC*ND*idx + ND*ic + id1])*migkg0[k][id1][id2]*psi[NC*ND*idx + ND*ic + id2]);
    };
  };
  
 // calculating the divergence of j
 
 double dj2  = 0;
 double aj2  = 0;
  
 for(int x = 0; x<VOL; x++)
  {
   double dj = 0;
   for(int mu = 0; mu<4; mu++)
   {
    dj += 0.5*(j[4*up(x,mu) + mu] - j[4*down(x,mu) + mu]);
    aj2+=j[4*x + mu]*j[4*x + mu];
   };    
   dj2+=dj*dj;
  };
  dj2 = sqrt(dj2/(double)VOL);
  aj2 = sqrt(aj2/(double)VOL);
 
  printf("Average square of the divergence is: %8.8lf\n",dj2);
  printf("Average square of the current    is: %8.8lf\n",aj2);
  
 // calculating the dispersions of j
 
 double a2j[4] = {0, 0, 0, 0};
 double a1j[4] = {0, 0, 0, 0};
 int    naj    = 0;
 
 for(int sidx = 0; sidx<LT; sidx++)
  for(int x = 0; x<sps; x++)
   for(int k = 0; k<4; k++)
   {
    a1j[k] += jsc*j[4*sps*sidx + 4*x + k];
    a2j[k] += jsc*jsc*j[4*sps*sidx + 4*x + k]*j[4*sps*sidx + 4*x + k];
    naj    ++;
   }; 
   
 printf("integral of j0: %4.10lf\n", a1j[0]);  
   
 for(int k = 0; k<4; k++)
  printf("%i: a1j = %4.10lf, a2j = %4.10lf\n", k, a1j[k]/(double)naj, a2j[k]/(double)naj);
   
 double j_blocked[LT][sps2][4];
 
 for(int sidx = 0; sidx < LT; sidx++)
  for(int xb1 = 0; xb1 < LS2; xb1++)
   for(int xb2 = 0; xb2 < LS2; xb2++)
    for(int xb3 = 0; xb3 < LS2; xb3++)
    {
     int xb = xb1 + LS2*xb2 + LS2*LS2*xb3;
     
     int x1 = 2*xb1;
     int x2 = 2*xb2;
     int x3 = 2*xb3;
     int x = x1 + LS*x2 + LS*LS*x3;
     j_blocked[sidx][xb][0] += 0.25*j[4*sps*sidx + 4*x + 0];
     j_blocked[sidx][xb][1] += 0.25*j[4*sps*sidx + 4*x + 1];
     j_blocked[sidx][xb][2] += 0.25*j[4*sps*sidx + 4*x + 2];
     j_blocked[sidx][xb][3] += 0.25*j[4*sps*sidx + 4*x + 3];
     
     x1 = 2*xb1 + 1;
     x2 = 2*xb2;
     x3 = 2*xb3;
     x  = x1 + LS*x2 + LS*LS*x3;
     j_blocked[sidx][xb][0] += 0.25*j[4*sps*sidx + 4*x + 0];
     j_blocked[sidx][xb][1] += 0.25*j[4*sps*sidx + 4*x + 1];
     j_blocked[sidx][xb][2] += 0.25*j[4*sps*sidx + 4*x + 2];
     j_blocked[sidx][xb][3] += 0.25*j[4*sps*sidx + 4*x + 3];
     
     x1 = 2*xb1;
     x2 = 2*xb2 + 1;
     x3 = 2*xb3;
     x  = x1 + LS*x2 + LS*LS*x3;
     j_blocked[sidx][xb][0] += 0.25*j[4*sps*sidx + 4*x + 0];
     j_blocked[sidx][xb][1] += 0.25*j[4*sps*sidx + 4*x + 1];
     j_blocked[sidx][xb][2] += 0.25*j[4*sps*sidx + 4*x + 2];
     j_blocked[sidx][xb][3] += 0.25*j[4*sps*sidx + 4*x + 3];
     
     
     x1 = 2*xb1;
     x2 = 2*xb2;
     x3 = 2*xb3 + 1;
     x  = x1 + LS*x2 + LS*LS*x3;
     j_blocked[sidx][xb][0] += 0.25*j[4*sps*sidx + 4*x + 0];
     j_blocked[sidx][xb][1] += 0.25*j[4*sps*sidx + 4*x + 1];
     j_blocked[sidx][xb][2] += 0.25*j[4*sps*sidx + 4*x + 2];
     j_blocked[sidx][xb][3] += 0.25*j[4*sps*sidx + 4*x + 3];
    };
 
 for(int sidx = 0; sidx < LT; sidx++)
 {
  sprintf(outfname, "./jslices/js%i.dat", sidx);
  sprintf(genfname, "./jslices/js%i.general", sidx);
  printf("outfname: %s, slice: %i\n", outfname, sidx);

  FILE* genfile = fopen(genfname, "w");
  fprintf(genfile, "file = js%i.dat\n", sidx);
  fprintf(genfile, "grid = %i x %i x %i\n", LS2, LS2, LS2);
  fprintf(genfile, "format = ascii\n");
  fprintf(genfile, "interleaving = field\n");
  fprintf(genfile, "majority = row\n");
  fprintf(genfile, "field = field0\n");
  fprintf(genfile, "structure = 3-vector\n");
  fprintf(genfile, "type = float\n");
  fprintf(genfile, "dependency = positions\n");
  fprintf(genfile, "positions = regular, regular, regular, 0, 1, 0, 1, 0, 1\n");
  fprintf(genfile, "\n");
  fprintf(genfile, "end\n");
  fclose(genfile);

  FILE *datfile = fopen(outfname, "w");
  for(int x = 0; x<sps2; x++)
   fprintf(datfile, "%4.2lf %4.2lf %4.2lf\n", jsc*j_blocked[sidx][x][1], jsc*j_blocked[sidx][x][2], jsc*j_blocked[sidx][x][3]);

  fclose(datfile);
 };
 
 printf("Norm of the density: %4.6lf\n", density_norm);
 
 strcpy(cmd, "zip -r -j ");
 strcat(cmd, argv[3]);
 strcat(cmd, " ./jslices");

 //system("rm -f ./evtmp.dat");
 system(cmd);
 system("rm -r -f ./jslices");
 
 free(psi);
 free(j);

 return EXIT_SUCCESS;
};
