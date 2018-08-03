#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define LS  	16
#define LT  	16
#define VOL 	(LS*LS*LS*LT)
#define LS2   LS 

float allmas[LT][LS][LS][LS];

int x,y,z,t;

int main(int argc, char *argv[])
{
  
  char* datafilename = "./instanton_chirality.dat";

  FILE *datfile = fopen(datafilename, "r");
  for(t = 0; t<LT; t++)				
  	  for(x = 0; x<LS; x++)
		for(y = 0; y<LS; y++)
			for(z = 0; z<LS; z++)
				fscanf(datfile, "%e", &allmas[t][x][y][z]);
  fclose(datfile);

 char cmd[500];

/*
 strcpy(cmd, evreader);
 strcat(cmd, " -i ");
 strcat(cmd, argv[2]);
 strcat(cmd, " -o evtmp1.dat");
 
 printf("Started processing the file %s, eigenvalue number is %i ...\n", argv[2], iev);
 
 system(cmd);

*/

 system("mkdir -p ./jslices");

 char outfname[300], genfname[300];
 int sidx;

 for(sidx = 0; sidx < LT; sidx++)
 {
  sprintf(outfname, "./jslices/js%i.dat", sidx);
  sprintf(genfname, "./jslices/js%i.general", sidx);
//  printf("outfname: %s, slice: %i\n", outfname, sidx);

  FILE* genfile = fopen(genfname, "w");
  fprintf(genfile, "file = js%i.dat\n", sidx);
  fprintf(genfile, "grid = %i x %i x %i\n", LS2, LS2, LS2);
  fprintf(genfile, "format = ascii\n");
  fprintf(genfile, "interleaving = field\n");
  fprintf(genfile, "majority = row\n");
  fprintf(genfile, "field = field0\n");
  fprintf(genfile, "structure = scalar\n");
  fprintf(genfile, "type = float\n");
  fprintf(genfile, "dependency = positions\n");
  fprintf(genfile, "positions = regular, regular, regular, 0, 1, 0, 1, 0, 1\n");
  fprintf(genfile, "\n");
  fprintf(genfile, "end\n");
  fclose(genfile);

  FILE* outfile = fopen(outfname, "w");
  for(x = 0; x<LS; x++)
	for(y = 0; y<LS; y++)
		for(z = 0; z<LS; z++)
		fprintf(outfile, "%4.6E\n",allmas[sidx][x][y][z]);
  fclose(outfile);

 };
 
 strcpy(cmd, "zip -r -j ");
 strcat(cmd, "kino");
 strcat(cmd, " ./jslices");

 system(cmd);
// system("rm -r -f ./jslices");
 
 return EXIT_SUCCESS;
};
