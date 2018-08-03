/*{{{*/
/*!
 * \file read_ev.c
 *
 * \brief Calculate spectrum of Free Wilson Fermions
 *
 *
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \date   Втр Сен 28 01:10:14 MSD 2004
 */
/*}}}*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include "defs.h"
#include "types.h"
#include "mt19937ar.h"
#include "lattice.h"
#include "arpack.h"
#include "dirac.h"
#include "panic.h"
#include "linal.h"

#ifdef SU2
#include "su2math.h"
#endif

#ifdef SU3
#include "su3math.h"
#endif

void printhelp(void) __attribute__ ((noreturn));

void printhelp(void)
{
	fprintf(stdout, "Reader of EV files\n");
	fprintf(stdout, "Usage:\t-i evfilename -o outfilename [-r] \n");
	fprintf(stdout, "      \t-r -- orthogonality chechk\n");
	exit(EXIT_SUCCESS);
}

/*{{{*/
/*!
 * \brief  main routine
 *
 * \param argc
 * \param argv[]
 *
 * \return int 
 */
/*}}}*/
int main (int argc, char * argv[]) /*{{{*/
{
	int c;
	char *iname = NULL;
	char *oname = NULL;
	t_op_ev evd;
	FILE * evdata;
	int i,j;
	int res;
	int vol;
	int ortho = 0;
	
	while ((c = getopt(argc, argv, "i:o:r")) != -1) /*{{{*/
	{
		switch (c) {
		case 'i':
			iname = malloc(sizeof(char)*(strlen(optarg)+1));
			memset(iname, 0, sizeof(char)*(strlen(optarg)+1));
			strcpy(iname, optarg);
			break;
		case 'o':
			oname = malloc(sizeof(char)*(strlen(optarg)+1));
			memset(oname, 0, sizeof(char)*(strlen(optarg)+1));
			strcpy(oname, optarg);
			break;
		case 'r':
			ortho = 1;
			break;
		case '?':
			if (isprint(optopt))
				fprintf(stderr, "Unknown option '-%c'.\n", optopt);
			else
				fprintf(stderr,
						"Unknown option character '\\x%x'.\n", optopt);
			printhelp();
		default:
			printhelp();
		}
	}/*}}}*/
	if (!iname || !oname)
		printhelp();
	fprintf(stdout,"Input file :\t%s\n", iname);
	fprintf(stdout,"Output file:\t%s\n", oname);
	evd.op = NULL;
	evd.nev = 0;
	evd.which[0] = 0;
	evd.which[1] = 0;
	evd.val_vec = 0;
	evd.tol = 0;
	evd.evals = NULL;
	evd.evecs = NULL;
	evd.data = 0;
	evdata = fopen(oname, "wb");

	res = read_ev(iname, &evd, SM_SECTION);

		for (vol = 1, i = 0; i < evd.dim-1; i++)
			vol *= evd.ls;
		vol *= evd.lt*evd.ndirac*evd.ncolors;


	if (!((res & E_READ_CANTREAD) || (res & E_READ_CANTOPEN) || (res & E_READ_SECNOTHERE))) {

		for (i = 0; i < evd.nev; i++){
			fprintf(evdata, "%12.10f\t%12.10f\n", CCALL(creal)(evd.evals[i]),
					CCALL(cimag)(evd.evals[i]));
			for (j = 0; j < vol; j++) 
			fprintf(evdata, "%12.10f\t%12.10f\n",CCALL(creal)(*(evd.evecs + i*vol + j)),
				CCALL(cimag)(*(evd.evecs + i*vol + j)));
			

		}
	}
	if (ortho) {

		check_ortho(evd.nev, vol, evd.evals, evd.evecs);
	}
	if (evd.evals != NULL) {
		free(evd.evals);
		evd.evals = NULL;
	}
	if (evd.evecs != NULL) {
		free(evd.evecs);
		evd.evecs = NULL;
	}

	fclose(evdata);

	return EXIT_SUCCESS;
}/*}}}*/
