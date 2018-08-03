#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
float parbeta=atof(argv[1]);
char fname[50], f2name[50];
sprintf(fname, "log%2.3f.txt", parbeta);
sprintf(f2name, "uzero%2.3f.txt", parbeta);
printf("\n%s\n", fname);
exit(EXIT_SUCCESS);
}