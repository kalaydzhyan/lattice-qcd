#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

extern void   rnd_init();
extern double rnd_get();
extern void rnd_seed_set();
#endif
