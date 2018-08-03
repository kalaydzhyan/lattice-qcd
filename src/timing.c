#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "timing.h"

long tim_ncalls[TIMING_MAX_ITEMS];
long long tim_timer[TIMING_MAX_ITEMS];

char tim_strs[TIMING_MAX_ITEMS][TIMING_STR_LEN] =
{
	"wdirac\t",
	"ov_dirac_mm",
	"h_wdirac",
	"h_wdirac_sq",
	"ov_dirac_ln",
	"h_wdirac_ln",
	"a_mul_cdsv",
	"ad_mul_cdsv"
};

#if (I686)

inline void start_time(int ti)
{
	__asm__ __volatile__(
			"pushl	%%eax\n\t"
			"pushl	%%edx\n\t"
			"rdtsc\n\t"
			"subl	%%eax, %0\n\t"
			"sbbl	%%edx, 4+%0\n\t"
			"addl	$1, %1\n\t"
			"popl	%%edx\n\t"
			"popl	%%eax"
			: "=m" (tim_timer[ti]), "=m" (tim_ncalls[ti])
			:
			: "%edx", "%eax");
}

inline void end_time(int ti)
{
	__asm__ __volatile__ (
			"pushl	%%eax\n\t"
			"pushl	%%edx\n\t"
			"rdtsc\n\t"
			"addl	%%eax, %0\n\t"
			"adcl	%%edx, 4+%0\n\t"
			"popl	%%edx\n\t"
			"popl	%%eax"
			: "=m" (tim_timer[ti])
			:
			: "%eax", "%edx");
}
#endif

/*****************************************************************************/
void init_time(void)
{
	int i;

	for (i = 0; i < TIMING_MAX_ITEMS; i++) {
		tim_ncalls[i] = 0;
		tim_timer[i] = 0;
	}
}

/*****************************************************************************/
void print_time(void)
{
	int i;

	for (i = 0; i < TIMING_MAX_ITEMS; i++)
		fprintf(stdout, "TIMING: %s\t%li\t%lli\t%f\n", &(tim_strs[i][0]),
				tim_ncalls[i], tim_timer[i], (tim_timer[i]*1.0/tim_ncalls[i]));
}
