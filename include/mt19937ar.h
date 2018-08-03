#ifndef _MT19937AR_H_
#define _MT19937AR_H_

#define MT_SEED_LEN	100

extern unsigned long mt[];
void init_genrand(unsigned long s);
void init_by_array(unsigned long init_key[], unsigned long key_length);
long genrand_int31(void);
unsigned long genrand_int32(void);
t_real genrand_real1(void);
t_real genrand_real2(void);
t_real genrand_real3(void);
double genrand_res53(void);

#endif
