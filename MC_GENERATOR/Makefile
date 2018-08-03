#include Makefile.var

CC	=	./gcc34
CFLAGS  =	-O2 -c -Wall -I ./lattice/  -std=c99
LDFLAGS =       -O2 -std=c99 -Wall -I./lattice/ -L./src -lblas -llapack -lm
SOURCES =
	
OBJECTS         =       $(SOURCES:.c=.o)

.c.o:
	        $(CC) $(CFLAGS) $< -o $@


meanplaq: $(OBJECTS) test2.c
	$(CC) $(LDFLAGS) test2.c src/liblattice.a $(OBJECTS) -o $@
	
alignment: $(OBJECTS) alignment.c
	$(CC) $(LDFLAGS) alignment.c src/liblattice.a $(OBJECTS) -o $@        
	
mc: $(OBJECTS) mc.c
	$(CC) $(LDFLAGS) mc.c src/liblattice.a $(OBJECTS) -o $@        
	