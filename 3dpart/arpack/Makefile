# quick hack over BLAS library. Routines not needed for ARPACK are removed.

include		../../Makefile.var

ARPACK_SRC	=	cgetv0.f cmout.f cnaitr.f cnapps.f cnaup2.f cnaupd.f cneigh.f \
				cneupd.f cngets.f csortc.f cstatn.f cvout.f debug.h dvout.f   \
				ivout.f second.f zgetv0.f zmout.f znaitr.f znapps.f znaup2.f  \
				znaupd.f zneigh.f zneupd.f zngets.f zsortc.f zstatn.f zvout.f \
				svout.f

ARPACK_OBJ	=	$(ARPACK_SRC:.f=.o)

###############################################################################

.SUFFIXES:
.SUFFIXES: .f .o

all: libarpack

clean: 
	rm -f *.o a.out core

libarpack: $(ARPACK_OBJ)
	@$(AR) $(AROPTS) $(ARPACK_LIB) $(ARPACK_OBJ)
	@$(RANLIB) $(ARPACK_LIB)
	@mv $(ARPACK_LIB) ../


.o:
	$(FC) $(FFLAGS) -c $< -o $@

