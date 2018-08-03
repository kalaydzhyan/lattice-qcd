# quick hack over BLAS library. Routines not needed for ARPACK are removed.

include		../../Makefile.var

LAPACK_SRC	=	cgeqr2.f clacpy.f clahqr.f clanhs.f clarnv.f clartg.f clascl.f \
				claset.f ctrevc.f ctrsen.f cunm2r.f dgeqr2.f dgeqrf.f dlabad.f \
				dlamch.f dlapy2.f dlarfb.f dlarf.f dlarfg.f dlarft.f dorm2r.f  \
				dormqr.f dtrtrs.f ieeeck.f ilaenv.f slabad.f slamch.f slapy2.f \
				xerbla.f zgeqr2.f zlacpy.f zlahqr.f zlanhs.f zlarnv.f zlartg.f \
				zlascl.f zlaset.f ztrevc.f ztrsen.f zunm2r.f lsame.f clarf.f   \
				clarfg.f zlarf.f zlarfg.f cladiv.f zladiv.f slapy3.f dlapy3.f  \
				classq.f zlassq.f clange.f zlange.f sladiv.f dladiv.f slaruv.f \
				dlaruv.f clatrs.f zlatrs.f ctrexc.f ztrexc.f clacon.f ctrsyl.f \
				zlacon.f ztrsyl.f crot.f zrot.f icmax1.f izmax1.f scsum1.f     \
				dzsum1.f
				
LAPACK_OBJ	=	$(LAPACK_SRC:.f=.o)

###############################################################################

.SUFFIXES:
.SUFFIXES: .f .o

all: liblapack

clean: 
	rm -f *.o a.out core

liblapack: $(LAPACK_OBJ)
	@$(AR) $(AROPTS) $(LAPACK_LIB) $(LAPACK_OBJ)
	@$(RANLIB) $(LAPACK_LIB)
	@mv $(LAPACK_LIB) ../


.o:
	$(FC) $(FFLAGS) -c $< -o $@

dlamch.o: dlamch.f
	$(FC) $(NOOPT) -c $< -o $@

slamch.o: slamch.f
	$(FC) $(NOOPT) -c $< -o $@
