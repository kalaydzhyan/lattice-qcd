# quick hack over BLAS library. Routines not needed for ARPACK are removed.

include		../../Makefile.var

BLAS_SRC	=	caxpy.f	ccopy.f	cdotc.f	cgemv.f	cgeru.f cscal.f csscal.f 	\
				cswap.f ctrmm.f dcopy.f dgemm.f dgemv.f dger.f dnrm2.f		\
				dscal.f	dtrmm.f	dtrmv.f	dtrsm.f	dznrm2.f scnrm2.f zaxpy.f	\
				zcopy.f	zdotc.f zdscal.f zgemv.f zgeru.f zscal.f zswap.f 	\
				ztrmm.f icamax.f izamax.f cgerc.f zgerc.f
BLAS_SRC	=	caxpy.f ccopy.f cdotc.f cdotu.f cgbmv.f cgemm.f cgemv.f cgerc.f\
				cgeru.f cscal.f csscal.f cswap.f ctbsv.f ctrmm.f ctrsm.f       \
				ctrsv.f dasum.f daxpy.f dcopy.f ddot.f dgbmv.f dgemm.f dgemv.f \
				dger.f dnrm2.f drot.f drotg.f dscal.f dswap.f dsymv.f dsyr2.f  \
				dtbsv.f dtrmm.f dtrmv.f dtrsm.f dzasum.f dznrm2.f icamax.f     \
				idamax.f isamax.f izamax.f sasum.f saxpy.f scasum.f scnrm2.f   \
				scopy.f sdot.f sgbmv.f sgemm.f sgemv.f sger.f snrm2.f srot.f   \
				srotg.f sscal.f sswap.f ssymv.f ssyr2.f stbsv.f strmm.f strsm.f\
				zaxpy.f zcopy.f zdotc.f zdotu.f zdscal.f zgbmv.f zgemm.f       \
				zgemv.f zgerc.f zgeru.f zscal.f zswap.f ztbsv.f ztrmm.f ztrsm.f\
				ztrsv.f

BLAS_OBJ	=	$(BLAS_SRC:.f=.o)

###############################################################################

.SUFFIXES:
.SUFFIXES: .f .o

all: libblas

clean: 
	rm -f *.o a.out core

libblas: $(BLAS_OBJ)
	@$(AR) $(AROPTS) $(BLAS_LIB) $(BLAS_OBJ)
	@$(RANLIB) $(BLAS_LIB)
	@mv $(BLAS_LIB) ../


.o:
	$(FC) $(FFLAGS) -c $< -o $@

