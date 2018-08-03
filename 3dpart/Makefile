# Makefile for 3dpart libraries

include		 ../Makefile.var

LIB_DIR		=	../lib

SUBDIRS		=	blas lapack arpack

all: dirs 

clean:
	rm -f *.o a.out core
	@list='$(SUBDIRS)'; for subdir in $$list; do \
		( cd $$subdir && $(MAKE) clean); \
	done

dirs:
	@list='$(SUBDIRS)'; for subdir in $$list; do \
		( cd $$subdir && $(MAKE) all); \
	done
	@mv $(BLAS_LIB) $(LIB_DIR)
	@mv $(LAPACK_LIB) $(LIB_DIR)
	@mv $(ARPACK_LIB) $(LIB_DIR)

