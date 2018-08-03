# VERBOSE
# DEBUG
###########################################################################
CFLAGS      = -Wall  -I ../lattice -I ~/.lib/fftw/include
CFLAGS     += -DVERBOSE
#CFLAGS     += -g -pg -DDEBUG
CFLAGS    += -O2

DEST	      = ../.

EXTHDRS	      =

HDRS	      =

INSTALL	      = cp

LD	      = $(CXX)

LDFLAGS	      = -march=i386 -O2

AR            = /usr/bin/ar

ARFLAGS       = cr

RANLIB        = /usr/bin/ranlib

RANLIBFLAGS   =

LIBS	      =

MAKEFILE      = Makefile

OBJS	      = MC-SU2.o \
		MC-SU_N.o \
		MC-U1.o \
		MC-align.o \
		MC-killing.o \
		SU2-utils.o \
		SU_N-utils.o \
		bessel.o \
		diagonalize.o \
		dedekind.o \
		geometry.o \
		lapack.o \
		names.o \
		newton-raphson.o \
		random.o \
		spinor.o \
		svd.o \
		triangulation.o \
		wilson-loops.o

PRINT	      = lpr

PROGRAM       = liblattice.a

SHELL	      = /bin/sh

SRCS	      = MC-SU2.c \
		MC-SU_N.c \
		MC-U1.c \
		MC-align.c \
		MC-killing.c \
		SU2-utils.c \
		SU_N-utils.c \
		bessel.c \
		diagonalize.c \
		dedekind.c \
		geometry.c \
		lapack.c \
		names.c \
		newton-raphson.c \
		random.c \
		spinor.c \
		svd.c \
		triangulation.c \
		wilson-loops.c \

SYSHDRS	      =

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS)
		@echo "Linking $(PROGRAM) ..."
		@echo $(AR) $(ARFLAGS) $(PROGRAM) $(OBJS)
		@$(AR) $(ARFLAGS) $(PROGRAM) $(OBJS)
		@echo $(RANLIB) $(RANLIBFLAGS) $(PROGRAM)
		@$(RANLIB) $(RANLIBFLAGS) $(PROGRAM)
		@echo "done"

clean:;		@rm -f $(OBJS) core

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags

depend:;	@mkmf -f $(MAKEFILE)

echo:;		@echo $(HDRS) $(SRCS)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) $(PROGRAM) $(DEST)); fi

print:;		@$(PRINT) $(HDRS) $(SRCS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)
