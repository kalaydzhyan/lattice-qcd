###############################################################################
# Makefile for the  iSPECdo ToolKit
###############################################################################

include			Makefile.var

PROJECT		=	iSPECdo
VERSION		=	0.25

HEADERS		=	defs.h \
			types.h \
			mt19937ar.h \
			arnoldi.h \
			lattice.h \
			statistics.h

ifeq ("$(NCOLORS)", "2")
HEADERS		+=	su2math.h
endif

ifeq ("$(NCOLORS)", "3")
HEADERS		+=	su3math.h
endif

SRC_C		=	lattice.c \
			dirac.c	\
			polacc.c \
			linal.c 
			
				
OBJ_NM		=	mt19937ar.o \
			panic.o	\
			arpack.o \
			arn_log.o \
			minmax.o \
			timing.o \
			statistics.o 

ifeq ("$(NCOLORS)", "2")
OBJ_NM		+=	su2math.o
endif

ifeq ("$(NCOLORS)", "3")
OBJ_NM		+=	su3math.o
endif

				
SRC_F		=	

OBJS		=	$(addprefix $(OBJ_DIR)/, $(notdir $(OBJ_NM)))
OBJS		+=	$(addprefix $(OBJ_DIR)/, $(notdir $(SRC_C:.c=$(SUFF).o)))
OBJS		+=	$(addprefix $(OBJ_DIR)/, $(notdir $(SRC_F:.f=$(SUFF).o)))
			
SUBDIRS		=

###############################################################################
# targets
###############################################################################

all: wd ovd prop meson

test:
	echo $(OBJS)

clean:
	rm -f -v *.o
	rm -f -v $(DOXYCONF)
	@( \
		cd $(OBJ_DIR); \
		rm -f -v *.o; \
		cd ..; \
	)

# TODO better cleaning
clobber:
	rm -f *.o
	@( \
		cd $(OBJ_DIR); \
		rm -f -v *.o; \
		cd ..; \
		cd $(BIN_DIR); \
		rm -v -i *; \
		cd ..; \
	)
	rm -rif ./doc/html ./doc/latex

wd:     $(OBJS)
	$(LD) -static -std=c99 -o $(BIN_DIR)/$@$(SUFF) $(SRC_DIR)/wd.c $(OBJS) $(LDFLAGS) $(CFLAGS)

ovd:    $(OBJS)
	$(LD) -static -std=c99 -o $(BIN_DIR)/$@$(SUFF) $(SRC_DIR)/ovd.c $(OBJS) $(LDFLAGS) $(CFLAGS)
	
prop:   $(OBJS)
	$(LD) -static -std=c99 -o $(BIN_DIR)/$@$(SUFF) $(SRC_DIR)/propagators.c $(OBJS) $(LDFLAGS) $(CFLAGS) 
	
meson:  $(OBJS)
	$(LD) -static -std=c99 -o $(BIN_DIR)/$@$(SUFF) $(SRC_DIR)/meson.c $(OBJS) $(LDFLAGS) $(CFLAGS) 



read_ev: $(OBJS)
	$(LD) -static -std=c99 -o $(BIN_DIR)/$@$(SUFF) $(SRC_DIR)/read_ev.c $(OBJS) $(LDFLAGS) $(CFLAGS)

mssc:	$(OBJS)
	$(LD) -static -std=c99 -o $(BIN_DIR)/$@$(SUFF) $(SRC_DIR)/mssc.c $(OBJS) $(LDFLAGS) $(CFLAGS)
	
curr:	$(OBJS)
	$(LD) -static -std=c99 -o $(BIN_DIR)/$@$(SUFF) $(SRC_DIR)/currents.c $(OBJS) $(LDFLAGS) $(CFLAGS)
mcurr:	$(OBJS)
	$(LD) -static -std=c99 -o $(BIN_DIR)/$@$(SUFF) $(SRC_DIR)/mcurrents.c $(OBJS) $(LDFLAGS) $(CFLAGS)
jxjy:	$(OBJS)
	$(LD) -static -std=c99 -o $(BIN_DIR)/$@$(SUFF) $(SRC_DIR)/jxjy.c $(OBJS) $(LDFLAGS) $(CFLAGS)
cr52j2:	$(OBJS)
	$(LD) -static -std=c99 -o $(BIN_DIR)/$@$(SUFF) $(SRC_DIR)/j2r2.c $(OBJS) $(LDFLAGS) $(CFLAGS)
tpc:	$(OBJS)
	$(LD) -static -std=c99 -o $(BIN_DIR)/$@$(SUFF) $(SRC_DIR)/tpc.c $(OBJS) $(LDFLAGS) $(CFLAGS)

mpl:	$(OBJS)
	$(LD) -static -std=c99 -o $(BIN_DIR)/$@$(SUFF) $(SRC_DIR)/mpl.c $(OBJS) $(LDFLAGS) $(CFLAGS)

artest:	$(OBJS)
	$(LD) -static -std=c99 -o $(BIN_DIR)/$@ $(SRC_DIR)/artest.c $(OBJS) $(LDFLAGS) $(CFLAGS)
	
shumr:  $(OBJS)
	$(LD) -static -std=c99 -o $(BIN_DIR)/$@ $(SRC_DIR)/shumr_test.c $(OBJS) $(LDFLAGS) $(CFLAGS)

ftest:	$(OBJS)
	$(LD) -static -std=c99 -o $(BIN_DIR)/$@ $(SRC_DIR)/ftest.c $(OBJS) $(LDFLAGS) $(CFLAGS)


dirs:
	@list='$(SUBDIRS)'; for subdir in $$list; do \
		( cd $$subdir && $(MAKE)  all); \
	done

.SUFFIXES:
.SUFFIXES: .c .o .f .h

$(OBJ_DIR)/%$(SUFF).o: $(SRC_DIR)/%.c
	$(CC) -c -o $@ $(CFLAGS) $<


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) -c -o $@ $(CFLAGS) $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) -c -o $@ $(FFLAGS) $<
	
%.h:	
	@echo $@

# code documentation parameters

DOXYGEN		= doxygen
DOXYCONF 	= doxygen.conf
DOCDIR 		= ./doc
DOCLANG 	= English

# doxygen configuration file

$(DOXYCONF):
	@echo -e "			\
PROJECT_NAME           = $(PROJECT)	\n\
PROJECT_NUMBER         = $(VERSION)	\n\
OUTPUT_DIRECTORY       = $(DOCDIR)	\n\
OUTPUT_LANGUAGE        = $(DOCLANG)	\n\
EXTRACT_ALL            = NO		\n\
SHOW_INCLUDE_FILES     = NO		\n\
OPTIMIZE_OUTPUT_FOR_C  = YES		\n\
INPUT                  = $(SRC_DIR) $(INC_DIR)		\n\
FILE_PATTERNS          = *.c *.h	\n\
RECURSIVE              = NO		\n\
IMAGE_PATH             =		\n\
SOURCE_BROWSER         = YES		\n\
HTML_HEADER            = 		\n\
HTML_FOOTER            = 		\n\
HTML_STYLESHEET        = 		\n\
GENERATE_TREEVIEW      = YES		\n\
COMPACT_LATEX          = YES		\n\
EXTRA_PACKAGES         = 		\n\
LATEX_HEADER           =		\n\
LATEX_HIDE_INDICES     = YES 		\n\
LATEX_HIDE_INDICES     = YES		\n\
IMAGE_PATH             = ./doc		\n\
TAB_SIZE               = 4"		\
			> $(DOXYCONF)
			
# doxygen target

doc: $(DOXYCONF) $(SRC_DIR)/*.c $(INC_DIR)/*.h
	$(DOXYGEN) $(DOXYCONF)
	cd $(DOCDIR)/latex && $(MAKE) ps
