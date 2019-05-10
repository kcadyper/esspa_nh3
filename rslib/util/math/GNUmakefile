#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#makefile for rs/util/math
#$Name$ 
#$Id$ 
#Copyright AER, Inc., 2002, 2003. All rights Reserved.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include ../../makefile.common

# CVS_TAG is to be filled with a delivery tag at the release
CVS_TAG = $Name$
TB_STATUS = $(findstring release, $(CVS_TAG))

LIBPWD := ..
LIB := ${LIBPWD}/bin/${PLATFORM}${DBG}/libmath.a 

ifeq ($(TB_STATUS),)
LIBSRC :=  \
	RandomMT.f90 \
	gauss.f90 \
	DataRegres.f90 \
	Invert.f90 \
	Quadrature.f90 \
	fft.f \
	recursiveFilter.f \
	compress_gen.f

LIBOBJ := \
	RandomMT.o \
	gauss.o \
	DataRegres.o \
	Invert.o \
	Quadrature.o \
	fft.o \
	recursiveFilter.o \
	compress_gen.o
else
LIBSRC :=  \
	RandomMT.f90 \
	gauss.f90 \
	Invert.f90

LIBOBJ := \
	RandomMT.o \
	gauss.o \
	Invert.o
endif

#############################
####### build targets #######
#############################

All: clean $(LIB)

$(LIB): $(LIBOBJ)
	test -d ${LIBPWD}/bin/${PLATFORM}${DBG} || mkdir ${LIBPWD}/bin/${PLATFORM}${DBG}
	ar -r $@ $?
	test  -d ${LIBPWD}/include/${PLATFORM}${DBG} || mkdir ${LIBPWD}/include/${PLATFORM}${DBG}
	find . -name '*.mod' -exec cp {} ${LIBPWD}/include/${PLATFORM}${DBG}/. \;

$(LIBOBJ): $(LIBSRC)

MODULES := $(wildcard *.mod)
ifeq ($(strip $(MODULES)),)
MOD :=
RMMOD :=
else
MOD := ${LIBPWD}/include/${PLATFORM}${DBG}/${MODULES};
RMMOD := (cd ${LIBPWD}/include/${PLATFORM}${DBG}; /bin/rm -f ${MODULES};);
endif

clean:
#	@echo "cleaning..." $(LIBOBJ) ${MODULES} ${MOD} ${LIB}
	(/bin/rm -f $(LIBOBJ) ${MODULES} ${LIB});
	${RMMOD}
