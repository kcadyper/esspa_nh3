#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Makefile for rslib/mw/core
#$Name$ 
#$Id$ 
#Copyright AER, Inc., 2002, 2003. All rights Reserved.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
include ../../makefile.common

LIBPWD := ..
LIB := $(LIBPWD)/bin/${PLATFORM}${DBG}/libretr.a 

ifeq (${BUILD_VERSION},AtmClass_BtDiff)
    ATMOS_CLASS_METHOD =_BtDiff
else
    ATMOS_CLASS_METHOD =
endif

f90SRC := \
	ChkValid.f90 \
	LimitsModule.f90 \
	LinInvert.f90 \
	MapInvert.f90 \
	setxPtrR.f90 \
	AtmosClassModule${ATMOS_CLASS_METHOD}.f90 \
	BackgroundModule.f90
fLIBOBJS :=$(patsubst %.f,%.o, $(wildcard *.f))
FLIBOBJS :=$(patsubst %.F,%.o, $(wildcard *.F))
f90LIBOBJS :=$(patsubst %.f90,%.o, $(f90SRC))

LIBOBJS:=$(fLIBOBJS) $(FLIBOBJS) $(f90LIBOBJS)

#############################
####### build targets #######
#############################

ALL: clean $(LIB)

$(LIB): $(LIBOBJS)
	test  -d ${LIBPWD}/bin/${PLATFORM}${DBG} || mkdir ${LIBPWD}/bin/${PLATFORM}${DBG}
	ar -r $@ $? 
	test  -d ${LIBPWD}/include/${PLATFORM}${DBG} || mkdir ${LIBPWD}/include/${PLATFORM}${DBG}
	find . -name '*.mod' -exec cp {} ${LIBPWD}/include/${PLATFORM}${DBG}/. \;

MODULES := $(wildcard *.mod)
ifeq ($(strip $(MODULES)),)
MOD :=
RMMOD :=
else
MOD := ${LIBPWD}/include/${PLATFORM}${DBG}/${MODULES};
RMMOD := (cd ${LIBPWD}/include/${PLATFORM}${DBG}; /bin/rm -f ${MODULES};);
endif

clean:
#	@echo "cleaning..." $(LIBOBJS) ${MODULES} ${MOD} ${LIB}
	(/bin/rm -f $(LIBOBJS) ${MODULES} ${LIB});
	${RMMOD}
