#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Makefile for rslib/ir/core
#$Name$
#$Id$
#Copyright AER, Inc., 2002, 2003. All rights Reserved.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
include ../../makefile.common

LIBPWD := ..
LIB := $(LIBPWD)/bin/${PLATFORM}${DBG}/libircore.a

f90SRC := \
	IRBkgPrepModule.f90 \
	IRMapInvert.f90 \
	RegrPrepModule.f90 \
	CloudMitigateModule.f90 \
	StatRetrModule.f90 \
	SecondaryRetrModule.f90 \
	PrimaryRetrModule.f90
fLIBOBJS :=$(patsubst %.f,%.o, $(wildcard *.f))
FLIBOBJS :=$(patsubst %.F,%.o, $(wildcard *.F))
#f90LIBOBJS :=$(patsubst %.f90,%.o, $(wildcard *.f90))
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
