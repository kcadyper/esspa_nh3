#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# $Id$ 
# Makefile for rslib/ir/rtm/
# Copyright AER, Inc., 2002, 2003. All rights Reserved.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
include ../../makefile.common

LIBPWD := ..
LIB := ${LIBPWD}/bin/${PLATFORM}${DBG}/libirrtm.a 

f90SRC := \
	FiniteDiffModule.f90 \
	IRRTmodule.f90 \
	ImgRTmodule.f90
fLIBOBJS := $(patsubst %.f,%.o, $(wildcard *.f))
FLIBOBJS := $(patsubst %.F,%.o, $(wildcard *.F))
#f90LIBOBJS := $(patsubst %.f90,%.o, $(wildcard *.f90))
f90LIBOBJS := $(patsubst %.f90,%.o, $(f90SRC))

LIBOBJS := $(fLIBOBJS) $(FLIBOBJS) $(f90LIBOBJS)

OSSIRPATH:=$(OSSPATH)/src/ir
OSSLIB:=$(OSSPATH)/lib/$(PLATFORM)/libfossir.a

#############################
####### build targets #######
#############################

ALL: clean $(OSSREF) ossready $(LIB)

$(LIB): $(LIBOBJS)
	test  -d ${LIBPWD}/bin/${PLATFORM}${DBG} || mkdir ${LIBPWD}/bin/${PLATFORM}${DBG}
	cp -f $(OSSLIB) $(LIB)
	ar -r $@ $? 
	test  -d ${LIBPWD}/include/${PLATFORM}${DBG} || mkdir ${LIBPWD}/include/${PLATFORM}${DBG}
	find . -name '*.mod' -exec cp {} ${LIBPWD}/include/${PLATFORM}${DBG}/. \;

$(OSSREF):
	@echo MSG: MISSING ossrtm - MAY NEED git submodule init AND git submodule update, OR USE --recursive WHEN CLONING
	test -e $(OSSREF)

ossready:
	@cd $(OSSIRPATH); $(MAKE) flib

MODULES := $(wildcard *.mod)
ifeq ($(strip $(MODULES)),)
MOD :=
RMMOD :=
else
MOD := ${LIBPWD}/include/${PLATFORM}${DBG}/${MODULES};
RMMOD := (cd ${LIBPWD}/include/${PLATFORM}${DBG}; /bin/rm -f ${MODULES};);
endif

clean:
	@echo "cleaning..." $(LIBOBJS) ${MODULES} ${MOD}
	(/bin/rm -f $(LIBOBJS) ${MODULES} ${LIB});
	${RMMOD}
