#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Makefile for rslib/mw/rtm/
#$Name$ 
#$Id$ 
#Copyright AER, Inc., 2002, 2003. All rights Reserved.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
include ../../makefile.common

LIBPWD := ..
LIB := ${LIBPWD}/bin/${PLATFORM}${DBG}/libmwrtm.a 

f90SRC := \
	PlanckTabMW.f90 \
	MWRTmodule.f90
fLIBOBJS := $(patsubst %.f,%.o, $(wildcard *.f))
FLIBOBJS := $(patsubst %.F,%.o, $(wildcard *.F))
#f90LIBOBJS := $(patsubst %.f90,%.o, $(wildcard *.f90))
f90LIBOBJS := $(patsubst %.f90,%.o, $(f90SRC))

LIBOBJS := $(fLIBOBJS) $(FLIBOBJS) $(f90LIBOBJS)

OSSMWPATH:=$(OSSPATH)/src/mw
OSSLIB:=$(OSSPATH)/lib/$(PLATFORM)/libfossmw.a

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
	cp -f mxdims.incl ${LIBPWD}/include/${PLATFORM}${DBG}/. 

# spemiss.f is used solely for Monte Carlo model and it is written 
# in fortran 77 specific language, and CFLAGS instead of FFLAGS
spemiss.o: spemiss.f
	${FC77} -c $(CFLAGS) spemiss.f

rough_ks.o: rough_ks.f
	${FC77} -c $(CFLAGS) rough_ks.f

$(OSSREF):
	@echo MSG: MISSING ossrtm - MAY NEED git submodule init AND git submodule update, OR USE --recursive WHEN CLONING
	test -e $(OSSREF)

ossready:
	@cd $(OSSMWPATH); $(MAKE) flib

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
