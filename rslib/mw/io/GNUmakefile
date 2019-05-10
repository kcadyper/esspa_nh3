#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#makefile for rslib/mw/io
#$Name$ 
#$Id$ 
#Copyright AER, Inc., 2002, 2003. All rights Reserved.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include ../../makefile.common

LIBPWD := ..
LIB := ${LIBPWD}/bin/${PLATFORM}${DBG}/libmwio.a 
IDLSCFREAD := ${LIBPWD}/bin/${PLATFORM}${DBG}/SCFread.so

fLIBOBJS := $(patsubst %.f,%.o, $(wildcard *.f))
FLIBOBJS := $(patsubst %.F,%.o, $(wildcard *.F))
f90LIBOBJS := $(patsubst %.f90,%.o, $(wildcard *.f90 rtm/*.f90))

LIBOBJS := $(fLIBOBJS) $(FLIBOBJS) $(f90LIBOBJS)

#############################
####### build targets #######
#############################

LIB: clean $(LIB)
IDL: ${IDLSCFREAD}

$(LIB): $(LIBOBJS)
	test  -d ${LIBPWD}/bin/${PLATFORM}${DBG} || mkdir ${LIBPWD}/bin/${PLATFORM}${DBG}
	ar -r $@ $? 
	test  -d ${LIBPWD}/include/${PLATFORM}${DBG} || mkdir ${LIBPWD}/include/${PLATFORM}${DBG}
	find . -name '*.mod' -exec cp {} ${LIBPWD}/include/${PLATFORM}${DBG}/. \;

${IDLSCFREAD}: SCFread.o ${RSLIB_HOME}/util/bin/${PLATFORM}${DBG}/*.a
	${LDR} ${FFLAGS} ${RS_INCLUDES} ${shared} SCFread.o ${RS_LIBS} -o ${IDLSCFREAD}

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
	(/bin/rm -f $(LIBOBJS) ${MODULES} ${LIB} $(IDLSCFREAD));
	${RMMOD}





