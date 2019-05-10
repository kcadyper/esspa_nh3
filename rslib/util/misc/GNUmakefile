#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Makefile for rslib/util/misc
#$Name$ 
#$Id$ 
#Copyright AER, Inc., 2002, 2003. All rights Reserved.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include ../../makefile.common

# CVS_TAG is to be filled with a delivery tag at the release
CVS_TAG = $Name$
TB_STATUS = $(findstring release, $(CVS_TAG))

LIBPWD := ..
LIB := ${LIBPWD}/bin/${PLATFORM}${DBG}/libmisc.a 
IDLIO := ${LIBPWD}/bin/${PLATFORM}${DBG}/idlio.so

ifeq ($(TB_STATUS),)
SRC = \
	VertCoord.f90 \
	StateIndexModule.f90 \
	LvlInterp.f90 \
	prof_int.f90 \
	date_module.f \
	angles.f \
	grav.f90 \
	planck.f \
	sunang.f \
	GeomSpheric.f90 \
        mapConvert.f90 \
	FOV2sfcgrd_module.f90 \
	SpectralOperations.f90 \
	MetFunctions.f90 \
	ToolboxModule.f90

IDLIOSRCS = \
	prof2CLW.f

LIBOBJS:= \
	VertCoord.o \
	StateIndexModule.o \
	LvlInterp.o \
	prof_int.o \
	date_module.o \
	angles.o \
	grav.o \
	planck.o \
	sunang.o \
	GeomSpheric.o \
        mapConvert.o \
	FOV2sfcgrd_module.o \
	SpectralOperations.o \
	MetFunctions.o \
	ToolboxModule.o

IDLIOOBJS = \
	prof2CLW.o
else
SRC = \
	VertCoord.f90 \
	StateIndexModule.f90 \
	LvlInterp.f90 \
	prof_int.f90 \
	date_module.f \
	grav.f90 \
	planck.f \
	GeomSpheric.f90 \
        mapConvert.f90 \
	SpectralOperations.f90 \
	MetFunctions.f90 \
	ToolboxModule.f90

LIBOBJS:= \
	VertCoord.o \
	StateIndexModule.o \
	LvlInterp.o \
	prof_int.o \
	date_module.o \
	grav.o \
	planck.o \
	GeomSpheric.o \
        mapConvert.o \
	SpectralOperations.o \
	MetFunctions.o \
	ToolboxModule.o
endif

#############################
####### build targets #######
#############################

LIB: clean $(LIB)
IDL: $(IDLIO)

$(LIB): $(LIBOBJS)
	test  -d ${LIBPWD}/bin/${PLATFORM}${DBG} || mkdir -p ${LIBPWD}/bin/${PLATFORM}${DBG}
	ar -r $@ $? 
	test  -d ${LIBPWD}/include/${PLATFORM}${DBG} || mkdir -p ${LIBPWD}/include/${PLATFORM}${DBG}
	find . -name '*.mod' -exec cp {} ${LIBPWD}/include/${PLATFORM}${DBG}/. \;

${IDLIO}: ${IDLIOOBJS} ${RSLIB_HOME}/util/bin/${PLATFORM}${DBG}/*.a
	${LDR} ${FFLAGS} ${RS_INCLUDES} ${shared} ${IDLIOOBJS} ${RS_LIBS} -o ${IDLIO}


MODULES := $(wildcard *.mod)
ifeq ($(strip $(MODULES)),)
MOD :=
RMMOD :=
else
MOD := ${LIBPWD}/include/${PLATFORM}${DBG}/${MODULES};
RMMOD := (cd ${LIBPWD}/include/${PLATFORM}${DBG}; /bin/rm -f ${MODULES};);
endif

clean:
	(/bin/rm -f *.o ${MODULES} ${LIB});
	${RMMOD}





