#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Makefile for rslib/mw/sensor_model
#$Name$
#$Id$
#Copyright AER, Inc., 2002, 2003. All rights Reserved.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
include ../../makefile.common

# CVS_TAG is to be filled with a delivery tag at the release
CVS_TAG = $Name$
TB_STATUS = $(findstring release, $(CVS_TAG))

vpath %.f90 ${RSLIB_HOME}/util/misc
vpath %.f ${RSLIB_HOME}/util/misc

LIBPWD := ..
LIB := $(LIBPWD)/bin/${PLATFORM}${DBG}/libinst.a
GETLATLON := ${LIBPWD}/bin/${PLATFORM}${DBG}/getLatLon.mexglx
TBSCAN2SDR := ${LIBPWD}/bin/${PLATFORM}${DBG}/tbScan2sdr.mexglx
MATLAB2SCENE := ${LIBPWD}/bin/${PLATFORM}${DBG}/matlab2scene.mexglx
IDLSCENE := $(LIBPWD)/bin/${PLATFORM}${DBG}/idl2scene.so

ifeq ($(TB_STATUS),)
SRCS =  \
	sfcgrdDBtype_module.f90	\
	bilintSintype_module.f90 \
	sfcgrd_DB_module.f90 \
	AccessSfcGrid.f90 \
	NWPprof.f90 \
	getNWPprof.f90 \
	getTopogr.f90 \
	gridProf.f90 \
	gtopo.f \
	ChannelSelectionModule.f90 \
	GeomSat.f90 \
	setEmisQuasiPol.f90 \
	readMM5Terr.f \
	InterpNWPtools.f90 \
	MeasErrStructure.f90

LIBOBJS = \
	sfcgrdDBtype_module.o	\
	bilintSintype_module.o \
	sfcgrd_DB_module.o \
	AccessSfcGrid.o \
	NWPprof.o \
	getNWPprof.o \
	getTopogr.o \
	gridProf.o \
	gtopo.o \
	ChannelSelectionModule.o \
	GeomSat.o \
	setEmisQuasiPol.o \
	readMM5Terr.o \
	InterpNWPtools.o \
	MeasErrStructure.o
else
SRCS =  \
	AccessSfcGrid.f90 \
	gridProf.f90 \
	InterpNWPtools.f90

LIBOBJS = \
	AccessSfcGrid.o \
	gridProf.o \
	InterpNWPtools.o
endif

#############################
####### build targets #######
#############################

LIB: clean $(LIB)
ifneq ($(MACHINE),ppc64)
APP: ${IDLSCENE} $(GETLATLON) $(TBSCAN2SDR) $(MATLAB2SCENE)
else
APP: ${IDLSCENE}
endif

$(LIB): $(LIBOBJS)
	test  -d ${LIBPWD}/bin/${PLATFORM}${DBG} || mkdir ${LIBPWD}/bin/${PLATFORM}${DBG}
	ar -r $@ $?
	test  -d ${LIBPWD}/include/${PLATFORM}${DBG} || mkdir ${LIBPWD}/include/${PLATFORM}${DBG}
	find . -name '*.mod' -exec cp {} ${LIBPWD}/include/${PLATFORM}${DBG}/. \;

${GETLATLON}: getLatLon.o
	@#create the command directory if it does not exist
	test  -d ${LIBPWD}/bin/${PLATFORM}${DBG} || mkdir ${LIBPWD}/bin/${PLATFORM}${DBG}
	${LDR} ${FFLAGS} ${MFLAGS} ${RS_INCLUDES} ${MATLABOBJS} getLatLon.o \
	${RS_LIBS} -o ${GETLATLON}

$(TBSCAN2SDR): tbScan2sdr.o
	@#create the command directory if it does not exist
	test  -d ${LIBPWD}/bin/${PLATFORM}${DBG} || mkdir ${LIBPWD}/bin/${PLATFORM}${DBG}
	${LDR} ${FFLAGS} ${MFLAGS} ${RS_INCLUDES} ${MATLABOBJS} tbScan2sdr.o \
	${RS_LIBS} -o $(TBSCAN2SDR)

$(MATLAB2SCENE): matlab2scene.o
	@#create the command directory if it does not exist
	test  -d ${LIBPWD}/bin/${PLATFORM}${DBG} || mkdir ${LIBPWD}/bin/${PLATFORM}${DBG}
	${LDR} ${FFLAGS} ${MFLAGS} ${RS_INCLUDES} ${MATLABOBJS} matlab2scene.o \
	${RS_LIBS} -o $(MATLAB2SCENE)

IDL2SCENEOBJS := \
	idl2scene.o \
	${RSLIB_HOME}/util/misc/prof2CLW.o

${IDLSCENE}:  ${IDL2SCENEOBJS} \
	${RSLIB_HOME}/util/bin/${PLATFORM}${DBG}/*.a
	${LDR} ${FFLAGS} ${RS_INCLUDES} ${shared} ${IDL2SCENEOBJS} \
	${RS_LIBS} -o ${IDLSCENE}

$(fLIBOBJS): $(fLIBSRC)
$(FLIBOBJS): $(FLIBSRC)
$(f90LIBOBJS): $(f90LIBSRC)

MODULES := $(wildcard *.mod)
ifeq ($(strip $(MODULES)),)
MOD :=
RMMOD :=
else
MOD := ${LIBPWD}/include/${PLATFORM}${DBG}/${MODULES};
RMMOD := (cd ${LIBPWD}/include/${PLATFORM}${DBG}; /bin/rm -f ${MODULES};);
endif

clean:
	(/bin/rm -f *.o ${MODULES} ${LIB} $(GETLATLON) $(TBSCAN2SDR) $(MATLAB2SCENE) $(IDLSCENE));
	${RMMOD}
