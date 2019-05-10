#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Makefile for rslib/mw/sensor_model/
#$Name$ 
#$Id$ 
#Copyright AER, Inc., 2002, 2003. All rights Reserved.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
include ../../makefile.common

# CVS_TAG is to be filled with a delivery tag at the release
CVS_TAG = $Name$
TB_STATUS = $(findstring release, $(CVS_TAG))

vpath %.f90 ${RSLIB_HOME}/mw/io
vpath %.f ${RSLIB_HOME}/mw/io

LIBPWD := ..
LIB := ${LIBPWD}/bin/${PLATFORM}${DBG}/libNoise.a 
GETCONST := ${LIBPWD}/bin/${PLATFORM}${DBG}/getConst.mexglx
IDLMV := ${LIBPWD}/bin/${PLATFORM}${DBG}/idlmvltlnd.so
IDLCONST := ${LIBPWD}/bin/${PLATFORM}${DBG}/idlConstInit.so
IDLCONSTFEW := ${LIBPWD}/bin/${PLATFORM}${DBG}/idlConstNoSCF.so
IDLPRESTOALT := $(LIBPWD)/bin/${PLATFORM}${DBG}/idlPresToAlt.so

ifeq ($(TB_STATUS),)
LIBSRCS := \
	SensorData.f90 \
	NoiseModule.f90 \
	CalibSim.f90 \
	PointSCFperturb.f90 \
	PolSim.f90 \
	RTFquadratic.f90 \
	typedEmis_module.f90 \
	MWmeasErrModule.f90

# IDL removed from LIB because apparently not needed in the library
# IDL not needed except for the IDL targets that make .so files. 
# IDL targets are used in IDL .pro files in sensormodel, but not in 
# rslib, mwsensors, or irsensors
#	idlmvltlnd.f \
#	idlConstInit.f \
#	idlConstNoSCF.f \
#	idlPresToAlt.f \

LIBOBJS := \
	SensorData.o \
	NoiseModule.o \
	RTFquadratic.o \
	CalibSim.o \
	PointSCFperturb.o \
	PolSim.o \
	typedEmis_module.o \
	MWmeasErrModule.o

#	idlmvltlnd.o \
#	idlConstInit.o \
#	idlConstNoSCF.o \
#	idlPresToAlt.o \

else
LIBSRCS := \
	SensorData.f90 \
	NoiseModule.f90

LIBOBJS := \
	SensorData.o \
	NoiseModule.o
endif

#############################
####### build targets #######
#############################

LIB: clean $(LIB) 
ifneq ($(MACHINE),ppc64)
APP: ${GETCONST} ${IDLMV} ${IDLCONST} ${IDLCONSTFEW} ${IDLPRESTOALT}
else
APP:  ${IDLMV} ${IDLCONST} ${IDLCONSTFEW} ${IDLPRESTOALT}
endif

$(LIB): $(LIBOBJS)
	test  -d ${LIBPWD}/bin/${PLATFORM}${DBG} || mkdir ${LIBPWD}/bin/${PLATFORM}${DBG}
	ar -r $@ $? 
	test  -d ${LIBPWD}/include/${PLATFORM}${DBG} || mkdir ${LIBPWD}/include/${PLATFORM}${DBG}
	find . -name '*.mod' -exec cp {} ${LIBPWD}/include/${PLATFORM}${DBG}/. \;

${GETCONST}: getConst.o
	@#create the command directory if it does not exist
	test  -d ${LIBPWD}/bin/${PLATFORM}${DBG} || mkdir ${LIBPWD}/bin/${PLATFORM}${DBG}
	${LDR} ${FFLAGS} ${MFLAGS} ${RS_INCLUDES} ${MATLABOBJS} getConst.o ${RS_LIBS} -o ${GETCONST}

${IDLMV}: idlmvltlnd.o \
	${RSLIB_HOME}/mw/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/util/bin/${PLATFORM}${DBG}/*.a
	${LDR} ${FFLAGS} ${RS_INCLUDES} ${shared} idlmvltlnd.o \
	${RS_LIBS} -o ${IDLMV}

${IDLCONST}: idlConstInit.o ${RSLIB_HOME}/mw/io/SCFread.o \
	${RSLIB_HOME}/mw/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/util/bin/${PLATFORM}${DBG}/*.a
	${LDR} ${FFLAGS} ${RS_INCLUDES} ${shared} idlConstInit.o \
	${RSLIB_HOME}/mw/io/SCFread.o ${RS_LIBS} -o ${IDLCONST}

${IDLCONSTFEW}: idlConstNoSCF.o \
	${RSLIB_HOME}/mw/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/util/bin/${PLATFORM}${DBG}/*.a
	${LDR} ${FFLAGS} ${RS_INCLUDES} ${shared} idlConstNoSCF.o \
	${RS_LIBS} -o ${IDLCONSTFEW}

${IDLPRESTOALT}: idlPresToAlt.o ${RSLIB_HOME}/mw/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/util/bin/${PLATFORM}${DBG}/*.a
	${LDR} ${FFLAGS} ${RS_INCLUDES} ${shared} idlPresToAlt.o ${RS_LIBS} \
	-o ${IDLPRESTOALT}

MODULES := $(wildcard *.mod)
ifeq ($(strip $(MODULES)),)
MOD :=
RMMOD :=
else
MOD := ${LIBPWD}/include/${PLATFORM}${DBG}/${MODULES};
RMMOD := (cd ${LIBPWD}/include/${PLATFORM}${DBG}; /bin/rm -f ${MODULES};);
endif

clean:
	(/bin/rm -f *.o ${MODULES} ${LIB} $(GETCONST) $(IDLMV) $(IDLCONST) $(IDLCONSTFEW) $(IDLPRESTOALT));
	${RMMOD}
