#-------------------------------------------------------------------------------
# $Name$
# $Id$
# Copyright AER, Inc., 2002, 2003. All rights Reserved.
#-------------------------------------------------------------------------------

#makefile for determining env data and land emissivity (GetEnvData and SetLandEmiss)

# load common rules and environment variables
include ${RSLIB_HOME}/makefile.common

COMMAND_PATH:=../bin/$(PLATFORM)
RETDIR:=./retr/
IODIR:=../../crimss/src/io
#append _debug=${DBG} to command name if compiled with the debug option

COMMAND:=${COMMAND_PATH}/getEmisData${DBG}

fSRC:=
f90SRC:= GetEmisData.f90

f90GSRC:=$(patsubst %.f90,%.o,$(GSRC))

RS_INCLUDES:= -I. $(RS_INCLUDES)

all: mkdir ${COMMAND}

mkdir:
	@#create the command directory if it does not exist
	test -d ${COMMAND_PATH} || mkdir -p ${COMMAND_PATH}

${COMMAND}: mkdir  \
	${RSLIB_HOME}/mw/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/util/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/ir/bin/${PLATFORM}${DBG}/*.a
	${LDR} ${FFLAGS} ${RS_INCLUDES} ${f90SRC} ${RS_LIBS} -o $@


*.o: ${RSLIB_HOME}/mw/include/${PLATFORM}${DBG}/*.mod \
	${RSLIB_HOME}/util/include/${PLATFORM}${DBG}/*.mod \
	${RSLIB_HOME}/ir/include/${PLATFORM}${DBG}/*.mod

clean:
	(find . \( -iname '*.o' -o -iname '*.mod' \) -delete)
	(\rm -f ${COMMAND} );

%.o: %.f90
	$(FC) -c $(FFLAGS) $(INCLUDES) $(RS_INCLUDES) $< -o $@


