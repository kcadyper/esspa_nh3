#-------------------------------------------------------------------------------
# $Name$ 
# $Id$ 
# Copyright AER, Inc., 2002, 2003. All rights Reserved.
#-------------------------------------------------------------------------------

#makefile for simulation part of retrieval code (sim)

# Notes: 
# 1. Need to replace getscene.f with getscene_old.f until I
#    integrate new version of getscene.f that utilizes new ncdf routines.
#    currently they reside on the same directory so I needed temporary way
#    to distinguish the two.

# load common rules and environment variables
include ${RSLIB_HOME}/makefile.common

vpath %.f90 ${GENERICSRC}
vpath %.f ${GENERICSRC}

COMMAND_PATH:=../../bin/$(PLATFORM)
#append _debug=${DBG} to command name if compiled with the debug option
COMMAND:=${COMMAND_PATH}/sim_modis${DBG}
COMMAND_SCT:=${COMMAND_PATH}/simScat_modis${DBG}

LSRC:= IRNoiseModule.f90 IRChannels.f90 	
LOBJ:= $(patsubst %.f90,%.o,$(patsubst %.f,%.o,$(LSRC)))

GSRC:= ${GENERICSRC}/sim.f90
GOBJ:= $(patsubst ${GENERICSRC}/%.f90,%.o,$(patsubst ${GENERICSRC}/%.f,%.o,$(GSRC)))

OBJ:= ${LOBJ} ${GOBJ}

GSRC_SCT:= ${GENERICSRC}/simScat.f90
GOBJ_SCT:= $(patsubst ${GENERICSRC}/%.f90,%.o,$(patsubst ${GENERICSRC}/%.f,%.o,$(GSRC_SCT)))

OBJ_SCT:= ${LOBJ} ${GOBJ_SCT}

ALL:  ${COMMAND} ${COMMAND_SCT} 

${COMMAND}: ${OBJ} ${RSLIB_HOME}/mw/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/util/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/ir/bin/${PLATFORM}${DBG}/*.a
	@#create the command directory if it does not exist
	test -d ${COMMAND_PATH} || mkdir -p ${COMMAND_PATH}
	${LDR} ${FFLAGS} ${RS_INCLUDES} ${OBJ} ${RS_LIBS} -o ${COMMAND}

${COMMAND_SCT}: ${OBJ_SCT} ${RSLIB_HOME}/mw/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/util/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/ir/bin/${PLATFORM}${DBG}/*.a
	@#create the command directory if it does not exist
	test -d ${COMMAND_PATH} || mkdir -p ${COMMAND_PATH}
	${LDR} ${FFLAGS} ${RS_INCLUDES} ${OBJ_SCT} ${RS_LIBS} -o ${COMMAND_SCT}

*.o:    ${RSLIB_HOME}/mw/include/${PLATFORM}${DBG}/*.mod \
	${RSLIB_HOME}/util/include/${PLATFORM}${DBG}/*.mod \
	${RSLIB_HOME}/ir/include/${PLATFORM}${DBG}/*.mod

${GOBJ} ${GOBJ_SCT}: ${LOBJ}


# clean up only the intermediate files
clean:
	(\rm  -f *.o *.mod);
