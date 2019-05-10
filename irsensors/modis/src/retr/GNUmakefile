#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# makefile for cmis/src/retr
#$Name$ 
#$Id$ 
#Copyright AER, Inc., 2002, 2003. All rights Reserved.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Notes: 
# 1. Need to replace getscene.f with getscene_old.f until I
#    integrate new version of getscene.f that utilizes new ncdf routines.
#    currently they reside on the same directory so I needed temporary way
#    to distinguish the two.

# load common rules
include ${RSLIB_HOME}/makefile.common

INSDIR:= ../sim

# Used for .mod files in compilation, preceding RS_INCLUDES:
INCLUDES:= -I. -I$(SENSINCL_PATH)/$(PLATFORM)${DBG}

vpath %.f90 ${GENERICSRC}
vpath %.f ${GENERICSRC}
vpath %.f90 ${INSDIR}
vpath %.f ${INSDIR}
vpath %.f90 core
vpath %.f core

COMMAND_PATH:=../../bin/$(PLATFORM)
#append _debug=${DBG} to command name if compiled with the debug option
COMMAND:=${COMMAND_PATH}/modis${DBG} 
COMMAND_LIN:=${COMMAND_PATH}/lin_modis${DBG} 
COMMAND_EMS:=${COMMAND_PATH}/GetEmisData${DBG} 

LSRC:= \
	core/classatm.f \
	core/setSfc.f90 \
	core/RadFile.f90 \
	core/CoefFile.f90 \
	core/OriginalFile.f90 \
	core/WriteFile.f90

ISRC:= \
        ${INSDIR}/IRNoiseModule.f90 \
        ${INSDIR}/IRChannels.f90 
GSRC:= ${GENERICSRC}/retr.f90
GSRC_LIN:= ${GENERICSRC}/LinRetrGroup.f90
GSRC_EMS:= ${GENERICSRC}/GetEmisData.f90

LOBJ:= $(patsubst core/%.f90,%.o,$(patsubst core/%.f,%.o,$(LSRC)))
IOBJ:= $(patsubst ${INSDIR}/%.f90,%.o,$(patsubst ${INSDIR}/%.f,%.o,$(ISRC)))
GOBJ:= $(patsubst ${GENERICSRC}/%.f90,%.o,$(patsubst ${GENERICSRC}/%.f,%.o,$(GSRC)))

OBJ:= $(LOBJ) $(IOBJ) $(GOBJ) 

GOBJ_LIN:= $(patsubst ${GENERICSRC}/%.f90,%.o,$(GSRC_LIN))
OBJ_LIN:= $(LOBJ) $(GOBJ_LIN)

GOBJ_EMS:= $(patsubst ${GENERICSRC}/%.f90,%.o,$(GSRC_EMS))
OBJ_EMS:= $(LOBJ) $(GOBJ_EMS)

ALL: ${COMMAND} ${COMMAND_LIN} ${COMMAND_EMS} sweep 


${COMMAND}: ${OBJ} ${RSLIB_HOME}/mw/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/util/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/ir/bin/${PLATFORM}${DBG}/*.a
	@#create the command directory if it does not exist
	test -d ${COMMAND_PATH} || mkdir -p ${COMMAND_PATH}
	${LDR} ${FFLAGS} ${RS_INCLUDES} ${OBJ} ${RS_LIBS} -o ${COMMAND} 

${COMMAND_LIN}: ${OBJ_LIN} ${RSLIB_HOME}/mw/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/util/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/ir/bin/${PLATFORM}${DBG}/*.a
	@#create the command directory if it does not exist
	test -d ${COMMAND_PATH} || mkdir -p ${COMMAND_PATH}
	${LDR} ${FFLAGS} ${RS_INCLUDES} ${OBJ_LIN} ${RS_LIBS} -o ${COMMAND_LIN} 

${COMMAND_EMS}: ${OBJ_EMS} ${RSLIB_HOME}/mw/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/util/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/ir/bin/${PLATFORM}${DBG}/*.a
	@#create the command directory if it does not exist
	test -d ${COMMAND_PATH} || mkdir -p ${COMMAND_PATH}
	${LDR} ${FFLAGS} ${RS_INCLUDES} ${OBJ_EMS} ${RS_LIBS} -o ${COMMAND_EMS} 

.PHONY: sweep
sweep:
	@echo +++ Compiled ${COMMAND_LIN}
	@#sweep the .mod files so another platform does not confuse them
	test -d ${SENSINCL_PATH}/${PLATFORM}${DBG} || mkdir -p ${SENSINCL_PATH}/${PLATFORM}${DBG}
	find . -name '*.mod' -exec mv {} ${SENSINCL_PATH}/${PLATFORM}${DBG}/. \;

*.o:    ${RSLIB_HOME}/mw/include/${PLATFORM}${DBG}/*.mod \
	${RSLIB_HOME}/util/include/${PLATFORM}${DBG}/*.mod \
	${RSLIB_HOME}/ir/include/${PLATFORM}${DBG}/*.mod

#Ensure that local .mod files are up to date
${GOBJ}: ${LOBJ} ${IOBJ}
${GOBJ_LIN}: ${LOBJ}
${GOBJ_EMS}: ${LOBJ}

MODULES:= $(wildcard ${SENSINCL_PATH}/${PLATFORM}${DBG}/*.mod)
ifeq ($(strip $(MODULES)),)
RMMOD :=
else
RMMOD := (cd ${SENSINCL_PATH}/${PLATFORM}${DBG}; /bin/rm -f *.mod;);
endif

# clean up only the intermediate files
clean:
	( /bin/rm  -f *.o $(wildcard *.mod) );
	${RMMOD}
