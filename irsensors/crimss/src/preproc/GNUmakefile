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

IODIR := ../io

# Used for .mod files in compilation, preceding RS_INCLUDES:
INCLUDES:= -I$(SENSINCL_PATH)/$(PLATFORM)${DBG}

vpath %.f90 ${PWD}
vpath %.f90 ${IODIR}
vpath %.f90 ${RETDIR}
vpath %.f90 ${GENERICSRC}
vpath %.f ${GENERICSRC}


COMMAND_PATH:=../../bin/$(PLATFORM)
COMMAND:=${COMMAND_PATH}/preprocessor${DBG}
COMMAND_ANC:=${COMMAND_PATH}/BuildAncillary${DBG}
COMMAND_TSKIN:=${COMMAND_PATH}/ModTskinAnc${DBG}

IDLLIB:=${RSLIB_HOME}/idl/util/misc
SET_IDL_PATH:=./set_IDL_search_path

ASRC:= ${RETDIR}/Dimensions.f90 \
       ${IODIR}/GranuleReadModule.f90 \
       ${GENERICSRC}/BuildAncillary.f90

TSRC:= ${RETDIR}/Dimensions.f90 \
       ${GENERICSRC}/ModTskinAnc.f90

IOSRC:= \
    ${IODIR}/MWsensorReader.f90 \
    ${IODIR}/IRsensorReader.f90

GSRC:= \
	${GENERICSRC}/Preprocessor.f90

IOOBJ:= $(patsubst $(IODIR)/%.f90,%.o,$(IOSRC))

AOBJ:= $(patsubst ${RETDIR}/%.f90,%.o,$(ASRC))
AOBJ:= $(patsubst ${IODIR}/%.f90,%.o,$(AOBJ))
AOBJ:= $(patsubst ${GENERICSRC}/%.f90,%.o,$(AOBJ))

TOBJ:= $(patsubst ${RETDIR}/%.f90,%.o,$(TSRC))
TOBJ:= $(patsubst ${IODIR}/%.f90,%.o,$(TOBJ))
TOBJ:= $(patsubst ${GENERICSRC}/%.f90,%.o,$(TOBJ))

$(info AOBJ $(AOBJ) )
$(info TOBJ $(TOBJ) )

OBJ:=  $(IOOBJ) $(GSRC)

ALL: ${COMMAND} ${COMMAND_ANC} ${COMMAND_TSKIN} sweep ${SET_IDL_PATH}

${COMMAND}: ${OBJ} ${RSLIB_HOME}/mw/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/util/bin/${PLATFORM}${DBG}/*.a \
	${RSLIB_HOME}/ir/bin/${PLATFORM}${DBG}/*.a
	@#create the command directory if it does not exist
	test -d ${COMMAND_PATH} || mkdir -p ${COMMAND_PATH}
	${LDR} ${FFLAGS} ${RS_INCLUDES} ${OBJ} ${INCLUDES} ${RS_LIBS} -o $@

${COMMAND_ANC}: ${AOBJ}
	test -d ${COMMAND_PATH} || mkdir -p ${COMMAND_PATH}
	${LDR} ${FFLAGS} ${RS_INCLUDES} ${AOBJ} ${RS_LIBS} -o $@

${COMMAND_TSKIN}: ${TOBJ}
	test -d ${COMMAND_PATH} || mkdir -p ${COMMAND_PATH}
	${LDR} ${FFLAGS} ${RS_INCLUDES} ${TOBJ} ${RS_LIBS} -o $@

.PHONY: sweep
sweep:
	@#sweep the .mod files so another platform does not confuse them
	test -d ${SENSINCL_PATH}/${PLATFORM}${DBG} || mkdir -p ${SENSINCL_PATH}/${PLATFORM}${DBG}
	find . -name '*.mod' -exec mv {} ${SENSINCL_PATH}/${PLATFORM}${DBG}/. \;

${SET_IDL_PATH}:
	echo -n "\\!PATH = '"${IDLLIB}"' + ':' + \\!PATH" > ${SET_IDL_PATH}

*.o:    ${RSLIB_HOME}/mw/include/${PLATFORM}${DBG}/*.mod \
	${RSLIB_HOME}/util/include/${PLATFORM}${DBG}/*.mod \
	${RSLIB_HOME}/ir/include/${PLATFORM}${DBG}/*.mod

#Ensure that local .mod files are up to date
${GOBJ}: ${LOBJ} ${IOBJ}

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
