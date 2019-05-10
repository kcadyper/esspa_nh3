#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Makefile for rslib/util/rtm
#$Name$ 
#$Id$ 
#Copyright AER, Inc., 2016-, All rights Reserved.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
include ../../makefile.common

LIBPWD := ..
LIB := ${LIBPWD}/bin/${PLATFORM}${DBG}/libstruct.a 

SRC = \
	ControlStructure.f90 \
	AncillaryStructure.f90 \
	ImgObsStructure.f90 \
	IRobsStructure.f90 \
	MWobsStructure.f90 \
	OutputStructure.f90

LIBOBJ := $(patsubst %.f90,%.o, $(SRC))

#############################
####### build targets #######
#############################

ALL: clean $(LIB)

$(LIB): $(LIBOBJ)
	test  -d ${LIBPWD}/bin/${PLATFORM}${DBG} || mkdir -p ${LIBPWD}/bin/${PLATFORM}${DBG}
	ar -r $@ $? 
	test  -d ${LIBPWD}/include/${PLATFORM}${DBG} || mkdir -p ${LIBPWD}/include/${PLATFORM}${DBG}
	find . -name '*.mod' -exec cp {} ${LIBPWD}/include/${PLATFORM}${DBG}/. \;

$(LIBOBJ): $(LIBSRC)

MODULES := $(wildcard *.mod)
ifeq ($(strip $(MODULES)),)
MOD :=
RMMOD :=
else
MOD := ${LIBPWD}/include/${PLATFORM}${DBG}/${MODULES};
RMMOD := (cd ${LIBPWD}/include/${PLATFORM}${DBG}; /bin/rm -f ${MODULES};);
endif

clean:
	(/bin/rm -f $(LIBOBJ) ${MODULES} ${LIB});
	${RMMOD}
