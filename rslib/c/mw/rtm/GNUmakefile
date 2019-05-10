#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Makefile for rslib/c/mw/rtm/
#$Name$ 
#$Id$ 
#Copyright AER, Inc., 2002, 2003. All rights Reserved.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
include ../../../makefile.common

LIBPWD := ..
LIB := ${LIBPWD}/bin/${PLATFORM}${DBG}/libmwrtm_c.a 

LIBOBJS := $(patsubst %.c,%.o, $(wildcard *.c))

#############################
####### build targets #######
#############################

ALL: clean $(LIB)

$(LIB): $(LIBOBJS)
	test  -d ${LIBPWD}/bin/${PLATFORM}${DBG} || mkdir ${LIBPWD}/bin/${PLATFORM}${DBG}
	ar -r $@ $? 

clean:
	(/bin/rm -f $(LIBOBJS) ${LIB});
