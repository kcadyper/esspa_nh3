#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Toplevel makefile for retr
#$Name$ 
#$Id$ 
#Copyright AER, Inc., 2002, 2003. All rights Reserved.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Important: you must set the base location of rslib.

RSLIB_HOME = /project/rs_sw/lib/rslib
export RSLIB_HOME

# load common rules
# force the retrieval code to use the file "makefile.common" from rslib

include ${RSLIB_HOME}/makefile.common

GENERICSRC = ../../../generic/src

SENSINCL_PATH:=$(PWD)/include
export SENSINCL_PATH

Default target = all: clean first--this prevents any unforseen dependencies
#from not being included.
all: aux sim retr

.PHONY: aux
aux: AUX

.PHONY: retr
retr: RETR

.PHONY: sim
sim: SIM

AUX:
	( cd ../generic/src; $(MAKE); )
RETR:
	( cd src/retr; $(MAKE) GENERICSRC=$(GENERICSRC); );

SIM:
	( cd src/sim; $(MAKE) GENERICSRC=$(GENERICSRC); );

auxclean: 
	( cd ../generic/src; $(MAKE) clean; );

simclean: 
	( cd src/sim; $(MAKE) clean; );

retrclean: 
	( cd src/retr; $(MAKE) clean; );

clean: auxclean retrclean simclean

local: 
	test -d src/lib || ln -s $(RSLIB) src/lib



