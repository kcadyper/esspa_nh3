#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Toplevel makefile for retr
#$Name$
#$Id$
#Copyright AER, Inc., 2002, 2003. All rights Reserved.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Important: you must set the base location of rslib.

RSLIB_HOME = //project/p1913/esspa11/rslib

export RSLIB_HOME

# load common rules
# force the retrieval code to use the file "makefile.common" from rslib

include ${RSLIB_HOME}/makefile.common

GENERICSRC = ${PWD}/../generic/src

SENSINCL_PATH:=$(PWD)/include
export SENSINCL_PATH

Default target = all: clean first--this prevents any unforseen dependencies
#from not being included.
all: sim retr directives preproc
	@if [ $(RSLIB_HOME) != "/project/rs_sw/lib/rslib" ]; \
	then echo ">>-------------------------------------<<" ; \
	echo ">> USING $(RSLIB_HOME)" ; \
	echo ">>-------------------------------------<<" ; fi

.PHONY: retr
retr: RETR

.PHONY: aux
aux: AUX

.PHONY: retr_sig
retr_sig: RETR_SIG

.PHONY: sim
sim: SIM

AUX:
	( cd ../generic/src; $(MAKE) );

AUXCLEAN:
	( cd ../generic/src; $(MAKE) clean; );

RETR:
	( cd src/retr; $(MAKE) GENERICSRC=$(GENERICSRC) RETDIR=$(GENERICSRC)/retr; );

RETR_SIG:
	( cd src/retr; $(MAKE) GENERICSRC=$(GENERICSRC) RETDIR=$(GENERICSRC)/retr -f GNUmakefile.retr_sig; );

SIM:
	( cd src/sim; $(MAKE) GENERICSRC=$(GENERICSRC); );

directives:
	python ./crimssrun --directives
	python ./crimssrun --ppdirectives

simclean:
	( cd src/sim; $(MAKE) clean; );

preproc:
	$(MAKE) -C ./src/preproc GENERICSRC=$(GENERICSRC) RETDIR=$(GENERICSRC)/retr

preprocclean:
	( $(MAKE) -C ./src/preproc  clean; );

retrclean:
	( cd src/retr; $(MAKE) clean; );

clean: retrclean simclean preprocclean AUXCLEAN
	(find . -iname '*.o' -delete)
	(rm -rfv  ./bin/$(PLATFORM)/*)

local:
	test -d src/lib || ln -s $(RSLIB) src/lib



