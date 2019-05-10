#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Toplevel makefile for rslib
#$Id$ 
#Copyright AER, Inc., 2002, 2003. All rights Reserved.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CVS_TAG is to be filled with a delivery tag at the release
CVS_TAG = $Name$
TB_STATUS = $(findstring release, $(CVS_TAG))

RSLIB_HOME = ${PWD}
export RSLIB_HOME

ifeq ($(TB_STATUS),)
  IREXTRALIB = IRRTM IRIO IRCORE IRNOISE
  IREXTRALIBCLEAN = irrtmclean ircoreclean irioclean irnoiseclean
  EXTRALIB = MWRTM_C
  EXTRALIBCLEAN = mwrtmclean_c
else
  IREXTRALIB = IRRTM IRIO IRCORE
  IREXTRALIBCLEAN = irrtmclean ircoreclean irioclean
  EXTRALIB = 
  EXTRALIBCLEAN = 
endif

platform :=$(shell uname)
ifeq ($(platform),IRIX64)
ALL: REF LAPACK MATH NCDF MISC STRUCT IO RETR INSTR MWIO MWRTM MWCORE NOISE \
	${IREXTRALIB} ${EXTRALIB}
CORE: REF LAPACK MATH NCDF MISC STRUCT IO RETR INSTR MWIO MWRTM MWCORE NOISE \
	${IREXTRALIB} ${EXTRALIB}
endif

ifeq ($(platform),Linux)
ALL: REF LAPACK MATH NCDF MISC STRUCT IO RETR INSTR MWIO MWRTM MWCORE NOISE \
	${IREXTRALIB} ${EXTRALIB} #\
#	MISCIDL INSTRAPP MWIOIDL NOISEAPP
CORE: REF LAPACK MATH NCDF MISC STRUCT IO RETR INSTR MWIO MWRTM MWCORE NOISE \
	${IREXTRALIB} ${EXTRALIB}
endif

.PHONY: struct
struct: STRUCT

.PHONY: ref
ref: REF

.PHONY: lapack
lapack: LAPACK

.PHONY: math
math: MATH

.PHONY: instr
instr: INSTR

.PHONY: instrapp
instrapp: INSTRAPP

.PHONY: io
io: IO

.PHONY: misc
misc: MISC

.PHONY: miscidl
miscidl: MISCIDL

.PHONY: ncdf
ncdf: NCDF

.PHONY: retr
retr: RETR

.PHONY: mwcore
mwcore: MWCORE

.PHONY: mwio
mwio: MWIO

.PHONY: mwioidl
mwioidl: MWIOIDL

.PHONY: mwrtm
mwrtm: MWRTM

.PHONY: mwrtm_c
mwrtm_c: MWRTM_C

.PHONY: irrtm
irrtm: IRRTM

.PHONY: ircore
ircore: IRCORE

.PHONY: irio
irio: IRIO

.PHONY: noise
noise: NOISE

.PHONY: noiseapp
noiseapp: NOISEAPP

.PHONY: irnoise
irnoise: IRNOISE

STRUCT:
	( cd util/struct ; $(MAKE); );

structclean:
	( cd util/struct ; $(MAKE) clean; );

REF:
	( cd util/reference ; $(MAKE); );

refclean:
	( cd util/reference ; $(MAKE) clean; );

LAPACK:
	( cd util/math/lapack_subset ; $(MAKE); );

lapackclean:
	( cd util/math/lapack_subset ; $(MAKE) clean; );

MATH:
	( cd util/math ; $(MAKE); );

mathclean:
	( cd util/math ; $(MAKE) clean; );

IO:
	( cd util/io ; $(MAKE); );

ioclean:
	( cd util/io ; $(MAKE) clean; );

INSTR:
	( cd util/sensor_model ; $(MAKE) LIB; );

instrclean:
	( cd util/sensor_model ; $(MAKE) clean; );

INSTRAPP:
	( cd util/sensor_model ; $(MAKE) APP; );

MISC:
	( cd util/misc ; $(MAKE) LIB; );

miscclean:
	( cd util/misc ; $(MAKE) clean; );

MISCIDL:
	( cd util/misc ; $(MAKE) IDL; );

NCDF:
	( cd util/netCDF ; $(MAKE); );

ncdfclean:
	( cd util/netCDF ; $(MAKE) clean; );

RETR:
	( cd util/retr ; $(MAKE); );

retrclean:
	( cd util/retr ; $(MAKE) clean; );

MWIO:
	( cd mw/io ; $(MAKE) LIB; );

mwioclean:
	( cd mw/io ; $(MAKE) clean; );

MWIOIDL:
	( cd mw/io ; $(MAKE) IDL; );

MWCORE:
	( cd mw/core ; $(MAKE); );

mwcoreclean:
	( cd mw/core ; $(MAKE) clean; );

MWRTM:
	( cd mw/rtm ; $(MAKE); );

mwrtmclean:
	( cd mw/rtm ; $(MAKE) clean; );

MWRTM_C:
	( cd c/mw/rtm ; $(MAKE); );

mwrtmclean_c:
	( cd c/mw/rtm ; $(MAKE) clean; );

IRRTM:
	( cd ir/rtm ; $(MAKE); );

irrtmclean:
	( cd ir/rtm ; $(MAKE) clean; );

IRCORE:
	( cd ir/core ; $(MAKE); );

ircoreclean:
	( cd ir/core ; $(MAKE) clean; );

IRIO:
	( cd ir/io ; $(MAKE); );

irioclean:
	( cd ir/io ; $(MAKE) clean; );

NOISE:
	( cd mw/sensor_model ; $(MAKE) LIB; );

noiseclean:
	( cd mw/sensor_model ; $(MAKE) clean; );

NOISEAPP:
	( cd mw/sensor_model ; $(MAKE) APP; );

IRNOISE:
	( cd ir/sensor_model ; $(MAKE); );

irnoiseclean:
	( cd ir/sensor_model ; $(MAKE) clean; );

clean:  refclean lapackclean mathclean miscclean structclean ncdfclean retrclean mwcoreclean \
	mwioclean mwrtmclean ioclean instrclean noiseclean $(EXTRALIBCLEAN) $(IREXTRALIBCLEAN)
