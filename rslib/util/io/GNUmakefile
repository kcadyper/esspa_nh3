#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#makefile for rslib/mw/io
#$Name$
#$Id$
#Copyright AER, Inc., 2002, 2003. All rights Reserved.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include ../../makefile.common

# CVS_TAG is to be filled with a delivery tag at the release
CVS_TAG = $Name$
TB_STATUS = $(findstring release, $(CVS_TAG))

LIBPWD := ..
LIB := ${LIBPWD}/bin/${PLATFORM}${DBG}/libiomod.a

ifeq ($(TB_STATUS),)
LIBSRC :=  \
	ncdfUtil.f90 \
	MapBeamDefaults.f90 \
	scene_io_module.f \
	rad_io_module.f \
	bkg_io_module.f \
	sfcgrd_io_module.f90 \
	RDR_io_module.f90 \
	RadFile.f90 \
	ClusterFileIO.f90 \
	IOutilities.f90

LIBOBJ := \
	ncdfUtil.o \
	MapBeamDefaults.o \
	scene_io_module.o \
	rad_io_module.o \
	bkg_io_module.o \
	sfcgrd_io_module.o \
	RDR_io_module.o \
	RadFile.o \
	ClusterFileIO.o \
	IOutilities.o
else
LIBSRC :=  \
	ncdfUtil.f90 \
	MapBeamDefaults.f90 \
	scene_io_module.f \
	rad_io_module.f \
	bkg_io_module.f \
	sfcgrd_io_module.f90 \
	ClusterFileIO.f90 \
	IOutilities.f90
LIBOBJ := \
	ncdfUtil.o \
	MapBeamDefaults.o \
	scene_io_module.o \
	rad_io_module.o \
	bkg_io_module.o \
	sfcgrd_io_module.o \
	ClusterFileIO.o \
	IOutilities.o
endif

CSRC := ncdfUtil_c.c
COBJ := $(patsubst %.c,%.o, $(CSRC))

#############################
####### build targets #######
#############################

All: clean $(LIB)

$(LIB): $(LIBOBJ) $(COBJ)
	test  -d ${LIBPWD}/bin/${PLATFORM}${DBG} || mkdir ${LIBPWD}/bin/${PLATFORM}${DBG}
	ar -r $@ $?
	test  -d ${LIBPWD}/include/${PLATFORM}${DBG} || mkdir ${LIBPWD}/include/${PLATFORM}${DBG}
	find . -name '*.mod' -exec cp {} ${LIBPWD}/include/${PLATFORM}${DBG}/. \;

$(LIBOBJ): $(LIBSRC)

${COBJ}: $(CSRC)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(C_INCLUDES) $<

MODULES := $(wildcard *.mod)
ifeq ($(strip $(MODULES)),)
MOD :=
RMMOD :=
else
MOD := ${LIBPWD}/include/${PLATFORM}${DBG}/${MODULES};
RMMOD := (cd ${LIBPWD}/include/${PLATFORM}${DBG}; /bin/rm -f ${MODULES};);
endif

clean:
#	@echo "cleaning..." $(LIBOBJ) ${MODULES} ${MOD} ${LIB}
	(/bin/rm -f $(LIBOBJ) ${MODULES} ${LIB} $(COBJ));
	${RMMOD}





