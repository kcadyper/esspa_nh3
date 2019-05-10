include ../../make.common

include ../common/make_list

COREPATH:=core
WRAPPATH:=wrapper

vpath %.f90 $(COREPATH)
vpath %.f90 $(WRAPPATH)

CORESRC:= \
	$(COREPATH)/netCDFhelper.f90 \
	$(COREPATH)/AngleUtil.f90 \
	$(COREPATH)/EmisTab.f90 \
	$(COREPATH)/OSSMWmoduleSubs.f90 \
	$(COREPATH)/oss_mw_module.f90 \
	$(COREPATH)/oss_mw_module_scat.f90 \
	$(COREPATH)/MwSfcModule.f90

WRAPCSRC:= 

COREOBJ:= $(patsubst $(COREPATH)/%.f90,%.o,$(CORESRC))
WRAPOBJ:= $(patsubst $(WRAPPATH)/%.f90,%.o,$(WRAPSRC))
OBJS:= ${COMOBJ} ${COREOBJ} ${WRAPOBJ}

LIB:=ossmw
FLIB := libf${LIB}.a
SLIB := lib${LIB}.so
LIBPATH:=../../lib/${PLATFORM}
INCLPATH:=../../include/${PLATFORM}
INCLUDES:= \
	   ${INCLUDES} \
	   -I${INCLPATH}

ALL: flib slib

flib: ${OBJS}
	ar -r ${FLIB} $? 
#	-ranlib ${RANOPTS} ${FLIB}
	test -d ${LIBPATH} || mkdir -p ${LIBPATH}
	cp -f ${FLIB} ${LIBPATH}
	test -d ${INCLPATH} || mkdir -p ${INCLPATH}
	cp -f *.mod ${INCLPATH}

slib: ${SOBJS}
	${LDR} -shared -fPIC -o ${SLIB} $? 
	test -d ${LIBPATH} || mkdir -p ${LIBPATH}
	cp -f ${SLIB} ${LIBPATH}
	test -d ${INCLPATH} || mkdir -p ${INCLPATH}
	cp -f *.mod ${INCLPATH}

#Dependencies
oss_mw_module.o oss_mw_module_scat.o: OSSMWmoduleSubs.o
${COREOBJ}: ${COMOBJ}
${WRAPOBJ}: ${COREOBJ}

# clean up only the intermediate files
clean:
	(\rm -f *.o *.mod *.a *.so *~);
