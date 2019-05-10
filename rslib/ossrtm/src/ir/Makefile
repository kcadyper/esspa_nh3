include ../../make.common

include ../common/make_list

COREPATH:=core
WRAPPATH:=wrapper

vpath %.f90 $(COREPATH)
vpath %.f90 $(WRAPPATH)

CORESRC:= \
	$(COREPATH)/asolar_source_function.f90 \
	$(COREPATH)/OSSIRmoduleSubs.f90 \
	$(COREPATH)/oss_ir_module.f90 \
	$(COREPATH)/oss_ir_module_scat.f90

WRAPCSRC:= 

COREOBJ:= $(patsubst $(COREPATH)/%.f90,%.o,$(CORESRC))
WRAPOBJ:= $(patsubst $(WRAPPATH)/%.f90,%.o,$(WRAPSRC))
OBJS:= ${COMOBJ} ${COREOBJ} ${WRAPOBJ}

LIB:=ossir
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
oss_ir_module.o: OSSIRmoduleSubs.o
oss_ir_module_scat.o : oss_ir_module.o OSSIRmoduleSubs.o
${COREOBJ}: ${COMOBJ}
${WRAPOBJ}: ${COREOBJ}

# clean up only the intermediate files
clean:
	(\rm -f *.o *.mod *.a *.so *~);
