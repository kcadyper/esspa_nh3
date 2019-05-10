
#
# Makefile to compile exec_oss
#
# Yuguang He, AER, September, 2015
#

include make.common

LIBS = \
	-L${netCDF}/lib -lnetcdf_c++ -lnetcdf \
	-L$(MAIN_DIR)/fortran -loss

LIBPATH:=./lib/${PLATFORM}
INCLPATH:=./include/${PLATFORM}

#
# Default target to build everything.
#
#all: os_check ossforirlib osscppirlib ossidlirlib \
#ossformwlib osscppmwlib ossidlmwlib \
#ossforirtest osscppirtest ossformwtest osscppmwtest
all: ossforirlib ossformwlib
for: ossforirlib ossformwlib
#for: os_check ossforirlib ossformwlib ossforirtest ossformwtest
#cpp: os_check osscppirlib osscppmwlib osscppirtest osscppmwtest
#idl: os_check ossidlirlib ossidlmwlib

#
# Check to make sure we are on a valid orperating system.
#
#commented out because there's a bug in the first expression
#os_check:
#	@if [ $(OSFOUND) != yes ]; \
#	then echo Unknown Operating System name: $(platform); \
#	exit 9; \
#	fi

#
# Build all object files and put them into libraries.
#
ossforirlib:
	@cd src/ir; $(MAKE) flib

ossformwlib:
	@cd src/mw; $(MAKE) flib

osscpp:
	@cd src/interfaces/cc; $(MAKE) clib

ossidl:
	@cd src/ir; $(MAKE) slib
	@cd src/mw; $(MAKE) slib

#
# Clean up the *.o, *.mod, *.a, and *.exe files.
#
clean:
	@cd src/ir; $(MAKE) clean
	@cd src/mw; $(MAKE) clean
	\rm -rf $(LIBPATH)
	\rm -rf $(INCLPATH)
