/*! \file
Read a simple file, with some of the features of netCDF-4.

This is a very simple example which demonstrates some of the new
features of netCDF-4.0.

This example reads a simple file created by simple_nc4_wr.c. This is
intended to illustrate the use of the netCDF-4 C API.

This is part of the netCDF package. Copyright 2006-2011 University
Corporation for Atmospheric Research/Unidata. See COPYRIGHT file for
conditions of use. Full documentation of the netCDF can be found at
http://www.unidata.ucar.edu/software/netcdf/docs.
*/
#include <stdio.h>
#include <netcdf.h>
#include <stdlib.h>
#include <string.h>

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
const int ERRCODE = 2;

void callNetCDF(int retval) {
   if (retval != NC_NOERR) {
   	printf("NETCDF Error: %s\n", nc_strerror(retval));
   	exit(ERRCODE);
   }
}

void readCharAttribute(int ncid, char* atrName, char* buf) {

   nc_type nAttrType = NC_NAT;
   size_t nAttrLen = 0;

   callNetCDF(nc_inq_att(ncid, NC_GLOBAL, atrName, &nAttrType, &nAttrLen));

   char *ppszTemp[1];

   callNetCDF(nc_get_att_string(ncid, NC_GLOBAL, atrName, ppszTemp));
   strcpy(buf, ppszTemp[0]);
   nc_free_string(nAttrLen, ppszTemp);
}

void readCharData(int ncid, char* varName, char* buf) {
   int varId;

   callNetCDF(nc_inq_varid(ncid, varName, &varId));
   callNetCDF(nc_get_var(ncid, varId, buf));
}

void readCharData1D(int ncid, int varId, char* buf) {

   callNetCDF(nc_get_var(ncid, varId, buf));
}

