#!/usr/bin/env python
#$Id$
############################################################
# Script for CrIMSS NMP
#
# Developed by: AER September 2017
# Copyright: Atmospheric and Environmental Research, Inc., 2017
############################################################
# The code has been tested only on Linux

import sys
import os
from netCDF4 import Dataset
import subprocess
import stat

def runCode(command):
    homeDir=os.getcwd()
    fname = homeDir + '/run_oss.tcsh'
    with open(fname, 'w+') as fid:
        fid.write('#!/bin/tcsh\n')
        fid.write(command + '\n')

    # add permission
    statFile = os.stat(fname).st_mode
    os.chmod(fname, statFile | stat.S_IEXEC)

    proc = subprocess.Popen(fname)
    dummy = proc.communicate()
    os.chdir(homeDir)
    if proc.returncode == 0:
        os.remove(fname)
    return proc.returncode



if __name__ == '__main__':
    gfcDir = None
    l1bFile = None
    nwpDir = None

    if len(sys.argv) < 3:
        print "err<crimsNMP>: CrISl1bFile GFSin [NWPout]"
        sys.exit(1)

    if len(sys.argv) > 1:
        l1bFile = sys.argv[1]

    if len(sys.argv) > 2:
        gfcDir = sys.argv[2]

    if len(sys.argv) > 3:
        nwpDir = sys.argv[3]
    else:
        nwpDir = gfcDir

    print 'l1bFile:', l1bFile
    print 'gfcDir:', gfcDir
    print 'nwpDir:', nwpDir
    homeDir=os.getcwd()

    workDir=homeDir + '/src/preproc'
    os.chdir(workDir)
    print 'workDir', workDir

    if (not os.path.isfile(l1bFile)):
        print "err<crimsNMP>: CrISl1bFile not FOUND: ", l1bFile
        sys.exit(1)

    if (not os.path.isdir(gfcDir)):
        print "err<crimsNMP>: GFSin not FOUND: ", gfcDir
        sys.exit(1)

    if (not os.path.isdir(nwpDir)):
        print "wrn<crimsNMP>: nwpPath benn created: ", nwpDir
        os.makedirs(nwpDir)


    fid = Dataset(l1bFile)
    timeStart = fid.time_coverage_start[:]
    timeEnd = fid.time_coverage_end[:]
    fid.close()

    # check gfc data
    # gfc dir year/month/day
    gfcName = '/gfsanl_4_'+timeStart[:4] + timeStart[5:7] + timeStart[8:10]
    gfcPath = gfcDir + '/' + timeStart[:4] + '/' + timeStart[5:7] + '/' + timeStart[8:10]
    nwpPath = nwpDir + '/' + timeStart[:4] + '/' + timeStart[5:7] + '/' + timeStart[8:10]
    nwpOutPath = nwpPath

    if (not os.path.isdir(gfcPath)):
        print "err<crimsNMP>: gfcPath not FOUND: ", gfcPath
        sys.exit(1)

    if (not os.path.isdir(nwpPath)):
        print "wrn<crimsNMP>: nwpPath benn created: ", nwpPath
        os.makedirs(nwpPath)
    gfcPath += gfcName
    nwpPath += gfcName

    # gfc file every 6 hours with name
    #   gfsanl_4_YYYYMMDD_0[0,6,12,18]00_000.grb2
    hs = int((int(timeStart[11:13]) + int(timeStart[14:16]) / 60.) / 6.) * 6
    he = int((int(timeEnd[11:13]) + int(timeEnd[14:16]) / 60.) / 6.) * 6

    timeList = set([hs, hs + 6, he, he + 6])

    fileInList = []
    fileOutList = []
    for idx in timeList:
        cFileName = gfcPath + '_%02d00_000.grb2' % (idx)
        if (not os.path.isfile(cFileName)):
            print "err<crimsNMP>: required gfc file not FOUND: ", cFileName
            sys.exit(1)

        oFileName = nwpPath + '%02d.nc' % (idx)
        if (os.path.isfile(oFileName)):
            print "wrn<crimsNMP>: nwp file FOUND: ", oFileName
            continue

        fileInList.append(cFileName)
        fileOutList.append(oFileName)
    if len(fileInList) == 0:
        print 'Nothing to process'
        sys.exit(0)

    fpath = 'set_IDL_search_path'
    ff=open(fpath,'r')
    pathstr = ff.readline()
    ff.close()

    for idx in range(len(fileInList)):
        command = 'idl -e "'+ pathstr +' & run_gfs" -args "' + fileInList[idx] +'" "'+ nwpOutPath + '"'
        runCode(command)

    os.chdir(homeDir)
