#!/bin/bash

echo "Starting job on " `date` #Only to display the starting of production date
echo "Running on " `uname -a` #Only to display the machine where the job is running
echo "System release " `cat /etc/redhat-release` #And the system release
echo "CMSSW on Condor"

# to get condor-chirp from CMSSW
PATH="/usr/libexec/condor:$PATH"
source /cvmfs/cms.cern.ch/cmsset_default.sh

tar -xzf ${4}.tar.gz
rm ${4}.tar.gz
cd ${4}/src
scram b ProjectRename
eval `scramv1 runtime -sh`

cd Phase2PixelSim/Phase2PixelStubs/python/

echo "cmsRun Phase2PixelStubs_cfg.py"
cmsRun Phase2PixelStubs_cfg.py  start=$1 finish=$2 geometry=$3

mv *.root ${_CONDOR_SCRATCH_DIR}
