#!/bin/bash

echo "Starting job on " `date` #Only to display the starting of production date
echo "Running on " `uname -a` #Only to display the machine where the job is running
echo "System release " `cat /etc/redhat-release` #And the system release
echo "CMSSW on Condor"

# to get condor-chirp from CMSSW
PATH="/usr/libexec/condor:$PATH"
source /cvmfs/cms.cern.ch/cmsset_default.sh

tar -xzvf CMSSW_9_2_3.tgz
rm CMSSW_9_2_3.tgz
cd CMSSW_9_2_3/src/
scram b ProjectRename
eval `scramv1 runtime -sh`

cd Phase2PixelSim/Phase2PixelStubs/python/

echo "cmsRun Phase2PixelStubs_cfg.py"
cmsRun Phase2PixelStubs_cfg.py -s $1 -e $2

pwd
ls -altrh

cd ${_CONDOR_SCRATCH_DIR}

pwd
ls -altrh

cp $CMSSW_BASE/src/Phase2PixelSim/Phase2PixelStubs/python/*.root .

