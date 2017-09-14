import sys
from os import system, environ

import optparse 
import subprocess

system("rm *.tar.gz")
system("tar -zcf ${CMSSW_VERSION}.tar.gz -C ${CMSSW_BASE}/.. ${CMSSW_VERSION} --exclude=Phase2PixelSim/.git --exclude=Phase2PixelSim/Phase2PixelStubs/python/condor/logs --exclude=Phase2PixelSim/Phase2PixelStubs/python/condor/finishedJobs")

parser = optparse.OptionParser("usage: %prog [options]\n")

parser.add_option ('-n',  dest='numfile', type='int', default = 10, help="number of files per job")
parser.add_option ('-c',  dest='noSubmit', action='store_true', default = False, help="Do not submit jobs. Only create condor_submit.txt.")
parser.add_option ('-g',  dest='geometry', type='string', default = '', help="Specify geometry to run.")

options, args = parser.parse_args()

submitFile = '''universe = vanilla
Executable = condorJobs.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = condorJobs.sh, ../Phase2PixelStubs_cfg.py, CMSSW_9_2_3.tar.gz
x509userproxy = $ENV(X509_USER_PROXY)

'''

fileLoc = environ["CMSSW_BASE"] + "/src/Phase2PixelSim/Phase2PixelStubs/python/GeomRootFiles/"
fileName = "OT613_200_IT4025_opt8s3l.txt"
filePath = fileLoc+fileName
file = open(filePath)
counter = 0

for line in file:
    counter += 1

prodName = "Phase2PixelStubs"
filesPerJob = options.numfile
geomSample = "OT613_200_IT4025_opt8s3l"
sample = "opt6s3l"
outFileList = [submitFile]

for startFileNum in xrange(0, counter, filesPerJob):
    endFileNumber = startFileNum+filesPerJob
    outFileList.append("Arguments = {0} {1} {2} $ENV(CMSSW_VERSION)\n".format(startFileNum, endFileNumber, options.geometry))
    outFileList.append("Output = logs/{0}_{1}_{2}.stdout\n".format(prodName, sample, startFileNum))
    outFileList.append("Error = logs/{0}_{1}_{2}.stderr\n".format(prodName, sample, startFileNum))
    outFileList.append("Log = logs/{0}_{1}_{2}.log\n".format(prodName, sample, startFileNum))
    outFileList.append("Queue\n\n")

file.close()

fout = open("condor_submit.txt", "w")
fout.write(''.join(outFileList))
fout.close()

if not options.noSubmit: 
    system('mkdir -p logs')
    system("echo 'condor_submit condor_submit.txt'")
    system('condor_submit condor_submit.txt')
