import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

import sys
import optparse 

#options for command line                                                                                                                                  
parser = optparse.OptionParser("usage: %prog [options]\n")

parser.add_option ('-s',  dest='startFile', type='int', default = 0, help="Indicates on which file from list to start running.")
parser.add_option ('-e',  dest='endFile', type='int', default = 5, help="Indicates on which file to stop.")

options, args = parser.parse_args()

#configuration file begin
process = cms.Process("FPIX")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')

process.TFileService = cms.Service("TFileService",
   fileName = cms.string("Phase2PixelStubs_{0}.root".format(options.startFile))
)

process.load("Phase2PixelSim.Phase2PixelStubs.Phase2PixelStubs_cfi")
process.load("Phase2PixelSim.Phase2PixelStubs.OT613_200_IT4025_opt8s3l_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#get root files from list
FileList = FileUtils.loadListFromFile('GeomRootFiles/OT613_200_IT4025_opt8s3l.txt')
size = len(FileList)

#verify user options input
if (options.startFile >= 0 and options.endFile <= size):
    startNum = options.startFile
    endNum = options.endFile
else:
    sys.exit("Input values out of bounds! Range of list goes from 0 to " + str(size) + ".")

#create partial list and run the corresponding jobs 
partialFileList = FileList[startNum:endNum]
readFiles = cms.untracked.vstring(*partialFileList) 
process.source = cms.Source("PoolSource", fileNames = readFiles)

Phase2PixelStubs = cms.EDAnalyzer("Phase2PixelStubs",
    src = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
    TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
    TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth"),
)

process.p = cms.Path(process.Phase2PixelStubs)

