import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

from os import system, environ
import sys
import optparse 

## --------------------------
## -- Geometry Dictionary ---
## --------------------------
dict = {}
dict['opt6s3l'] = environ["CMSSW_BASE"]+"/src/Phase2PixelSim/Phase2PixelStubs/python/GeomRootFiles/OT613_200_IT4025_opt6s3l.txt"
dict['opt7s4l'] = environ["CMSSW_BASE"]+"/src/Phase2PixelSim/Phase2PixelStubs/python/GeomRootFiles/OT613_200_IT4025_opt7s4l.txt"
dict['opt8s3l'] = environ["CMSSW_BASE"]+"/src/Phase2PixelSim/Phase2PixelStubs/python/GeomRootFiles/OT613_200_IT4025_opt8s3l.txt"

## --------------------------
## -- Command line options --
## --------------------------                                                                                                                                
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('standard')

options.register('startFile', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Indicates on which file to start.")
options.register('endFile', 10, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Indicates on which file to stop.")
options.register('geometry', 'opt8s3l', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Indicate geometry to run on.")

options.parseArguments()
options._tagOrder =[]

## --------------------------
## --- Config File Begin ----
## --------------------------
process = cms.Process("FPIX")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')

process.TFileService = cms.Service("TFileService",
   fileName = cms.string("Phase2PixelStubs_{0}_{1}.root".format(options.geometry, options.startFile))
)

process.load("Phase2PixelSim.Phase2PixelStubs.Phase2PixelStubs_cfi")
process.load("Phase2PixelSim.Phase2PixelStubs.OT613_200_IT4025_{0}_cfi".format(options.geometry))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#get root files from list
FileList = FileUtils.loadListFromFile(dict[options.geometry])
size = len(FileList)

#verify user options input
if (options.startFile >= 0 and options.endFile <= size):
    startNum = options.startFile
    endNum = options.endFile
elif (options.startFile == size - options.endFile and options.endFile < size):
    startNum = options.startFile
    endNum = size - options.endFile
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

