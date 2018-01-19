import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

from geomDict import dict

from os import system, environ
import sys
import optparse 

## --------------------------
## -- Command line options --
## --------------------------                                                                                                            
                   
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('standard')

options.register('start', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Indicates on which file to start.")
options.register('finish', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Indicates on which file to stop.")
options.register('geometry', 'opt8s4l', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Indicate geometry to run on.")

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
   fileName = cms.string("Phase2PixelStubs_{0}_{1}.root".format(options.geometry, options.start))
)

process.load("Phase2PixelSim.Phase2PixelStubs.Phase2PixelStubs_cfi")
process.load("Phase2PixelSim.Phase2PixelStubs.GeomConf.OT613_200_IT4025_{0}_cfi".format(options.geometry))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#get root files from list
FileList = FileUtils.loadListFromFile(dict[options.geometry])
size = len(FileList)

startNum = options.start
endNum = options.finish

#create partial list and run the corresponding jobs 
partialFileList = FileList[startNum:endNum]
readFiles = cms.untracked.vstring(*partialFileList) 
process.source = cms.Source("PoolSource", fileNames = readFiles)
process.p = cms.Path(process.Phase2PixelStubs)

