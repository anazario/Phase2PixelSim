import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("FPIX")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')

process.TFileService = cms.Service("TFileService",
   fileName = cms.string('phase2PixelStub.root')
)

process.load("Phase2PixelSim.Phase2PixelStubs.Phase2PixelStubs_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

mylist = FileUtils.loadListFromFile('rootfiles.txt')
mylist.extend ( FileUtils.loadListFromFile('rootfiles.txt'))
readFiles = cms.untracked.vstring( *mylist) 

process.source = cms.Source("PoolSource", fileNames = readFiles)

Phase2PixelStubs = cms.EDAnalyzer("Phase2PixelStubs",
    src = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted")
)

process.p = cms.Path(process.Phase2PixelStubs)
