import FWCore.ParameterSet.Config as cms

process = cms.Process("FPIX")

process.load("FWCore.MessageService.MessageLogger_cfi")


process.TFileService = cms.Service("TFileService",
   fileName = cms.string('phase2PixelStub.root')
)

process.load("Phase2PixelSim.Phase2PixelStubs.Phase2PixelStubs_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",

    fileNames = cms.untracked.vstring(
        '/store/user/lpcfpix/TTbar_14TeV_923_OT613_200_IT4025_opt8s3l/step3_PU200/170626_081843/0000/step3_762.root'
    )
)

Phase2PixelStubs = cms.EDAnalyzer("Phase2PixelStubs",
    src = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted")
)

process.p = cms.Path(process.Phase2PixelStubs)

