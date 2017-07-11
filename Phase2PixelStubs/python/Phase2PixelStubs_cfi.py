import FWCore.ParameterSet.Config as cms

Phase2PixelStubs = cms.EDAnalyzer('Phase2PixelStubs',

  TTStubs = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
)
