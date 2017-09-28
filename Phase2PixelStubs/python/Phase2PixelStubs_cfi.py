import FWCore.ParameterSet.Config as cms

Phase2PixelStubs = cms.EDAnalyzer('Phase2PixelStubs',

  TTStubs = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
  TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
  TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth"),
  L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),               ## TTTrack input
  MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"), ## MCTruth input 
)


