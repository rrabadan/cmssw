import FWCore.ParameterSet.Config as cms

hltntupler = cms.EDAnalyzer('MuMuTrkHLT',
  #tagtriggerResult      = cms.untracked.InputTag("TriggerResults::HLT"),
  #tagtriggerSummary     = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
  triggerResult         = cms.untracked.InputTag("TriggerResults::TEST"),
  triggerSummary        = cms.untracked.InputTag("hltTriggerSummaryAOD::TEST"),
  #puInfoTag             = cms.untracked.InputTag("addPileupInfo"),
  #offlineVtx            = cms.InputTag("offlinePrimaryVertices"),
  #beamspot              = cms.InputTag("offlineBeamSpot"),
  genParticlesTag       = cms.untracked.InputTag("genParticles"),
  L3CandidatesTag       = cms.InputTag("hltL3MuonCandidates"), 
  TkCandidatesTag       = cms.InputTag("hltJpsiTkAllConeTracksIter"), 
  MuMuVtxTag            = cms.untracked.InputTag("hltDisplacedmumuVtxProducerDoubleMu4Jpsi"),
  V0VtxTag              = cms.untracked.InputTag("hltDisplacedtktkVtxProducerV0")
)
