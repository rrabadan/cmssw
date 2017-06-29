import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('HLTbits',
  triggerLabel = cms.InputTag("TriggerResults","","HLT")
)
