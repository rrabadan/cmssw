import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing


#----------------------------------------------------------------------------

# Setup Settings for ONIA FOREST:

HLTProcess     = "HLT"    # Name of HLT process 
isMC           = False # if input is MONTECARLO: True or if it's DATA: False
muonSelection  = "Trk"    # Single muon selection: Glb(isGlobal), GlbTrk(isGlobal&&isTracker), Trk(isTracker) are availale

triggerList    = {
		# Double Muon Trigger List
		'DoubleMuonTrigger' : cms.vstring(
         'HLT_HIUPC_DoubleMu0_BptxAND_MaxPixelTrack_v1',
         'HLT_HIUPC_DoubleMu0_NotMBHF2AND_MaxPixelTrack_v1',
         'HLT_HIUPC_DoubleMu0_NotMBHF2AND_v1',
         'HLT_HIUPC_DoubleMu0_NotMBHF2OR_MaxPixelTrack_v1',
         'HLT_HIUPC_DoubleMu0_NotMBHF2OR_v1',
         'HLT_HIUPC_DoubleMuOpen_BptxAND_MaxPixelTrack_v1',
         'HLT_HIUPC_DoubleMuOpen_NotMBHF2OR_MaxPixelTrack_v1',
         'HLT_HIUPC_DoubleMuOpen_NotMBHF2OR_v1',
),
		# Double Muon Filter List
		'DoubleMuonFilter'  : cms.vstring(
         'hltMaxPixelTrackForUPC',
         'hltMaxPixelTrackForUPC',
         'hltL1sDoubleMu0NotMBHF2AND',
         'hltMaxPixelTrackForUPC',
         'hltL1sDoubleMu0NotMBHF2OR',
         'hltMaxPixelTrackForUPC',
         'hltMaxPixelTrackForUPC',
         'hltL1sDoubleMuOpenNotMBHF2OR',
			),
		# Single Muon Trigger List
	'SingleMuonTrigger' : cms.vstring(
        'HLT_HIUPC_SingleMu0_BptxAND_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleMu0_NotMBHF2AND_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleMu0_NotMBHF2AND_v1',
        'HLT_HIUPC_SingleMu0_NotMBHF2OR_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleMu0_NotMBHF2OR_v1',
        'HLT_HIUPC_SingleMu3_BptxAND_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleMu3_NotMBHF2OR_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleMu3_NotMBHF2OR_v1',
        'HLT_HIUPC_SingleMuOpen_BptxAND_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleMuOpen_NotMBHF2AND_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v1',
        'HLT_HIUPC_SingleMuOpen_NotMBHF2OR_MaxPixelTrack_v1',
        'HLT_HIUPC_SingleMuOpen_NotMBHF2OR_v1',
			),
	# Single Muon Filter List
	'SingleMuonFilter'  : cms.vstring(
        'hltMaxPixelTrackForUPC',
        'hltMaxPixelTrackForUPC',
        'hltL1sSingleMu0NotMBHF2AND',
        'hltMaxPixelTrackForUPC',
        'hltL1sSingleMu0NotMBHF2OR',
        'hltMaxPixelTrackForUPC',
        'hltMaxPixelTrackForUPC',
        'hltL1sSingleMu3NotMBHF2OR',
        'hltMaxPixelTrackForUPC',
        'hltMaxPixelTrackForUPC',
        'hltL1sSingleMuOpenNotMBHF2AND',
        'hltMaxPixelTrackForUPC',
        'hltL1sSingleMuOpenNotMBHF2OR',
			)
	}

if isMC:
	#globalTag = 'auto:run2_mc_GRun'
  globalTag = '103X_upgrade2018_realistic_HI_v6'
else:
  #globalTag = '101X_dataRun2_HLT_frozen_v6'
  globalTag = '103X_dataRun2_Express_v2'

#----------------------------------------------------------------------------


# set up process
process = cms.Process("TriggerAnalysis")

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# Input and Output File Names
options.outputFile = "OniaForest.root"
options.secondaryOutputFile = "Jpsi_DataSet.root"
options.inputFiles = [
      #'/store/group/phys_heavyions/jaebeom/JpsiMM_0_15/ReEmul_HLTMenuV35_JpsiLP_EALB/181023_181845/0009/step3_mc_9652.root'
      #'/store/express/HIRun2018/HIExpressPhysics/FEVT/Express-v1/000/326/262/00000/FEA00CBC-C1DD-1943-B92F-3F65F644A101.root'
      '/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/501/00000/37EF0469-9F24-0D45-9799-326C05EB3FEE.root'
      ]
options.maxEvents = -1 # -1 means all events
options.maxEvents = 5000 # -1 means all events

# Get and parse the command line arguments
options.parseArguments()

# load the Geometry and Magnetic Field
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')

# Global Tag:
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTag, '')

#----------------------------------------------------------------------------

# For OniaTree Analyzer
from HiAnalysis.HiOnia.oniaTreeAnalyzer_cff import oniaTreeAnalyzer
oniaTreeAnalyzer(process, muonTriggerList=triggerList, HLTProName=HLTProcess, muonSelection=muonSelection, useL1Stage2=True, isMC=isMC, outputFileName=options.outputFile)

#----------------------------------------------------------------------------

# For HLTBitAnalyzer
process.load("HLTrigger.HLTanalyzers.HLTBitAnalyser_cfi")
process.hltbitanalysis.HLTProcessName              = HLTProcess
process.hltbitanalysis.hltresults                  = cms.InputTag("TriggerResults","",HLTProcess)
process.hltbitanalysis.l1tAlgBlkInputTag           = cms.InputTag("hltGtStage2Digis","",HLTProcess)
process.hltbitanalysis.l1tExtBlkInputTag           = cms.InputTag("hltGtStage2Digis","",HLTProcess)
process.hltbitanalysis.gObjectMapRecord            = cms.InputTag("hltGtStage2ObjectMap","",HLTProcess)
process.hltbitanalysis.gmtStage2Digis              = cms.string("hltGtStage2Digis")
process.hltbitanalysis.caloStage2Digis             = cms.string("hltGtStage2Digis")
process.hltbitanalysis.UseL1Stage2                 = cms.untracked.bool(True)
process.hltbitanalysis.getPrescales                = cms.untracked.bool(False)
process.hltbitanalysis.getL1InfoFromEventSetup     = cms.untracked.bool(False)
process.hltbitanalysis.UseTFileService             = cms.untracked.bool(True)
process.hltbitanalysis.RunParameters.HistogramFile = cms.untracked.string(options.outputFile)
process.hltbitanalysis.RunParameters.isData        = cms.untracked.bool(not isMC)
process.hltbitanalysis.RunParameters.Monte         = cms.bool(isMC)
process.hltbitanalysis.RunParameters.GenTracks     = cms.bool(False)
#process.hltBitAna = cms.EndPath(process.hltbitanalysis)
if (HLTProcess == "HLT") :
	process.hltbitanalysis.l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis","","RECO")
	process.hltbitanalysis.l1tExtBlkInputTag = cms.InputTag("gtStage2Digis","","RECO")
	process.hltbitanalysis.gmtStage2Digis    = cms.string("gtStage2Digis")
	process.hltbitanalysis.caloStage2Digis   = cms.string("gtStage2Digis")

#----------------------------------------------------------------------------

# For HLTObject Analyzer
process.load("HeavyIonsAnalysis.EventAnalysis.hltobject_cfi")
process.hltobject.processName = cms.string(HLTProcess)
process.hltobject.treeName = cms.string(options.outputFile)
process.hltobject.loadTriggersFromHLT = cms.untracked.bool(False)
process.hltobject.triggerNames = triggerList['DoubleMuonTrigger'] + triggerList['SingleMuonTrigger']
process.hltobject.triggerResults = cms.InputTag("TriggerResults","",HLTProcess)
process.hltobject.triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","",HLTProcess)
#process.hltObjectAna = cms.EndPath(process.hltobject)

#----------------------------------------------------------------------------
process.load("HiAnalysis.HiOnia.HFMaxCaloTowerProducer_cfi")
process.hfmax     = cms.Path(process.hfmaxcalotower)

#Options:
process.source = cms.Source("PoolSource",
		fileNames = cms.untracked.vstring( options.inputFiles ),
		)
process.TFileService = cms.Service("TFileService", 
		fileName = cms.string( options.outputFile )
		)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
#process.schedule  = cms.Schedule( process.oniaTreeAna , process.hltBitAna , process.hltObjectAna )
process.schedule  = cms.Schedule( process.hfmax, process.oniaTreeAna )
