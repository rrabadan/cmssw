import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing


#----------------------------------------------------------------------------

# Setup Settings for ONIA FOREST:

HLTProcess     = "MyWonderfulHLT"    # Name of HLT process 
isMC           = False # if input is MONTECARLO: True or if it's DATA: False
muonSelection  = "Glb"    # Single muon selection: Glb(isGlobal), GlbTrk(isGlobal&&isTracker), Trk(isTracker) are availale

triggerList    = {
    # Double Muon Trigger List
    'DoubleMuonTrigger' : cms.vstring(
        "HLT_HIL1DoubleMu0_v1",
        "HLT_HIL2DoubleMu0_v1",
        "HLT_HIL3DoubleMu0_v2",
        "HLT_L1DoubleMu0_v1",
        ),
    # Double Muon Filter List
    'DoubleMuonFilter'  : cms.vstring(
        "hltL1fL1sDoubleMu0L1Filtered0",
        "hltL2fL1sDoubleMu0L1f0L2Filtered0",
        "hltL3fL1sDoubleMu0L1f0L2f0L3Filtered0",
        "hltDoubleMu0L1Filtered",
        ),
    # Single Muon Trigger List
    'SingleMuonTrigger' : cms.vstring(
        "HLT_L1SingleMuOpen_v1",
        "HLT_L1SingleMu3_v1",
        "HLT_L1SingleMu5_v1",
        "HLT_L1SingleMu7_v1",
        "HLT_HIL2Mu7_v1",
        "HLT_HIL3Mu3_v1",
        "HLT_HIL3Mu5_v1",
        "HLT_HIL3Mu7_v2",
        ),
    # Single Muon Filter List
    'SingleMuonFilter'  : cms.vstring(
        "hltL1MuOpenL1Filtered0",
        "hltL1fL1sMu3L1Filtered0",
        "hltL1fL1sMu5L1Filtered0",
        "hltL1fL1sMu7L1Filtered0",
        "hltL2fL1sSingleMu3OR5L1f0L2Filtered7",
        "hltL3fL1sSingleMu3L1f0L2f0L3Filtered3",
        "hltL3fL1sSingleMu3OR5L1f0L2f0L3Filtered5",
        "hltL3fL1sSingleMu3OR5L1f0L2f0L3Filtered7",
        )
    }

if isMC:
    globalTag = '101X_upgrade2018_realistic_v8'
else:
    globalTag = '100X_dataRun2_v1'

#----------------------------------------------------------------------------


# set up process
process = cms.Process("TriggerAnalysis")

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# Input and Output File Names
options.outputFile = "Data_XeXe_OniaForest.root"
options.secondaryOutputFile = "Jpsi_DataSet.root"
options.inputFiles = [
#'/store/user/bdiab/SingleMuPt_0_30/SingleMuPt_0_30_CMSSW_10_1_7_RECO/180710_155328/0000/step2_RAW2DIGI_L1Reco_RECO_99.root'
'file:step3_XeXe_v7.root'
]
#options.inputFiles = '/store/user/bdiab/SingleMuPt_0_30/SingleMuPt_0_30_CMSSW_10_1_7_RECO/180710_155328/0000/step2_RAW2DIGI_L1Reco_RECO_99.root'
#options.inputFiles = '/store/data/Run2018B/DoubleMuon/AOD/PromptReco-v2/000/318/953/00000/6C6B2831-D97E-E811-8643-FA163E56EEC9.root'
options.maxEvents = -1 # -1 means all events

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
process.hltBitAna = cms.EndPath(process.hltbitanalysis)
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
process.hltObjectAna = cms.EndPath(process.hltobject)

#----------------------------------------------------------------------------

#Options:
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring( options.inputFiles ),
                            )
process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string( options.outputFile )
                                   )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.schedule  = cms.Schedule( process.oniaTreeAna , process.hltBitAna , process.hltObjectAna )
