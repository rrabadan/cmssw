import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing


#----------------------------------------------------------------------------

# Setup Settings for ONIA FOREST:

HLTProcess     = "MyWonderfulHLT"    # Name of HLT process 
isMC           = True # if input is MONTECARLO: True or if it's DATA: False
muonSelection  = "Glb"    # Single muon selection: Glb(isGlobal), GlbTrk(isGlobal&&isTracker), Trk(isTracker) are availale

triggerList    = {
    # Double Muon Trigger List
    'DoubleMuonTrigger' : cms.vstring(
       "HLT_HIL1DoubleMuOpen_v1",
       "HLT_HIL1DoubleMuOpen_Centrality_30_100_v1",
       "HLT_HIL1DoubleMuOpen_OS_Centrality_30_100_v1",
       "HLT_HIL1DoubleMuOpen_Centrality_40_100_v1",
       "HLT_HIL1DoubleMuOpen_OS_Centrality_40_100_v1",
       "HLT_HIL1DoubleMuOpen_Centrality_50_100_v1",
       "HLT_HIL3Mu0_L2Mu0_v1",
        ),
    # Double Muon Filter List
    'DoubleMuonFilter'  : cms.vstring(
        "hltL1fL1sDoubleMuOpenL1Filtered0",
        "hltL1fL1sL1DoubleMuOpenCentrality30100L1Filtered0",
        "hltL1fL1sL1DoubleMuOpenOSCentrality30100L1Filtered0",
        "hltL1fL1sL1DoubleMuOpenCentrality40100L1Filtered0",
        "hltL1fL1sL1DoubleMuOpenOSCentrality40100L1Filtered0",
        "hltL1fL1sL1DoubleMuOpenCentrality50100L1Filtered0",
        "hltL3f0L3Mu0L2Mu0Filtered0",
        ),
    # Single Muon Trigger List
    'SingleMuonTrigger' : cms.vstring(
        "HLT_HIL1MuOpen_Centrality_70_100_v1",
        "HLT_HIL1MuOpen_Centrality_80_100_v1",
        "HLT_HIL1Mu3_Centrality_70_100_v1",
        "HLT_HIL1Mu3_Centrality_80_100_v1",
        ),
    # Single Muon Filter List
    'SingleMuonFilter'  : cms.vstring(
        "hltL1fL1sL1MuOpenCentrality70100L1Filtered0",
        "hltL1fL1sL1MuOpenCentrality80100L1Filtered0",
        "hltL1fL1sL1Mu3Centrality70100L1Filtered0",
        "hltL1fL1sL1Mu3Centrality80100L1Filtered0",
        )
    }

if isMC:
  #globalTag = 'auto:run2_mc_GRun'
  globalTag = '103X_upgrade2018_realistic_HI_v6'
else:
    globalTag = '101X_dataRun2_HLT_frozen_v6'

#----------------------------------------------------------------------------


# set up process
process = cms.Process("TriggerAnalysis")

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# Input and Output File Names
options.outputFile = "OniaForest.root"
options.secondaryOutputFile = "Jpsi_DataSet.root"
options.inputFiles = [
#'/store/user/bdiab/SingleMuPt_0_30/SingleMuPt_0_30_CMSSW_10_1_7_RECO/180710_155328/0000/step2_RAW2DIGI_L1Reco_RECO_99.root'
'/store/group/phys_heavyions/jaebeom/JpsiMM_0_30_emb_Gen/ReEmul_HLTMenuV27_EmbJpsi/181017_212930/0000/step3_mc_131.root'
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
#process.schedule  = cms.Schedule( process.oniaTreeAna , process.hltBitAna , process.hltObjectAna )
process.schedule  = cms.Schedule( process.oniaTreeAna )
