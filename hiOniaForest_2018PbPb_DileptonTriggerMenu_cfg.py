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
			"HLT_HIL1DoubleMuOpen_OS_Centrality_40_100_v1", 
			"HLT_HIL1DoubleMuOpen_Centrality_50_100_v1", 
			"HLT_HIL1DoubleMu10_v1", 
			"HLT_HIL2_L1DoubleMu10_v1", 
			"HLT_HIL3_L1DoubleMu10_v1", 
			"HLT_HIL2DoubleMuOpen_v1", 
			"HLT_HIL3DoubleMuOpen_v1", 
			"HLT_HIL3DoubleMuOpen_M60120_v1", 
			"HLT_HIL3DoubleMuOpen_JpsiPsi_v1", 
			"HLT_HIL3DoubleMuOpen_Upsi_v1", 
			"HLT_HIL3Mu0_L2Mu0_v1", 
			"HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1",
			"HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v1",
			"HLT_HIL3Mu3_L1TripleMuOpen_v1"
),
		# Double Muon Filter List
		'DoubleMuonFilter'  : cms.vstring(
			"hltL1fL1sL1DoubleMuOpenL1Filtered0",
			"hltL1fL1sL1DoubleMuOpenOSCentrality40100L1Filtered0",
			"hltL1fL1sL1DoubleMuOpenCentrality50100L1Filtered0",
			"hltL1fL1sL1DoubleMu10L1Filtered0",
			"hltL2fL1sL1DoubleMu10L1f0L2Filtered0",
			"hltDoubleMuOpenL1DoubleMu10Filtered",
			"hltL2fL1sL1DoubleMuOpenL1f0L2Filtered0",
			"hltL3fL1DoubleMuOpenL3Filtered0",
			"hltL3fL1DoubleMuOpenL3FilteredM60120",
			"hltL3fL1DoubleMuOpenL3FilteredPsi",
			"hltL3fL1DoubleMuOpenL3FilteredUpsi",
			"hltL3f0L3Mu0L2Mu0Filtered0",
			"hltL3f0L3Mu0NHitQ10L2Mu0FilteredM1to5",
			"hltL3f0L3Mu2p5NHitQ10L2Mu2FilteredM7toinf",
			"hltL3fL1sL1DoubleMuOpenL1fN3L2f0L3Filtered3"
			),
		# Single Muon Trigger List
	'SingleMuonTrigger' : cms.vstring(
        "HLT_HIL1MuOpen_Centrality_70_100_v1",
        "HLT_HIL1MuOpen_Centrality_80_100_v1",
        "HLT_HIL2Mu3_NHitQ15_v1",
        "HLT_HIL2Mu5_NHitQ15_v1",
        "HLT_HIL2Mu7_NHitQ15_v1",
        "HLT_HIL3Mu3_NHitQ10_v1",
        "HLT_HIL3Mu5_NHitQ10_v1",
        "HLT_HIL3Mu7_NHitQ10_v1",
        "HLT_HIL3Mu12_v1",
        "HLT_HIL3Mu15_v1",
        "HLT_HIL3Mu20_v1",
			),
	# Single Muon Filter List
	'SingleMuonFilter'  : cms.vstring(
        "hltL1fL1sL1MuOpenCentrality70100L1Filtered0",
        "hltL1fL1sL1MuOpenCentrality80100L1Filtered0",
        "hltL2fL1sMuOpenL1f0L2Filtered3NHitQ15",
        "hltL2fL1sMuOpenL1f0L2Filtered5NHitQ15",
        "hltL2fL1sMuOpenL1f0L2Filtered7NHitQ15",
        "hltL3fL1sL1SingleMuOpenL1f0L2f0L3Filtered3NHitQ10",
        "hltL3fL1sL1SingleMuOpenL1f0L2f0L3Filtered5NHitQ10",
        "hltL3fL1sL1SingleMuOpenL1f0L2f0L3Filtered7NHitQ10",
        "hltL3fL1sL1SingleMuOpenL1f0L2f0L3Filtered12",
        "hltL3fL1sL1SingleMuOpenL1f0L2f0L3Filtered15",
        "hltL3fL1sL1SingleMuOpenL1f0L2f0L3Filtered20",
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
options.inputFiles = ['/store/group/phys_heavyions/jaebeom/JpsiMM_0_15/ReEmul_HLTMenuV35_JpsiLP_EALB/181023_181845/0009/step3_mc_9652.root']
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
