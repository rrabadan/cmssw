# Auto generated configuration file
# using: 
# Revision: 1.334 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: test -n 2 -s RAW2DIGI,RECO,DQM,ALCA:HcalCalMinBias,REPACK:DigiToSplitRawRepack --eventcontent RECO,REPACKRAW,DQM --datatier RECO,REPACKRAW,DQM --data --scenario HeavyIons --no_exec --conditions GR_R_39X_V6B::All --filein /store/hidata/HIRun2010/HICorePhysics/RAW/v1/000/151/353/8AA9D775-F6F1-DF11-8658-001D09F2AF96.root --fileout output_rereco_test.root
import FWCore.ParameterSet.Config as cms
import sys
process = cms.Process('RECO')
#print(sys.argv[1])
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('DQMOffline.Configuration.DQMOfflineHeavyIons_cff')
process.load('Configuration.StandardSequences.AlCaRecoStreamsHeavyIons_cff')
process.load('Configuration.StandardSequences.DigiToRaw_Repack_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoLocalTracker.SiStripZeroSuppression.SiStripBaselineAnalyzer_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
#    fileNames = cms.untracked.vstring("file:" + sys.argv[2])
    fileNames = cms.untracked.vstring(
'/store/hidata/HIRun2010/HICorePhysics/RAW/v1/000/151/077/0475D7D4-3FF0-DF11-B171-003048F11C5C.root',
'/store/hidata/HIRun2010/HICorePhysics/RAW/v1/000/151/077/2C11C83E-35F0-DF11-B42E-0030486780EC.root',
'/store/hidata/HIRun2010/HICorePhysics/RAW/v1/000/151/077/00CC1018-25F0-DF11-AB5F-001D09F23D1D.root',
'/store/hidata/HIRun2010/HICorePhysics/RAW/v1/000/151/077/022817E8-27F0-DF11-B132-0030487C2B86.root',
'/store/hidata/HIRun2010/HICorePhysics/RAW/v1/000/151/077/0268E489-2DF0-DF11-8903-0030487C8CB8.root',
'/store/hidata/HIRun2010/HICorePhysics/RAW/v1/000/151/077/0296E43C-60F0-DF11-8CF8-003048D2C16E.root',
'/store/hidata/HIRun2010/HICorePhysics/RAW/v1/000/151/077/02A2D6F0-35F0-DF11-86EE-001D09F28EC1.root',
'/store/hidata/HIRun2010/HICorePhysics/RAW/v1/000/151/077/02A6F4BE-31F0-DF11-B292-001D09F26C5C.root',
'/store/hidata/HIRun2010/HICorePhysics/RAW/v1/000/151/077/040B4769-64F0-DF11-B328-001D09F2932B.root',
'/store/hidata/HIRun2010/HICorePhysics/RAW/v1/000/151/077/040C1ED3-2CF0-DF11-AD9F-003048F11942.root',
'/store/hidata/HIRun2010/HICorePhysics/RAW/v1/000/151/077/042CD873-39F0-DF11-BE6E-001617C3B69C.root',
'/store/hidata/HIRun2010/HICorePhysics/RAW/v1/000/151/077/04333B77-5DF0-DF11-9C29-003048F118AC.root'
#'/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/532/0E7F8B9C-190D-E111-A451-BCAEC5329702.root',
#'/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/532/147A0425-210D-E111-9341-BCAEC5329715.root',
#'/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/532/26FBC3D0-1D0D-E111-8848-BCAEC518FF68.root',
#'/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/532/78D524F2-230D-E111-9FE9-BCAEC518FF44.root',
#'/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/532/98C61024-1D0D-E111-9F44-003048F11114.root',
#'/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/532/AE499CFE-130D-E111-9C0A-BCAEC532971E.root',
#'/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/532/D21E0979-100D-E111-AEB7-0025901D62A6.root',
#'/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/532/E26095D0-160D-E111-8E0C-BCAEC5364CED.root',
#'/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/532/EA961950-0C0D-E111-A698-0025901D5D78.root'
#'/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/530/A8D6061E-030D-E111-A482-BCAEC532971A.root',
#'/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/181/530/5433308B-040D-E111-8004-003048F118C4.root',
#'/store/hidata/HIRun2011/HIHighPt/RAW/v1/000/182/838/9429025A-F51C-E111-B165-0025901D624E.root'
#'/store/hidata/HIRun2010/HIAllPhysics/RAW/v1/000/152/698/CE0B4D68-D9FB-DF11-87D0-001D09F29619.root'
#'/store/hidata/HIRun2010/HIAllPhysics/RAW/v1/000/152/698/26852261-DAFB-DF11-8B44-001D09F232B9.root'             
#'file:SD_Cen0_2p5_119_1_D2P.root'
#'/store/hidata/HIRun2010/HIAllPhysics/RAW/v1/000/152/698/00EF189B-CAFB-DF11-B3E7-003048F024FE.root'
#'/store/hidata/HIRun2010/HIAllPhysics/RAW/v1/000/152/698/FEE8B867-D9FB-DF11-8648-001D09F28F0C.root',
#'/store/hidata/HIRun2010/HIAllPhysics/RAW/v1/000/152/698/FEA10166-D9FB-DF11-A90C-0019B9F72F97.root',
#'/store/hidata/HIRun2010/HIAllPhysics/RAW/v1/000/152/698/CAF8525F-DAFB-DF11-93B4-0019B9F70607.root'
)
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.334 $'),
    annotation = cms.untracked.string('test nevts:2'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RECOoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
   # outputCommands = process.RECOEventContent.outputCommands,
    fileName = cms.untracked.string('output_rereco_test_1_DF_12_10_1_ME.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('RECO')
    )
)
        
#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string('StripHisto181969_300ev_DF_12_6_1.root')
#)

# Additional output definition

# Other statements
#process.GlobalTag.globaltag = 'GR_P_V43D::All'
process.GlobalTag.globaltag = 'GR_P_V46::All'

#process.SiStripBaselineAnalyzer.srcProcessedRawDigi =  cms.InputTag('siStripVRDigis','VirginRaw')
#process.SiStripBaselineAnalyzer.nModuletoDisplay = cms.uint32(10000)
#process.SiStripBaselineAnalyzer.plotClusters = cms.bool(True)
#process.SiStripBaselineAnalyzer.plotBaseline = cms.bool(False)
#process.SiStripBaselineAnalyzer.plotBaselinePoints = cms.bool(False)
#process.SiStripBaselineAnalyzer.plotRawDigi	= cms.bool(True)
#process.SiStripBaselineAnalyzer.plotAPVCM = cms.bool(False)
#process.SiStripBaselineAnalyzer.plotPedestals = cms.bool(True)
#process.siStripZeroSuppression.Algorithms.APVInspectMode = cms.string("BaselineFollower")
process.siStripZeroSuppression.Algorithms.APVInspectMode = cms.string("ForceAllInspect")
#process.siStripZeroSuppression.Algorithms.APVRestoreMode = cms.string("BaselineFollower")
process.siStripZeroSuppression.Algorithms.APVRestoreMode = cms.string("DerivativeFollower")
process.siStripZeroSuppression.Algorithms.lastGradient = cms.int32(10)

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstructionHeavyIons)
#process.reconstruction_step = cms.Path(process.siStripZeroSuppression)
#process.analyzer_step = cms.Path(process.SiStripBaselineAnalyzer)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOoutput_step = cms.EndPath(process.RECOoutput)




# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step, process.digi2repack_step, process.endjob_step, process.REPACKRAWoutput_step,process.analyzer_step)#,process.RECOoutput_step, process.analyzer_step,)
#
#process.schedule = cms.Schedule(process.raw2digi_step, process.reconstruction_step, process.endjob_step,process.RECOoutput_step, process.analyzer_step)

process.schedule = cms.Schedule(process.raw2digi_step, process.reconstruction_step, process.endjob_step,process.RECOoutput_step)
