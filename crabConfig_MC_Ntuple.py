from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.requestName = 'NameOfYourProject'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'NameOfYourConfigFile'
config.JobType.outputFiles = ['OniaForest.root']
config.JobType.maxMemoryMB = 2500

config.section_('Data')
config.Data.inputDataset =''
config.Data.inputDBS = 'phys03'
config.Data.unitsPerJob = 1
config.Data.splitting = 'FileBased'

#config.Data.runRange = '304899-304907'
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/DCSOnly/json_DCSONLY.txt'


config.Data.outLFNDirBase = '/store/user/%s/NameOfYourProject' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName

config.section_('Site')
config.Data.ignoreLocality = False 
config.Site.storageSite = 'T2_CH_CERN'
