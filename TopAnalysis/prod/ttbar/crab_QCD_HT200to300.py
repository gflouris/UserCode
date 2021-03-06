from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()
config.General.requestName = 'QCD_HT200to300'
config.General.transferOutputs = True
config.General.transferLogs = False
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'flat-QCD-cfg.py'
config.JobType.maxJobRuntimeMin = 2800
#config.JobType.inputFiles = ['Summer15_50nsV2_MC.db']
config.Data.inputDataset = '/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/group/cmst3/user/kkousour/ttbar/'
config.Data.publication = False
config.Site.storageSite = 'T2_CH_CERN'
