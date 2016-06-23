from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'LQ_CH_v2_June2016_MC_DYPowheg_1500analysis'
config.General.workArea = 'crab_projectsMC'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'rootTupleMaker_CRAB_MC_2012_53X_cfg.py'

config.Data.inputDataset = '/DYToMuMu_M-1500_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C2-v1/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'CRAB3_LQ_ISR_MC_2016_June_DYToMuMu_M-1500_CT10_TuneZ2star_8TeV-powheg-pythia6'
config.Site.blacklist = ['T2_ES_IFCA','T2_IT_Legnaro']
config.Site.storageSite = 'T3_KR_KISTI'
