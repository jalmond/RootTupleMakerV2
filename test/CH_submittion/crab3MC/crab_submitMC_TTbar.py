from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'LQ_CH_v4_June2015_MC_analysis'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'rootTupleMaker_MC_2012_53X_chs_cfg.py'

config.Data.inputDataset = '/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM'   
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'CRAB3_LQ_CH_MC_2016_June_TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola'
config.Site.storageSite = 'T3_KR_KISTI'


