from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

data_type="ttbar"
version="v6"
datasetname="TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola"
dataset_postfix='/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM'

config.General.requestName = 'LQ_CH_'+version+'_June2016_MC_'+data_type+'analysis'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'rootTupleMaker_MC_2012_53X_chs_cfg.py'

config.Data.inputDataset = '/'+datasetname+dataset_postfix
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'CRAB3_LQ_CH_MC_2016_nopdfv2_June_'+datasetname
config.Site.storageSite = 'T3_KR_KISTI'


