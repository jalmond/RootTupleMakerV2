def make_submitfile(jobid, dataset, dataset_postfix, version, tag):
    

    config=""
    config+="from CRABClient.UserUtilities import config, getUsernameFromSiteDB\n"
    config+="config = config()\n"
    config+="\n"
    config+="config.General.requestName = 'LQ_CH_"+version+"_June2016_MC_"+ jobid+"analysis'\n"
    config+="config.General.workArea = 'crab_projectsMC'\n"
    config+="config.General.transferOutputs = True\n"
    config+="config.General.transferLogs = False\n"
    config+="\n"
    config+="config.JobType.pluginName = 'Analysis'\n"
    config+="config.JobType.psetName = 'rootTupleMaker_MC_2012_53X_chs_cfg.py'\n"
    config+="\n"
    config+="config.Data.inputDataset = '/"+dataset+dataset_postfix+"'\n"
    config+="config.Data.inputDBS = 'global'\n"
    config+="config.Data.splitting = 'FileBased'\n"
    config+='config.Data.unitsPerJob = 10\n'
    config+="config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())\n"
    config+='config.Data.publication = True\n'
    config+="config.Data.outputDatasetTag = '"+tag+"_"+dataset+"'\n"
    config+="config.Site.storageSite = 'T3_KR_KISTI'\n"
    return config
