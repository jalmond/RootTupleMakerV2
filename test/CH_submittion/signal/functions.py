def makeConfig(samplelist, mass, outputdir, job):

    ntuplemaker=open("CH_submit_mc_53X.py","r")
    

    config=""
    fread = open(samplelist, 'r')
    counter = 0
    for samples in fread:
        counter +=1
        if counter == job:
            for line in ntuplemaker:
            
                if "fileName = cms.string" in line:
                    config+='fileName = cms.string("' + outputdir + '/ntuple' + str(job) +'.root"),'     
                elif "file:/afs/" in line:
                    config+= "                         fileNames=cms.untracked.vstring('file:/data4/DATA/AOD/CH/step2/" + mass + samples.strip() + "'      ),"
                else:
                    config+=line +"\n"
            
            return config
