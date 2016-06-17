import os
from functions import *

samples_create = [ "WZ", "ZZ","WW"]#, "TTZ","TTH","TTW"]
#samples_create= ["DY10to50" , "DY50plus", "WJets", "Wbb", "Ztobb"]

#samples_create= ["st_sch", "stbar_sch" , "st_tch", "stbar_tch", "st_tW" ,"stbar_tW"]
samples_create = ["Wbb"]

### name of dir in kisti (part of name)
tag="CRAB3_LQ_CH_MC_2016_June"

list_to_submit=samples_create

### version of local submittion
version="v2"
for i in list_to_submit:

    dataset=""
    dataset_postfix=""
    filename = 'datasetlist.txt'
    for line in open(filename, 'r'):
        if not line.startswith("#"):
            entries = line.split()
            if len(entries)==3:
                if i == entries[0]:
                    dataset=entries[1]
                    dataset_postfix=entries[2]
    
                    
    cfgfile='crab_submitMC_'+i+'.py'
    configfile=open(cfgfile,'w')
    configfile.write(make_submitfile(i,dataset,dataset_postfix,version, tag) )
    
    print "dataset=" + dataset           
    print "submit with : crab submit -c " + cfgfile
    print "check status with : crab status -d crab_projectsMC/LQ_CH_"+version+"_June2016_MC_"+i+"analysis"
