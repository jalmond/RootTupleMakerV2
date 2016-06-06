import os
from functions import *

VV = [ "WZ"]#,"ZZ","WW", "TTZ","TTH","TTW"]

### name of dir in kisti (part of name)
tag="CRAB3_LQ_CH_MC_2016_June"

list_to_submit=VV

### version of local submittion
version="v1"
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
    
