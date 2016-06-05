import os
from functions import *

mass="CH100/"
os.system("mkdir " + mass)
os.system("ls /data4/DATA/AOD/CH/step2/" + mass +" > " + mass+"/samplelog.txt")

fread = open(mass+'samplelog.txt', 'r')

samplecounter=0
for samples in fread:
    if ".root" in samples:
        samplecounter+=1
        if not os.path.exists("/data4/DATA/LQNtuples/LQNtuples_5_3_14_snu27_PF2PAT/SIGNAL/"+mass ):
            os.system("mkdir " + "/data4/DATA/LQNtuples/LQNtuples_5_3_14_snu27_PF2PAT/SIGNAL/"+mass )
        configfile=open(mass+"CH_submit_mc"+ str(samplecounter) +"_53X.py","w")
        configfile.write(makeConfig(mass+"samplelog.txt", mass, "/data4/DATA/LQNtuples/LQNtuples_5_3_14_snu27_PF2PAT/SIGNAL/"+mass , samplecounter ))
        configfile.close()
        runcommand = "qsub -V " + mass+"CH_submit_mc"+ str(samplecounter) +"_53X.py"
        os.system(runcommand)
