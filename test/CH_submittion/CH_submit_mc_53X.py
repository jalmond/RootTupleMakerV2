##############################################################################################################
#### configuration
##############################################################################################################
print '*' * 60   
print "Configuration"
print '*' * 60   
##### specify if is Data
isData=False
runOnMC= not isData
# list of input files                                                                                                                                           
inputFiles = ['file:/afs/cern.ch/work/j/jalmond/FE4C2F81-D0E1-E111-9080-0030487E0A2D.root'] # overwritten, if "useRelVals" is 'True'                           
maxEvents = 20 # reduce for testing, should be -1
skipEvents=0

outputFile = 'CH_submit_mc_53X.root'
if isData:
  outputFile = 'CH_submit_data_53X.root'

##############################################################################################################                                                             
##### PRINT CONGURATION

print '=' * 60  
if isData:
  print "Job is configured to run on data"
else:
  print"Job is configured to run on MC"
print '=' * 60  

#### Load job settings + cuts
from Leptoquarks.RootTupleMakerV2.patRefSel_refMuJets import *


# output file

GLOBALTAG='START53_V27::All'

print '=' * 40  
print "Using Global Tag: "+GLOBALTAG


# maximum number of events
print "Running on "+str(maxEvents)+" events"

### Particle flow                                                                                                                                                      
postfix = "PF"

typeIMetCorrections = True

# event frequency of Fwk report
fwkReportEvery = 1000

# switch for 'TrigReport'/'TimeReport' at job end
wantSummary = True

################################################################################
###                              End of constants                            ###
###                                                                          ###
### ======================================================================== ###
###
### Basic PAT/JOB configuration
###

from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.GlobalTag.globaltag = GLOBALTAG
# Input files                                                                                                                                                             
process.maxEvents.input = maxEvents
#process.source.fileNames = inputFiles

#### Load module for good vertex 
process.load( "Leptoquarks.RootTupleMakerV2.patRefSel_goodVertex_cfi" )
########################################################################
### Check JECs
########################################################################
jecLevels = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
if not runOnMC:
  jecLevels = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])

print "JEC = " + str(jecLevels)
print '=' * 40


### Add in type - corrections
##____________________________________________________________________________||

print '-' * 40
process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")
print "Changing producePatPFMETCorrections: Adding type0PFMEtCorrection + patPFMETtype0Corr to sequence. This is needed to apply type 0 corrections "
print '-' * 40

process.producePatPFMETCorrections.replace(
  process.pfCandMETcorr,
  process.type0PFMEtCorrection *
  process.patPFMETtype0Corr *
  process.pfCandMETcorr 
  )

print '*' * 60
print "Setting up usePF2PAT"
print '*' * 60

from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
usePF2PAT( process
         , runPF2PAT           = True
         , runOnMC             = runOnMC
         , jetAlgo             = jetAlgo
         , postfix             = postfix
         , jetCorrections      = ( 
                                 jecLevels
                                 )
         , typeIMetCorrections = typeIMetCorrections
         , pvCollection        = cms.InputTag( pfVertices )
         )
print "\n" *2
print "usePF2PAT( process , runPF2PAT           = True         , runOnMC             = runOnMC         , jetAlgo             ="+ jetAlgo +"        , postfix             = "+postfix +"        , jetCorrections      = ( "+                             str(jecLevels)  +"                               )         , typeIMetCorrections = "+str(typeIMetCorrections)      +"   , pvCollection        = cms.InputTag( "+pfVertices+" )         )"
print "\n" *2
print '*' * 60
print "Finished setting up usePF2PAT"
print '*' * 60

postfixLoose='PFLoose'
from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
usePF2PAT( process
         , runPF2PAT           = True
         , runOnMC             = runOnMC
         , jetAlgo             = jetAlgo
         , postfix             = postfixLoose
         , jetCorrections      = (
                                 jecLevels
                                 )
         , typeIMetCorrections = typeIMetCorrections
         , pvCollection        = cms.InputTag( pfVertices )
         )
print "\n" *2
print "usePF2PAT( process , runPF2PAT           = True         , runOnMC             = runOnMC         , jetAlgo             ="+ jetAlgo +"        , postfix             = "+postfixLoose +"        , jetCorrections      = ( "+                             str(jecLevels)  +"                               )         , typeIMetCorrections = "+str(typeIMetCorrections)      +"   , pvCollection        = cms.InputTag( "+pfVertices+" )         )"
print "\n" *2
print '*' * 60
print "Finished setting up usePF2PATLoose"
print '*' * 60



## pile up corrections
from CommonTools.ParticleFlow.Tools.enablePileUpCorrection import *
#enablePileUpCorrectionInPF2PAT( process, postfix, sequence = "patPF2PATSequence"+postfix)
enablePileUpCorrection( process, postfix, sequence = "patPF2PATSequence"+postfix)

print '=' * 60
print "Ran enablePileUpCorrection: 'CommonTools.ParticleFlow.Tools.enablePileUpCorrection'"
print '=' * 60
################################################
################################################   ADD TOP PROJECTION CUTS TO PF Objects
################################################

print '=' * 60
print "----- pfMuons: setting cuts for top projection: Warning need to loosen isolation in PFLoose collection for analysis"
print "+++++++ "
print "process.pfSelectedMuonsPF.cut = 'abs(eta)<2.5 && pt>10.'"
print "process.pfIsolatedMuonsPF.doDeltaBetaCorrection = True"
print "process.pfIsolatedMuonsPF.deltaBetaFactor = -0.5"
print "process.pfIsolatedMuonsPF.isolationCut = 0.20"
print '=' * 60

#### for top projection
getattr( process, 'pfSelectedMuons' + postfix ).cut = 'abs(eta)<2.5 && pt>10.'
getattr( process, 'pfIsolatedMuons' + postfix ).doDeltaBetaCorrection = True
getattr( process, 'pfIsolatedMuons' + postfix ).deltaBetaFactor = -0.5
getattr( process, 'pfIsolatedMuons' + postfix ).isolationCut = 0.20
### Loose
print "process.pfSelectedMuonsPFLoose.cut = 'abs(eta)<2.5 && pt>10.'"
print "process.pfIsolatedMuonsPFLoose.doDeltaBetaCorrection = True"
print "process.pfIsolatedMuonsPFLoose.deltaBetaFactor = -0.5"
print "process.pfIsolatedMuonsPFLoose.isolationCut = 0.50"
getattr( process, 'pfSelectedMuons' + postfixLoose ).cut = 'abs(eta)<2.5 && pt>10.'
getattr( process, 'pfIsolatedMuons' + postfixLoose ).doDeltaBetaCorrection = True
getattr( process, 'pfIsolatedMuons' + postfixLoose ).deltaBetaFactor = -0.5
getattr( process, 'pfIsolatedMuons' + postfixLoose ).isolationCut = 0.50


print '=' * 60
print "----- pfElectrons: setting cuts for top projections"
print "WARNING: These cuts are too tight for analysis. Need to make new patElectronLoose collection with looser pt and isolation"
print '=' * 60

process.load( "EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi" )
process.load('EgammaAnalysis.ElectronTools.electronIsolatorFromEffectiveArea_cfi')

process.pfIdentifiedElectronsPF = cms.EDFilter("ElectronIDPFCandidateSelector",
                                              recoGsfElectrons = cms.InputTag("gsfElectrons"),
                                              electronIdMap = cms.InputTag("mvaTrigV0"),
                                              electronIdCut = cms.double(0.0),
                                              src = cms.InputTag("pfElectronsFromVertex"+postfix))

getattr( process, 'pfSelectedElectrons' + postfix).src = 'pfIdentifiedElectrons'+postfix
getattr( process, 'pfSelectedElectrons' + postfix).cut = 'abs(eta)<2.5 && pt>20. && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits<2'
print "pfSelectedElectronsPF.cut = 'abs(eta)<2.5 && pt>20. && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits<2'"
getattr( process, 'pfIsolatedElectrons' + postfix).isolationValueMapsCharged  = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFId"+postfix))
getattr( process, 'pfIsolatedElectrons' + postfix).isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFId"+postfix),
                                                                      cms.InputTag("elPFIsoValueGamma03PFId"+postfix))

process.elPFIsoValueEA03.pfElectrons = cms.InputTag('pfSelectedElectrons' + postfix)
getattr( process, 'pfIsolatedElectrons' + postfix).deltaBetaIsolationValueMap = 'elPFIsoValueEA03' # new Effective Area edm::ValueMap 
getattr( process, 'pfIsolatedElectrons' + postfix).doDeltaBetaCorrection = True                    # not really a 'deltaBeta' correction, but it serves  
getattr( process, 'pfIsolatedElectrons' + postfix).deltaBetaFactor = -1.0
getattr( process, 'pfIsolatedElectrons' + postfix).isolationCut = 0.15
print "pfIsolatedElectronsPF.isolationCut = 0.15"

process.patPF2PATSequencePF.replace( process.pfSelectedElectronsPF,
                                   process.pfIdentifiedElectronsPF + process.pfSelectedElectronsPF  + process.elPFIsoValueEA03)
#### Loose

process.pfIdentifiedElectronsPFLoose = cms.EDFilter("ElectronIDPFCandidateSelector",
                                              recoGsfElectrons = cms.InputTag("gsfElectrons"),
                                              electronIdMap = cms.InputTag("mvaTrigV0"),
                                              electronIdCut = cms.double(0.0),
                                              src = cms.InputTag("pfElectronsFromVertex"+postfixLoose))

getattr( process, 'pfSelectedElectrons' + postfixLoose ).src = 'pfIdentifiedElectrons'+postfixLoose
getattr( process, 'pfSelectedElectrons' + postfixLoose ).cut = 'abs(eta)<2.5 && pt>10. && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits<2'

getattr( process, 'pfIsolatedElectrons' + postfixLoose ).isolationValueMapsCharged  = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFId"+postfixLoose))
getattr( process, 'pfIsolatedElectrons' + postfixLoose ).isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFId"+postfixLoose),
                                                                      cms.InputTag("elPFIsoValueGamma03PFId"+postfixLoose))

process.elPFIsoValueEA03Loose = process.elPFIsoValueEA03.clone()
process.elPFIsoValueEA03Loose.pfElectrons = cms.InputTag('pfSelectedElectrons' + postfixLoose)
getattr( process, 'pfIsolatedElectrons' + postfixLoose ).deltaBetaIsolationValueMap = 'elPFIsoValueEA03Loose' # new Effective Area edm::ValueMap    
getattr( process, 'pfIsolatedElectrons' + postfixLoose ).doDeltaBetaCorrection = True                    # not really a 'deltaBeta' correction, but it serves 
getattr( process, 'pfIsolatedElectrons' + postfixLoose ).deltaBetaFactor = -1.0
getattr( process, 'pfIsolatedElectrons' + postfixLoose ).isolationCut = 0.5
print "pfIsolatedElectronsPFLoose.isolationCut = 0.5"


process.patPF2PATSequencePFLoose.replace( process.pfSelectedElectronsPFLoose,
                                   process.pfIdentifiedElectronsPFLoose + process.pfSelectedElectronsPFLoose + process.elPFIsoValueEA03Loose )


print '=' * 60
print "----- PATElectron isolation: Settinig 03 cone isolation as default"
print '=' * 60
getattr( process, 'patElectrons' + postfix ).isolationValues.pfNeutralHadrons   = cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
getattr( process, 'patElectrons' + postfix ).isolationValues.pfChargedAll       = cms.InputTag( 'elPFIsoValueChargedAll03PFId' + postfix )
getattr( process, 'patElectrons' + postfix ).isolationValues.pfPUChargedHadrons = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
getattr( process, 'patElectrons' + postfix ).isolationValues.pfPhotons          = cms.InputTag( 'elPFIsoValueGamma03PFId' + postfix )
getattr( process, 'patElectrons' + postfix ).isolationValues.pfChargedHadrons   = cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix )
getattr( process, 'patElectrons' + postfix ).isolationValues.user =   cms.VInputTag(cms.InputTag("elPFIsoValueEA03"),"elPFIsoValueCharged04PFIdPF","elPFIsoValueNeutral04PFIdPF","elPFIsoValueGamma04PFIdPF","elPFIsoValuePU04PFIdPF","elPFIsoValueChargedAll04PFIdPF")

getattr( process, 'patElectrons' + postfixLoose ).isolationValues.pfNeutralHadrons   = cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfixLoose )
getattr( process, 'patElectrons' + postfixLoose ).isolationValues.pfChargedAll       = cms.InputTag( 'elPFIsoValueChargedAll03PFId' + postfixLoose )
getattr( process, 'patElectrons' + postfixLoose ).isolationValues.pfPUChargedHadrons = cms.InputTag( 'elPFIsoValuePU03PFId' + postfixLoose )
getattr( process, 'patElectrons' + postfixLoose ).isolationValues.pfPhotons          = cms.InputTag( 'elPFIsoValueGamma03PFId' + postfixLoose )
getattr( process, 'patElectrons' + postfixLoose ).isolationValues.pfChargedHadrons   = cms.InputTag( 'elPFIsoValueCharged03PFId' + postfixLoose )
getattr( process, 'patElectrons' + postfixLoose ).isolationValues.user =   cms.VInputTag(cms.InputTag("elPFIsoValueEA03Loose"),"elPFIsoValueCharged04PFIdPFLoose","elPFIsoValueNeutral04PFIdPFLoose","elPFIsoValueGamma04PFIdPFLoose","elPFIsoValuePU04PFIdPFLoose","elPFIsoValueChargedAll04PFIdPFLoose")


print '=' * 60
print "----- Set up electron PATConversionProducer + MVA electron ID"
print '=' * 60
process.patConversions = cms.EDProducer("PATConversionProducer",
    electronSource = cms.InputTag('selectedPatElectrons' + postfix)
)

# MVA electron ID                                                                                                                                                       
process.eidMVASequence = cms.Sequence(
  process.mvaTrigV0
+ process.mvaNonTrigV0
)

print '=' * 60
print "----- Turn on Top Projection"
print '=' * 60
# top projections in PF2PAT:
getattr(process,"pfNoPileUp"+postfix).enable = True
getattr(process,"pfNoMuon"+postfix).enable = True
getattr(process,"pfNoElectron"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True
# verbose flags for the PF2PAT modules
getattr(process,"pfNoMuon"+postfix).verbose = False

print 'getattr(process,"pfNoPileUp"+postfix).enable = True'
print 'getattr(process,"pfNoMuon"+postfix).enable = True'
print 'getattr(process,"pfNoElectron"+postfix).enable = True'
print 'getattr(process,"pfNoTau"+postfix).enable = False'
print 'getattr(process,"pfNoJet"+postfix).enable = True'




## LOADING Ntuple modules
process.load('Leptoquarks.RootTupleMakerV2.Ntuple_cff')

# Output ROOT file                                                                                                                                                        
process.TFileService = cms.Service("TFileService",
    fileName = cms.string( outputFile)
)

# load the PAT trigger Python tools                                                                                                                                            
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTriggerMatching( process, sequence = "patPF2PATSequence"+postfix , triggerMatchers = [
     'cleanElectronTriggerMatchHLTSingleElectron8',
        'cleanElectronTriggerMatchHLTSingleElectron17',
     'cleanElectronTriggerMatchHLTSingleElectronWP80',
        'cleanElectronTriggerMatchHLTDoubleElectron',
     'cleanMuonTriggerMatchHLTSingleMuon',
        'cleanMuonTriggerMatchHLTSingleMuon5',
        'cleanMuonTriggerMatchHLTSingleMuon8',
        'cleanMuonTriggerMatchHLTSingleMuon12',
        'cleanMuonTriggerMatchHLTSingleMuon17',
        'cleanMuonTriggerMatchHLTSingleMuon24',
        'cleanMuonTriggerMatchHLTDoubleMuon',
        'cleanMuonTriggerMatchHLTSingleIsoMuon'] )


from PhysicsTools.PatAlgos.tools.coreTools import *


##### LOAD predefined selections cuts
from Leptoquarks.RootTupleMakerV2.patRefSel_ch_cfi import *

# remove MC matching, object cleaning, objects etc.
if isData:
  runOnData( process
           , names = [ 'PFAll' ]
           , postfix = postfix
           )


# additional event content has to be (re-)added _after_ the call to 'removeCleaning()':
process.out.outputCommands += [ 'keep edmTriggerResults_*_*_*'
                                ,                                 'keep *_elPFIsoValueEA03_*_*'
                              , 'keep *_hltTriggerSummaryAOD_*_*'
                              # vertices and beam spot
                              , 'keep *_offlineBeamSpot_*_*'
                              , 'keep *_offlinePrimaryVertices*_*_*'
                              , 'keep *_goodOfflinePrimaryVertices*_*_*'
                              ]
if runOnMC:
  process.out.outputCommands += [ 'keep GenEventInfoProduct_*_*_*'
                                , 'keep recoGenParticles_*_*_*'
                                , 'keep *_addPileupInfo_*_*'
                                ]
###
### Additional configuration
###


### ELectrons (not really needed any more)
intermediatePatElectrons.src = cms.InputTag( 'selectedPatElectrons' + postfixLoose )
setattr( process, 'intermediatePatElectrons' + postfix, intermediatePatElectrons )
### Muons
process.out.outputCommands.append( 'keep *_pfMuonsFromVertex*_*_*' )

###
### Selection configuration
###

### Muons

# Use the good primary vertices everywhere. (This is done a few times), but no harm in taht
for imod in [process.patMuonsPF,
             process.patMuonsPFLoose,
             process.patElectronsPFLoose,
             process.patElectronsPF]:
    imod.pvSrc = "goodOfflinePrimaryVertices"
    imod.embedTrack = True
    imod.usePV = False

### Basic pt and eta cut for us in jet smearing tool 
getattr( process, 'selectedPatMuons' + postfix ).cut = LQmuonCut
getattr( process, 'selectedPatMuons' + postfixLoose ).cut = LQmuonCut

### Add user isolation (so we have both dr cones)
getattr( process, 'patMuons' + postfix ).user = cms.VInputTag("muPFIsoValueCharged03PF","muPFIsoValueNeutral03PF","muPFIsoValueGamma03PF","muPFIsoValuePU03PF","muPFIsoValueChargedAll03PF")
getattr( process, 'patMuons' + postfixLoose ).user = cms.VInputTag("muPFIsoValueCharged03PF","muPFIsoValueNeutral03PF","muPFIsoValueGamma03PF","muPFIsoValuePU03PF","muPFIsoValueChargedAll03PF")

### Jets Apply loose pt cuts on jets
getattr( process, 'selectedPatJets'  + postfix ).cut = LQJetCut

### Electrons
getattr( process, 'selectedPatElectrons' + postfix ).cut = LQelectronCut
getattr( process, 'selectedPatElectrons' + postfixLoose ).cut = LQelectronCut

# change PV in Jet corrector                                                                                                                                     
for module in [
               process.patJetCorrFactorsPF
               ]:
    module.primaryVertices = "goodOfflinePrimaryVertices"


# Options and Output Report                                                                                                                                      
process.options.wantSummary = wantSummary

import os

############## IMPORTANT ########################################                                                                                                
# If you run over many samples and you save the log, remember to reduce                                                                                          
# the size of the output by prescaling the report of the event number                                                                                            
process.MessageLogger.cerr.FwkReport.reportEvery = fwkReportEvery
process.MessageLogger.cerr.default.limit = 10
################################################################
process.source = cms.Source ("PoolSource",
                             fileNames=cms.untracked.vstring('file:/afs/cern.ch/work/j/jalmond/FE4C2F81-D0E1-E111-9080-0030487E0A2D.root'      ),
)

process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJetsForHEEPIsolation = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForHEEPIsolation.Rho_EtaMax = cms.double(2.5)

#----------------------------------------------------------------------------------------------------                                                            
# Add MET filters                                                                                                                                                
#----------------------------------------------------------------------------------------------------                                                            

process.load("Leptoquarks.RootTupleMakerV2.metFilters_cfi")


#----------------------------------------------------------------------------------------------------                                                            

# Modify cleanPatTaus (HPS Taus) - loosen up a bit                                                                                                               
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatAlgos/python/cleaningLayer1/tauCleaner_cfi.py?revision=1.11&view=markup                       
#----------------------------------------------------------------------------------------------------                                                            

getattr( process, 'selectedPatTaus' + postfix ).cut = cms.string("pt > 10. && tauID('decayModeFinding')> 0.5")
getattr( process, 'selectedPatTaus' + postfix ).preselection = cms.string(' tauID("decayModeFinding") > 0.5 ')
getattr( process, 'selectedPatTaus' + postfix ).finalCut     = cms.string(' pt > 15.0 & abs(eta) < 2.5      ')

#----------------------------------------------------------------------------------------------------                                                            
# Add Tau ID sources (HPS Taus)                                                                                                                                  
#----------------------------------------------------------------------------------------------------                                                            

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

#----------------------------------------------------------------------------------------------------                                                            
# Add the HEEP ID bit to the electrons                                                                                                                           
#----------------------------------------------------------------------------------------------------                                                            

from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
process.HEEPId = cms.EDProducer("HEEPIdValueMapProducer",
                                eleLabel = cms.InputTag("gsfElectrons"),
                                barrelCuts = cms.PSet(heepBarrelCuts),
                                endcapCuts = cms.PSet(heepEndcapCuts),
                                eleIsolEffectiveAreas = cms.PSet(heepEffectiveAreas),
                                eleRhoCorrLabel = cms.InputTag("kt6PFJetsForIsolation","rho"),
                                applyRhoCorrToEleIsol = cms.bool(True),
                                verticesLabel = cms.InputTag("offlinePrimaryVerticesWithBS"),
                                writeIdAsInt =cms.bool(True)
                                )
getattr( process, 'patElectrons' + postfix ).userData.userInts.src = cms.VInputTag('HEEPId')
getattr( process, 'patElectrons' + postfix ).electronIDSources = electronIDSources
getattr( process, 'patElectrons' + postfixLoose ).electronIDSources = electronIDSources
getattr( process, 'patElectrons' + postfixLoose ).userData.userInts.src = cms.VInputTag('HEEPId')

#----------------------------------------------------------------------------------------------------  
# Make analysisPatTaus and add them to the cleanPatCandidates sequence  
#----------------------------------------------------------------------------------------------------                                                            

process.analysisPatTaus = process.selectedPatTausPF.clone()
process.analysisPatTaus.preselection = cms.string(
    'tauID("decayModeFinding") > 0.5 &'
    ' tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 &'
    ' tauID("againstMuonLoose3") > 0.5 &'
    ' tauID("againstElectronLooseMVA3") > 0.5'
)
process.analysisPatTaus.finalCut = cms.string('pt > 20. & abs(eta) < 2.3')

process.patPF2PATSequencePF.replace( process.selectedPatTausPF, process.selectedPatTausPF +process.analysisPatTaus)



### setting up MET                                                                                                                                               
getattr(process,'patPFMet'+postfix).addGenMET = cms.bool(not isData)
#process.patPFMet.addGenMET = cms.bool(not isData)

#----------------------------------------------------------------------------------------------------                                                            
# Define the systematic shift correction                                                                                                                         
#----------------------------------------------------------------------------------------------------                                                            
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")
from PhysicsTools.PatAlgos.tools.helpers  import cloneProcessingSnippet
cloneProcessingSnippet(process, process.pfMEtSysShiftCorrSequence, postfix)
process.pfMEtSysShiftCorrPF.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc
process.pfMEtSysShiftCorrPF.src =  cms.InputTag('pfMet')

#----------------------------------------------------------------------------------------------------                                                            
###### MET                                                                                                                                                               
print '=' * 60
print "----- Add patPFMETtype0Corr to patType1CorrectedPFMet"
print '=' * 60

process.metUncertaintySequence = cms.Sequence()
#### patMETs changed to == type 1 corrected MET. 
getattr(process,'patMETs'+postfix).metSource = cms.InputTag("patType1CorrectedPFMet"+postfix)
getattr(process,'patType1CorrectedPFMet'+postfix).srcType1Corrections = cms.VInputTag(
    cms.InputTag("patPFJetMETtype1p2Corr"+postfix,"type1")
)


jetSmearFileName='PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'
jetSmearHistogram='pfJetResolutionMCtoDataCorrLUT'

import RecoMET.METProducers.METSigParams_cfi as jetResolutions
process.smearedPatPFMetSequencePF = cms.Sequence()
if not isData:
  SmearMET=True ### smear raw met
  process.smearedAnalysisPatJets = cms.EDProducer("SmearedPATJetProducer",
                                                  src = cms.InputTag('selectedPatJets'     + postfix),
                                                  dRmaxGenJetMatch = cms.string('TMath::Min(0.5, 0.1 + 0.3*TMath::Exp(-0.05*(genJetPt - 10.)))'),
                                                  sigmaMaxGenJetMatch = cms.double(5.),
                                                  inputFileName = cms.FileInPath(jetSmearFileName),
                                                  lutName = cms.string(jetSmearHistogram),
                                                  jetResolutions = jetResolutions.METSignificance_params,
                                                  skipJetSelection = cms.string(
      'jecSetsAvailable & abs(energy - correctedP4("Uncorrected").energy) > (5.*min(energy, correctedP4("Uncorrected").energy))'
      ),
                                                  skipRawJetPtThreshold = cms.double(10.), # GeV                                                                                       
                                                  skipCorrJetPtThreshold = cms.double(1.e-2)
                                                  )

  if SmearMET:
    setattr(process, "patPFMetForMEtUncertainty"+postfix, getattr(process, "patPFMet"+postfix).clone())
    process.smearedPatPFMetSequencePF += getattr(process, "patPFMetForMEtUncertainty"+postfix)
    setattr(process, "patPFMETcorrJetSmearing"+postfix, cms.EDProducer("ShiftedParticleMETcorrInputProducer",
                                                                       srcOriginal = cms.InputTag("selectedPatJets"+postfix),
                                                                       srcShifted = cms.InputTag("smearedAnalysisPatJets")
                                                                       ))
    process.smearedPatPFMetSequencePF += getattr(process, "patPFMETcorrJetSmearing"+postfix)
    getattr(process, "producePatPFMETCorrections"+postfix).replace(getattr(process, "patPFMet"+postfix), process.smearedAnalysisPatJets+ process.smearedPatPFMetSequencePF)
    setattr(process, "patPFMet"+postfix, getattr(process, "patType1CorrectedPFMet"+postfix).clone(
        src = cms.InputTag('patPFMetForMEtUncertainty'+postfix),
        srcType1Corrections = cms.VInputTag(
          cms.InputTag('patPFMETcorrJetSmearing'+postfix)
          )
        ))
    process.smearedPatPFMetSequencePF += getattr(process, "patPFMet"+postfix)

### FIX run metuncertaunty
### patType1CorrectedPFMetType1OnlyPF == patMETsPF
process.patType1CorrectedPFMetType1OnlyPF = process.patType1CorrectedPFMetPF.clone()
process.patType1CorrectedPFMetType1OnlyPF.srcType1Corrections = cms.VInputTag(
    cms.InputTag("patPFJetMETtype1p2CorrPF","type1"),
)

### patType1CorrectedPFMetType01OnlyPF = Type 0 + 1 
process.patType1CorrectedPFMetType01OnlyPF = process.patType1CorrectedPFMetPF.clone()
process.patType1CorrectedPFMetType01OnlyPF.srcType1Corrections = cms.VInputTag(
    cms.InputTag("patPFJetMETtype1p2CorrPF","type1"),
    cms.InputTag("patPFMETtype0CorrPF"),
    # cms.InputTag("pfMEtSysShiftCorr")                                                                                                                                                                          
)

### patType1CorrectedPFMetType01XYOnlyPF = Type 0 + 1   + phi corrected
process.patType1CorrectedPFMetType01XYOnlyPF = process.patType1CorrectedPFMetPF.clone()
process.patType1CorrectedPFMetType01XYOnlyPF.srcType1Corrections = cms.VInputTag(
    cms.InputTag("patPFJetMETtype1p2CorrPF","type1"),
    cms.InputTag("patPFMETtype0CorrPF"),
    cms.InputTag("pfMEtSysShiftCorrPF"), 
)

print '=' * 60
print "----- MET unclustered corrections: For systematics"
print '=' * 60

unclEnMETcorrections = [
  [ 'pfCandMETcorr'+postfix, [ '' ] ],
  [ 'patPFJetMETtype1p2Corr'+postfix, [ 'type2', 'offset' ] ],
  [ 'patPFJetMETtype2Corr'+postfix, [ 'type2' ] ],
  ]
unclEnMETcorrectionsUp = [ 
  cms.InputTag("patPFJetMETtype1p2CorrPF","type1"),
  cms.InputTag("patPFMETtype0CorrPF"),
  cms.InputTag("pfMEtSysShiftCorrPF"),]
unclEnMETcorrectionsDown = [ 
  cms.InputTag("patPFJetMETtype1p2CorrPF","type1"),
  cms.InputTag("patPFMETtype0CorrPF"),
  cms.InputTag("pfMEtSysShiftCorrPF"),]

varyByNsigmas=1.
basename=""
for srcUnclEnMETcorr in unclEnMETcorrections:

  moduleUnclEnMETcorrUp = cms.EDProducer("ShiftedMETcorrInputProducer",
                                         src = cms.VInputTag(
      [ cms.InputTag(srcUnclEnMETcorr[0], instanceLabel) for instanceLabel in srcUnclEnMETcorr[1] ]
      ),
                                         uncertainty = cms.double(0.10),
                                         shiftBy = cms.double(+1.*varyByNsigmas)
                                         )
  baseName = srcUnclEnMETcorr[0]
  
  if postfix != "":
    if baseName[-len(postfix):] == postfix:
      baseName = baseName[0:-len(postfix)]
    else:
      raise StandardError("Tried to remove postfix %s from label %s, but it wasn't there" % (postfix, baseName))
    moduleUnclEnMETcorrUpName = "%sUnclusteredEnUp" % baseName
    moduleUnclEnMETcorrUpName += postfix

    setattr(process, moduleUnclEnMETcorrUpName, moduleUnclEnMETcorrUp)

    process.metUncertaintySequence += getattr(process, moduleUnclEnMETcorrUpName)

    for instanceLabel in srcUnclEnMETcorr[1]:
      unclEnMETcorrectionsUp.extend([ cms.InputTag(moduleUnclEnMETcorrUpName, instanceLabel)
                                    for instanceLabel in srcUnclEnMETcorr[1] ] )
    moduleUnclEnMETcorrDown = moduleUnclEnMETcorrUp.clone(
      shiftBy = cms.double(-1.*varyByNsigmas)
      )
    moduleUnclEnMETcorrDownName = "%sUnclusteredEnDown" % baseName
    moduleUnclEnMETcorrDownName += postfix
    setattr(process, moduleUnclEnMETcorrDownName, moduleUnclEnMETcorrDown)
    process.metUncertaintySequence += getattr(process, moduleUnclEnMETcorrDownName)

    for instanceLabel in srcUnclEnMETcorr[1]:
      unclEnMETcorrectionsDown.extend([ cms.InputTag(moduleUnclEnMETcorrDownName, instanceLabel)
                                      for instanceLabel in srcUnclEnMETcorr[1] ] )


#### patType1CorrectedPFMetType01XYOnlyPFUnclustX == corrected Tpye01XY
process.patType1CorrectedPFMetType01XYOnlyPFUnclustup = process.patType1CorrectedPFMetType01XYOnlyPF.clone()
process.patType1CorrectedPFMetType01XYOnlyPFUnclustup.srcType1Corrections = cms.VInputTag(unclEnMETcorrectionsUp)
process.patType1CorrectedPFMetType01XYOnlyPFUnclustdown = process.patType1CorrectedPFMetType01XYOnlyPF.clone()
process.patType1CorrectedPFMetType01XYOnlyPFUnclustdown.srcType1Corrections = cms.VInputTag(unclEnMETcorrectionsDown)


from PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi import patMETs

### Create Raw MET collection
process.RawMET = patMETs.clone(
    metSource = cms.InputTag('patPFMetPF'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue')
)

#### Noew create Jet smeared and shifted collections
#jetSmearFileName='PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'
#jetSmearHistogram='pfJetResolutionMCtoDataCorrLUT'

#import RecoMET.METProducers.METSigParams_cfi as jetResolutions

#process.smearedAnalysisPatJets = cms.EDProducer("SmearedPATJetProducer",
#            src = cms.InputTag('selectedPatJets'     + postfix),
#            dRmaxGenJetMatch = cms.string('TMath::Min(0.5, 0.1 + 0.3*TMath::Exp(-0.05*(genJetPt - 10.)))'),
#            sigmaMaxGenJetMatch = cms.double(5.),
#            inputFileName = cms.FileInPath(jetSmearFileName),
#            lutName = cms.string(jetSmearHistogram),
#            jetResolutions = jetResolutions.METSignificance_params,
#            skipJetSelection = cms.string(
#        'jecSetsAvailable & abs(energy - correctedP4("Uncorrected").energy) > (5.*min(energy, correctedP4("Uncorrected").energy))'
#        ),
#                                         skipRawJetPtThreshold = cms.double(10.), # GeV                                                             
#                                         skipCorrJetPtThreshold = cms.double(1.e-2)
#                                               )

varyByNsigmas = 1.0
process.smearedJetsResDown = process.smearedAnalysisPatJets.clone(
  shiftBy = cms.double(-1.*varyByNsigmas)
  )
process.smearedJetsResUp = process.smearedAnalysisPatJets.clone(
  shiftBy = cms.double(+1.*varyByNsigmas)
  )

process.shiftedJetsEnUp = cms.EDProducer("ShiftedPATJetProducer",
                                         src = cms.InputTag('selectedPatJets'   + postfix),
                                         jetCorrPayloadName = cms.string('AK5PFchs'),
                                         jetCorrUncertaintyTag = cms.string("Uncertainty"),
                                         addResidualJES = cms.bool(True),
                                         jetCorrLabelUpToL3 = cms.string("ak5PFL1FastL2L3"),
            jetCorrLabelUpToL3Res = cms.string("ak5PFL1FastL2L3Residual"),                               
                                         shiftBy = cms.double(+1.*varyByNsigmas)
                                         )
process.shiftedJetsEnDown = process.shiftedJetsEnUp.clone(
  shiftBy = cms.double(-1.*varyByNsigmas)
  )

####### For systematics create shifted Lepton collections
### Muons
process.muonsEnUp = cms.EDProducer("ShiftedPATMuonProducer",
                           src =  cms.InputTag('selectedPatMuons' + postfixLoose),
                           uncertainty = cms.double(0.002),
                           shiftBy = cms.double(+1.*varyByNsigmas)
                           )
process.muonsEnDown = cms.EDProducer("ShiftedPATMuonProducer",
                           src =  cms.InputTag('selectedPatMuons' + postfixLoose),
                           uncertainty = cms.double(0.002),
                           shiftBy = cms.double(-1.*varyByNsigmas)
                           )

### Electrons
process.electronsEnUp = cms.EDProducer("ShiftedPATElectronProducer",
                src = cms.InputTag('selectedPatElectrons' + postfixLoose),
                binning = cms.VPSet(
                    cms.PSet(
                        binSelection = cms.string('isEB'),
                        binUncertainty = cms.double(0.006)
                    ),
                    cms.PSet(
                        binSelection = cms.string('!isEB'),
                        binUncertainty = cms.double(0.015)
                    ),
                ),      
                shiftBy = cms.double(+1.*varyByNsigmas)
            )

process.electronsEnDown = cms.EDProducer("ShiftedPATElectronProducer",
                src = cms.InputTag('selectedPatElectrons' + postfixLoose),
                binning = cms.VPSet(
                    cms.PSet(
                        binSelection = cms.string('isEB'),
                        binUncertainty = cms.double(0.006)
                    ),
                    cms.PSet(
                        binSelection = cms.string('!isEB'),
                        binUncertainty = cms.double(0.015)
                    ),
                ),      
                shiftBy = cms.double(-1.*varyByNsigmas)
            )

### Muons
process.rootTupleMuons.InputTag = cms.InputTag('selectedPatMuons'+ postfix)
process.rootTupleMuons.InputTagEnUp  = cms.InputTag('selectedPatMuons'+ postfix)
process.rootTupleMuons.InputTagEnDown  = cms.InputTag('selectedPatMuons'+ postfix)
### Loose Muons
process.rootTupleMuonsLoose.InputTag = cms.InputTag('selectedPatMuons'+ postfixLoose)
process.rootTupleMuonsLoose.InputTagEnUp  = cms.InputTag('muonsEnUp')
process.rootTupleMuonsLoose.InputTagEnDown  = cms.InputTag('muonsEnDown')

### Electrons
process.rootTupleElectrons.InputTag = cms.InputTag('selectedPatElectrons' + postfix)
process.rootTupleElectrons.InputTagEnUp  =  cms.InputTag('selectedPatElectrons' + postfix)
process.rootTupleElectrons.InputTagEnDown =  cms.InputTag('selectedPatElectrons' + postfix)
### Loose Electtrons
process.rootTupleElectronsLoose.InputTag = cms.InputTag('selectedPatElectrons' + postfixLoose)
process.rootTupleElectronsLoose.InputTagEnUp  =  cms.InputTag('electronsEnUp')
process.rootTupleElectronsLoose.InputTagEnDown =  cms.InputTag('electronsEnDown')

### JETS
process.rootTuplePFJets.InputTag = cms.InputTag('smearedAnalysisPatJets')
process.rootTuplePFJets.InputTagSmearedUp   = cms.InputTag('smearedJetsResUp')
process.rootTuplePFJets.InputTagSmearedDown = cms.InputTag('smearedJetsResDown')
process.rootTuplePFJets.InputTagScaledUp    = cms.InputTag('shiftedJetsEnUp')
process.rootTuplePFJets.InputTagScaledDown  = cms.InputTag('shiftedJetsEnDown')

#----------------------------------------------------------------------------------------------------                                                            
# Lepton + Jets filter                                                                                                                                           
#----------------------------------------------------------------------------------------------------                                                            

process.load("Leptoquarks.LeptonJetFilter.leptonjetfilter_cfi")

#### Shared Muon/Electron/Tau Skim                                                                                                                               
process.LJFilter.tauLabel  = cms.InputTag('selectedPatTaus' + postfix)
process.LJFilter.muLabel   = cms.InputTag('selectedPatMuons' + postfixLoose)
process.LJFilter.elecLabel = cms.InputTag('selectedPatElectrons' + postfixLoose)
process.LJFilter.jetLabel  = cms.InputTag('smearedAnalysisPatJets')
process.LJFilter.muonsMin = 1
process.LJFilter.muPT     = 15.0
process.LJFilter.electronsMin = 1
process.LJFilter.elecPT       = 15.0
process.LJFilter.tausMin = 1
process.LJFilter.tauPT   = 15.0
process.LJFilter.jetsMin = 0
process.LJFilter.jetPT   = 20.0
process.LJFilter.counteitherleptontype = True
process.LJFilter.customfilterEMuTauJet2012 = False


process.pdfWeights = cms.EDProducer("PdfWeightProducer",
        # Fix POWHEG if buggy (this PDF set will also appear on output,                                                                                          
        # so only two more PDF sets can be added in PdfSetNames if not "")                                                                                       
        #FixPOWHEG = cms.untracked.string("CT10.LHgrid"),                                                                                                        
        # GenTag = cms.untracked.InputTag("genParticles"),                                                                                                       
        useFirstAsDefault = cms.untracked.bool(True),
        PdfInfoTag = cms.untracked.InputTag("generator"),
        PdfSetNames = cms.untracked.vstring(
           "CT10.LHgrid",
            "MSTW2008nlo68cl.LHgrid",
           "NNPDF20_100.LHgrid"
        )
)

#----------------------------------------------------------------------------------------------------                                                            
# Define the output tree for RootTupleMakerV2                                                                                                                    
#----------------------------------------------------------------------------------------------------                                                            

# RootTupleMakerV2 tree                                                                                                                                          
process.rootTupleTree = cms.EDAnalyzer("RootTupleMakerV2_Tree",
    outputCommands = cms.untracked.vstring(
        'drop *',
        # Event information                                                                                                                                      
        'keep *_rootTupleEvent_*_*',
        'keep *_rootTupleEventSelection_*_*',
        # Single objects                                                                                                                                         
        'keep *_rootTuplePFCandidates_*_*',
        'keep *_rootTuplePFJets_*_*',
        'keep *_rootTupleElectrons_*_*',
        'keep *_rootTupleElectronsLoose_*_*',
        'keep *_rootTupleMuons_*_*',
        'keep *_rootTupleMuonsLoose_*_*',
        'keep *_rootTupleHPSTaus_*_*',
        'keep *_rootTuplePhotons_*_*',
        'keep *_rootTupleVertex_*_*',
        # MET objects for analysis                                                                                                                               
        'keep *_rootTupleCaloMETType1Cor_*_*',
        'keep *_rootTuplePFMET_*_*',
        'keep *_rootTuplePFMETType1Cor_*_*',
        'keep *_rootTuplePFMETType01Cor_*_*',
        'keep *_rootTuplePFMETType01XYCor_*_*',
        'keep *_rootTuplePFMETType01XYCor_*_*',
        # pdf weights                                                                                                                                            
        'keep *_rootTuplePFMETType01XYCorUnclusteredUp_*_*',
        'keep *_rootTuplePFMETType01XYCorUnclusteredDown_*_*',
        'keep *_rootTuplePFMETType01XYCorElectronEnUp_*_*',
        'keep *_rootTuplePFMETType01XYCorElectronEnDown_*_*',
        'keep *_rootTuplePFMETType01XYCorMuonEnUp_*_*',
        
        'keep *_rootTuplePFMETType01XYCorMuonEnDown_*_*',
        'keep *_rootTuplePFMETType01XYCorTauEnUp_*_*',
        'keep *_rootTuplePFMETType01XYCorTauEnDown_*_*',
        'keep *_rootTuplePFMETType01XYCorJetResUp_*_*',
        'keep *_rootTuplePFMETType01XYCorJetResDown_*_*',
        'keep *_rootTuplePFMETType01XYCorJetEnUp_*_*',
        'keep *_rootTuplePFMETType01XYCorJetEnDown_*_*',
        # Trigger objects                                                                                                                                        
        'keep *_rootTupleTrigger_*_*',
        'keep *_rootTupleTriggerObjects_*_*',
        # GEN objects                               
 'keep *_rootTupleGenEventInfo_*_*',
        'keep *_rootTupleGenParticles_*_*',
        'keep *_rootTupleGenJets_*_*',
        'keep *_rootTupleGenElectronsFromWs_*_*',
        'keep *_rootTupleGenElectronsFromZs_*_*',
        'keep *_rootTupleGenMuonsFromWs_*_*',
        'keep *_rootTupleGenMuonsFromZs_*_*',
        'keep *_rootTupleGenTausFromWs_*_*',
        'keep *_rootTupleGenTausFromZs_*_*',
        'keep *_rootTupleGenMETTrue_*_*'
    )
)

#----------------------------------------------------------------------------------------------------                                                            
# Define GEN particle skimmer modules                                                                                                                            
#----------------------------------------------------------------------------------------------------                                                            

process.load ('Leptoquarks.LeptonJetGenTools.genTauMuElFromZs_cfi')
process.load ('Leptoquarks.LeptonJetGenTools.genTauMuElFromWs_cfi')

#----------------------------------------------------------------------------------------------------                                                            
# Define the path                                                                                                                                                
#----------------------------------------------------------------------------------------------------       

#patAddOnSequence = cms.Sequence(
#  getattr( process, 'intermediatePatElectrons' + postfix )
#)
#setattr( process, 'patAddOnSequence' + postfix, patAddOnSequence )

# The paths      


process.p = cms.Path(
  process.goodOfflinePrimaryVertices*
  process.eidMVASequence*
  process.genTausFromWs*
  process.genMuonsFromWs*
  process.genElectronsFromWs*
  process.genTausFromZs*
  process.genMuonsFromZs*
  process.genElectronsFromZs*
  # pdf weights                                                                                                                                                
  process.pdfWeights*
  ## HEEP electron ID                                                                                                                                          
  process.HEEPId*
  # MVA electron ID                                                                                                                                            
  process.eidMVASequence*
  # HEEP rho for isolation correction                                                                                                                          
  process.kt6PFJetsForHEEPIsolation*
  # Good vertices                                                                                                                                              
  process.goodVertices*
  # PFMET corrections                                                                                                                                          
  #process.type0PFMEtCorrection*
  #process.patPFMETtype0Corr*
 # process.producePFMETCorrections*
  #process.pfMEtSysShiftCorrSequencePF*
  # MET filters (required):                                                                                                                                    
  process.EcalDeadCellTriggerPrimitiveFilter*
  process.EcalDeadCellBoundaryEnergyFilter*
  process.HBHENoiseFilterResultProducer*
  process.trackingFailureFilter*
  process.eeBadScFilter*
  process.ecalLaserCorrFilter*
  getattr( process, 'patPF2PATSequence' + postfix )*
  getattr( process, 'patPF2PATSequence' + postfixLoose )*
  #process.looseLeptonSequence * 
  getattr( process, 'intermediatePatElectrons' + postfix )*
  #process.smearedAnalysisPatJets*
  #process.smearedPatPFMetSequencePF*
  process.smearedJetsResDown*
  process.smearedJetsResUp*
  process.shiftedJetsEnDown*
  process.shiftedJetsEnUp*
  process.electronsEnUp*
  process.electronsEnDown*
  process.muonsEnUp*
  process.muonsEnDown*
  process.RawMET*
  process.pfMEtSysShiftCorrSequencePF*
  #process.smearedPatPFMetSequencePF*
  process.patType1CorrectedPFMetType1OnlyPF*
  process.patType1CorrectedPFMetType01OnlyPF*
  process.patType1CorrectedPFMetType01XYOnlyPF*
  process.metUncertaintySequence*
  process.patType1CorrectedPFMetType01XYOnlyPFUnclustup*
  process.patType1CorrectedPFMetType01XYOnlyPFUnclustdown*
  process.LJFilter*
    # Run PAT conversions for electrons                                                                                                                          
  process.patConversions*
    # Re-run full HPS sequence to fully profit from the fix of high pT taus                                                                                      
    process.recoTauClassicHPSSequence*
  (
   # Event information                                                                                                                                          
    process.rootTupleEvent+
    process.rootTupleEventSelection+
    # Single objects                                                                                                                                             
    process.rootTuplePFCandidates+
    process.rootTuplePFJets+
    process.rootTupleElectrons+
    process.rootTupleElectronsLoose+
    process.rootTupleMuons+
    process.rootTupleMuonsLoose+
    #process.rootTupleHPSTaus+
    #process.rootTuplePhotons+
    process.rootTupleVertex+
    # MET objects for analysis                                                                                                                                   
    process.rootTuplePFMET+
    process.rootTuplePFMETType1Cor+
    process.rootTuplePFMETType01Cor+
    process.rootTuplePFMETType01XYCor+
    # MET objects for systematics                                                                                                                                
    process.rootTuplePFMETType01XYCorUnclusteredUp+
    process.rootTuplePFMETType01XYCorUnclusteredDown+
    #process.rootTuplePFMETType01XYCorElectronEnUp+
    #process.rootTuplePFMETType01XYCorElectronEnDown+
    #process.rootTuplePFMETType01XYCorMuonEnDown+
    #process.rootTuplePFMETType01XYCorMuonEnUp+
    #process.rootTuplePFMETType01XYCorTauEnUp+
    #process.rootTuplePFMETType01XYCorTauEnDown+
    #process.rootTuplePFMETType01XYCorJetResUp+
    #process.rootTuplePFMETType01XYCorJetResDown+
    #process.rootTuplePFMETType01XYCorJetEnUp+
    #process.rootTuplePFMETType01XYCorJetEnDown+
    ## Trigger objects                                                                                                                                           
    process.rootTupleTrigger+
    process.rootTupleTriggerObjects+
    # GEN objects                                                                                                                                                
    process.rootTupleGenEventInfo+
    process.rootTupleGenParticles+
    process.rootTupleGenJets+
    process.rootTupleGenElectronsFromWs+
    process.rootTupleGenElectronsFromZs+
    process.rootTupleGenMuonsFromWs+
    process.rootTupleGenMuonsFromZs+
    process.rootTupleGenTausFromWs+
    process.rootTupleGenTausFromZs+
    process.rootTupleGenMETTrue+
    process.rootTupleGenMETCalo
    )*
 process.rootTupleTree
    #process.dump                                                                                                                                                
)
del process.out
del process.outpath

process.schedule = cms.Schedule(process.p)
