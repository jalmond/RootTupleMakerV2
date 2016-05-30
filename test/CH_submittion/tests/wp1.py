##### specify if is Data
isData=False

runOnMC= not isData
##############################################################################################################
#### configuration
##############################################################################################################

from Leptoquarks.RootTupleMakerV2.patRefSel_refMuJets import *

# list of input files
inputFiles = ['file:/afs/cern.ch/work/j/jalmond/FE4C2F81-D0E1-E111-9080-0030487E0A2D.root'] # overwritten, if "useRelVals" is 'True'

# output file
outputFile = 'CH_submit_mc_53X.root'
GLOBALTAG='START53_V27::All'
if isData:
  outputFile = 'CH_submit_data_53X.root'
  GLOBALTAG='FT53_V21A_AN6::All'



# maximum number of events
maxEvents = 1 # reduce for testiang


### Particle flow                                                                                                                                                      
postfix = "PF"

typeIMetCorrections = True

# event frequency of Fwk report
fwkReportEvery = 1000

# switch for 'TrigReport'/'TimeReport' at job end
wantSummary = True


###                              End of constants                            ###
###                                                                          ###
### ======================================================================== ###


###
### Basic configuration
###

from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.GlobalTag.globaltag = GLOBALTAG
# Input files                                                                                                                                                             
process.maxEvents.input = maxEvents
process.source.fileNames = inputFiles

process.load( "Leptoquarks.RootTupleMakerV2.patRefSel_goodVertex_cfi" )

### Check JECs
jecLevels = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
if not runOnMC:
  jecLevels = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])

### Switch configuration

### Add in type - corrections
##____________________________________________________________________________||
process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")

process.producePatPFMETCorrections.replace(
  process.pfCandMETcorr,
  process.type0PFMEtCorrection *
  process.patPFMETtype0Corr *
  process.pfCandMETcorr 
  )

print '=' * 60
print "Setting up usePF2PAT"
print '=' * 60

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
print '=' * 60
print "Finished setting up usePF2PAT"
print '=' * 60


## pile up corrections
from CommonTools.ParticleFlow.Tools.enablePileUpCorrection import *
#enablePileUpCorrectionInPF2PAT( process, postfix, sequence = "patPF2PATSequence"+postfix)
enablePileUpCorrection( process, postfix, sequence = "patPF2PATSequence"+postfix)

print '=' * 60
print "Ran enablePileUpCorrection"
print '=' * 60
################################################
################################################   ADD TOP PROJECTION CUTS TO PF Objects
################################################

print '=' * 60
print "Setting up top projection"
print '=' * 60

print '=' * 60
print "----- pfMuons"
print '=' * 60
process.pfSelectedMuonsPF.cut = 'abs(eta)<2.5 && pt>10.'
process.pfIsolatedMuonsPF.doDeltaBetaCorrection = True
process.pfIsolatedMuonsPF.deltaBetaFactor = -0.5
process.pfIsolatedMuonsPF.isolationCut = 0.20

print '=' * 60
print "----- pfElectrons"
print '=' * 60
process.pfIsolatedElectronsPF.doDeltaBetaCorrection = True
# Turn on the delta-beta corrections

process.load( "EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi" )
process.load('EgammaAnalysis.ElectronTools.electronIsolatorFromEffectiveArea_cfi')
process.pfIdentifiedElectronsPF = cms.EDFilter("ElectronIDPFCandidateSelector",
                                              recoGsfElectrons = cms.InputTag("gsfElectrons"),
                                              electronIdMap = cms.InputTag("mvaTrigV0"),
                                              electronIdCut = cms.double(0.0),
                                              src = cms.InputTag("pfElectronsFromVertex"+postfix))

process.pfSelectedElectronsPF.src = 'pfIdentifiedElectronsPF'
process.pfSelectedElectronsPF.cut = 'abs(eta)<2.5 && pt>20. && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits<2'

process.pfIsolatedElectronsPF.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFId"+postfix))
process.pfIsolatedElectronsPF.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFId"+postfix), 
                                                                      cms.InputTag("elPFIsoValueGamma03PFId"+postfix))

process.elPFIsoValueEA03.pfElectrons = cms.InputTag('pfSelectedElectrons' + postfix)
process.pfIsolatedElectronsPF.deltaBetaIsolationValueMap = 'elPFIsoValueEA03' # new Effective Area edm::ValueMap
process.pfIsolatedElectronsPF.doDeltaBetaCorrection = True                    # not really a 'deltaBeta' correction, but it serves
process.pfIsolatedElectronsPF.deltaBetaFactor = -1.0
process.pfIsolatedElectronsPF.isolationCut = 0.15

process.patPF2PATSequencePF.replace( process.pfSelectedElectronsPF,
                                   process.pfIdentifiedElectronsPF + process.pfSelectedElectronsPF + process.elPFIsoValueEA03 )


print '=' * 60
print "----- PATElectron isolation"
print '=' * 60
getattr( process, 'patElectrons' + postfix ).isolationValues.pfNeutralHadrons   = cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
getattr( process, 'patElectrons' + postfix ).isolationValues.pfChargedAll       = cms.InputTag( 'elPFIsoValueChargedAll03PFId' + postfix )
getattr( process, 'patElectrons' + postfix ).isolationValues.pfPUChargedHadrons = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
getattr( process, 'patElectrons' + postfix ).isolationValues.pfPhotons          = cms.InputTag( 'elPFIsoValueGamma03PFId' + postfix )
getattr( process, 'patElectrons' + postfix ).isolationValues.pfChargedHadrons   = cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix )
getattr( process, 'patElectrons' + postfix ).isolationValues.user =   cms.VInputTag("elPFIsoValueCharged03PFIdPF","elPFIsoValueNeutral03PFIdPF","elPFIsoValueGamma03PFIdPF","elPFIsoValuePU03PFIdPF","elPFIsoValueChargedAll03PFIdPF")        



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
# enable delta beta correction for muon selection in PF2PAT?
#getattr(process,"pfIsolatedMuons"+postfix).doDeltaBetaCorrection = False

print '=' * 60
print "----- Add patPFMETtype0Corr to patType1CorrectedPFMet"
print '=' * 60



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

from Leptoquarks.RootTupleMakerV2.patRefSel_ch_cfi import *

# remove MC matching, object cleaning, objects etc.
if isData:
  runOnData( process
           , names = [ 'PFAll' ]
           , postfix = postfix
           )

from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')


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


### ELectrons
intermediatePatElectrons.src = cms.InputTag( 'selectedPatElectrons' + postfix )
setattr( process, 'intermediatePatElectrons' + postfix, intermediatePatElectrons )
### Muons
process.out.outputCommands.append( 'keep *_pfMuonsFromVertex*_*_*' )

###
### Selection configuration
###

### Muons

# Use the good primary vertices everywhere.                                                                                                                               
for imod in [process.patMuonsPF,
             process.patElectronsPF]:
    imod.pvSrc = "goodOfflinePrimaryVertices"
    imod.embedTrack = True
    imod.usePV = False

getattr( process, 'selectedPatMuons' + postfix ).cut = LQmuonCut
getattr( process, 'patMuons' + postfix ).user = cms.VInputTag("muPFIsoValueCharged03PF","muPFIsoValueNeutral03PF","muPFIsoValueGamma03PF","muPFIsoValuePU03PF","muPFIsoValueChargedAll03PF")

### Jets
getattr( process, 'selectedPatJets'  + postfix ).cut = LQJetCut

### Electrons
getattr( process, 'selectedPatElectrons' + postfix ).cut = LQelectronCut

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
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.default.limit = 10
################################################################
process.source = cms.Source ("PoolSource",
                             fileNames=cms.untracked.vstring('file:/afs/cern.ch/work/j/jalmond/FE4C2F81-D0E1-E111-9080-0030487E0A2D.root'      ),
)

process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJetsForHEEPIsolation = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForHEEPIsolation.Rho_EtaMax = cms.double(2.5)

#----------------------------------------------------------------------------------------------------                                                            
# Add PFMET and TCMET                                                                                                                                            
# See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#MET_Tools                                                                                   
#----------------------------------------------------------------------------------------------------                                                            

#from PhysicsTools.PatAlgos.tools.metTools import *
#addTcMET(process, 'TC') # This adds raw TC MET                                                                                                                   


#----------------------------------------------------------------------------------------------------                                                            
# Add MET filters                                                                                                                                                
#----------------------------------------------------------------------------------------------------                                                            

process.load("Leptoquarks.RootTupleMakerV2.metFilters_cfi")


#----------------------------------------------------------------------------------------------------                                                            
# Rerun full HPS sequence to fully profit from the fix of high pT taus                                                                                           
#----------------------------------------------------------------------------------------------------                                                            

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
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

process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc
#process.pfMEtSysShiftCorrPF.src =  cms.InputTag('pfMet')


#---------------------------------------------------------------------------------------------------- 
# Use the runMetUncertainties tool here                                                                                                                          
# See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#MET_Systematics_Tools                                                                       
#----------------------------------------------------------------------------------------------------                                                            

from PhysicsTools.PatUtils.tools.metUncertaintyTools import runMEtUncertainties

runMEtUncertainties(
    process,
    jetCollection           = cms.InputTag('selectedPatJets' + postfix),
    doApplySysShiftCorr     = False,  # Apply correction for systematic x/y shift in MET                                                                          
    doApplyType0corr        = False,  # Apply correction for pileup                                                                                               
    makeType1corrPFMEt      = False,  # Apply correction for jet energy scale                                                                                     

    makeType1p2corrPFMEt    = False, # DO NOT apply correction for unclustered energy (degrades MET resolution)                                                  
    makePFMEtByMVA          = False, # We don't use MVA PFMET                                                                                                    
    doSmearJets             = False,  # Very important to smear the pfjets (MC ONLY)                                                                              
    addToPatDefaultSequence = False,  # Add this to the PAT sequence                                                                                             
    electronCollection      = cms.InputTag('selectedPatElectrons' + postfix),
    tauCollection           = cms.InputTag('selectedPatTaus' + postfix),
    muonCollection          = cms.InputTag('selectedPatMuons' + postfix),
    sysShiftCorrParameter   = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc,
    postfix='PF'
)


###### MET                                                                                                                                                               
print '=' * 60
print "----- Add patPFMETtype0Corr to patType1CorrectedPFMet"
print '=' * 60

getattr(process,'patMETs'+postfix).metSource = cms.InputTag("patType1CorrectedPFMet"+postfix)
getattr(process,'patType1CorrectedPFMet'+postfix).srcType1Corrections = cms.VInputTag(
    cms.InputTag("patPFJetMETtype1p2Corr"+postfix,"type1")
    #cms.InputTag("patPFMETtype0Corr"+postfix)                                                                                                                          
)

### FIX run metuncertaunty                                                                                                                                               
#getattr(process,'patType1CorrectedPFMet').srcType1Corrections = cms.VInputTag(                                                                                        
#    cms.InputTag("patPFJetMETtype1p2Corr","type1")                                                                                                                     
#    #cms.InputTag("patPFMETtype0Corr"+postfix)                                                                                                                       
#                                                                                                                                                                       
#)                                                                                                                                                                       

process.patType1CorrectedPFMetType01OnlyPF = process.patType1CorrectedPFMetPF.clone()
process.patType1CorrectedPFMetType01OnlyPF.srcType1Corrections = cms.VInputTag(
    cms.InputTag("patPFJetMETtype1p2CorrPF","type1"),
    cms.InputTag("patPFMETtype0CorrPF"),
    # cms.InputTag("pfMEtSysShiftCorr")                                                                                                                                                                          
)

from PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi import patMETs

process.apatPFMet = patMETs.clone(
    metSource = cms.InputTag('pfMet'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue')
)

process.patPFMet1 = patMETs.clone(
    metSource = cms.InputTag('patType1CorrectedPFMetPF'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue')
)


process.patPFMet1b = patMETs.clone(
    metSource = cms.InputTag('patType1CorrectedPFMetPF'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue')
)



#these flags are false for '+postfix' mets by default, but true for non-postfix ones!                                                                           \


jetSmearFileName='PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'
jetSmearHistogram='pfJetResolutionMCtoDataCorrLUT'

import RecoMET.METProducers.METSigParams_cfi as jetResolutions

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

process.rootTuplePFJets.InputTag = cms.InputTag('smearedAnalysisPatJets')
process.rootTupleElectrons.InputTag = cms.InputTag('selectedPatElectrons' + postfix)
process.rootTupleMuons.InputTag = cms.InputTag('selectedPatMuons'+ postfix)
process.rootTupleMuons.InputTagEnUp  = cms.InputTag('selectedPatMuons' + postfix)
process.rootTupleMuons.InputTagEnDown  = cms.InputTag('selectedPatMuons' + postfix)
process.rootTupleElectrons.InputTagEnUp  =  cms.InputTag('selectedPatElectrons' + postfix)
process.rootTupleElectrons.InputTagEnDown =  cms.InputTag('selectedPatElectrons' + postfix)

process.rootTuplePFJets.InputTagSmearedUp   = cms.InputTag('smearedJetsResUp')
process.rootTuplePFJets.InputTagSmearedDown = cms.InputTag('smearedJetsResDown')
process.rootTuplePFJets.InputTagScaledUp    = cms.InputTag('shiftedAnalysisPatJetsAK5PFenUpForCorrMEt')
process.rootTuplePFJets.InputTagScaledDown  = cms.InputTag('shiftedAnalysisPatJetsAK5PFenDownForCorrMEt')

#----------------------------------------------------------------------------------------------------                                                            
# Lepton + Jets filter                                                                                                                                           
#----------------------------------------------------------------------------------------------------                                                            

process.load("Leptoquarks.LeptonJetFilter.leptonjetfilter_cfi")

#### Shared Muon/Electron/Tau Skim                                                                                                                               
process.LJFilter.tauLabel  = cms.InputTag('selectedPatTaus' + postfix)
process.LJFilter.muLabel   = cms.InputTag('selectedPatMuons' + postfix)
process.LJFilter.elecLabel = cms.InputTag('selectedPatElectrons' + postfix)
process.LJFilter.jetLabel  = cms.InputTag('smearedAnalysisPatJets')
process.LJFilter.muonsMin = 0
process.LJFilter.muPT     = 1.0
process.LJFilter.electronsMin = 0
process.LJFilter.elecPT       = 1.0
process.LJFilter.tausMin = 0
process.LJFilter.tauPT   = 1.0
process.LJFilter.jetsMin = 0
process.LJFilter.jetPT   = 1.0
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
        'keep *_rootTupleMuons_*_*',
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
  process.pfMEtSysShiftCorrSequence*
  # MET filters (required):                                                                                                                                    
  process.EcalDeadCellTriggerPrimitiveFilter*
  process.EcalDeadCellBoundaryEnergyFilter*
  process.HBHENoiseFilterResultProducer*
  process.trackingFailureFilter*
  process.eeBadScFilter*
  process.ecalLaserCorrFilter*
  getattr( process, 'patPF2PATSequence' + postfix )*
  #process.patDefaultSequence*
  #getattr( process, 'patAddOnSequence' + postfix )*
  getattr( process, 'intermediatePatElectrons' + postfix )*
  process.metUncertaintySequencePF*
  process.smearedAnalysisPatJets*
  process.smearedJetsResDown*
  process.smearedJetsResUp*
  process.shiftedJetsEnDown*
  process.shiftedJetsEnUp*
  process.apatPFMet*
  process.patPFMet1*
  process.patPFMet1b*
#  process.patMETsRawCalo*
#  process.patMETsRawPF*
  #process.patType1CorrectedPFMetType1OnlyPF*
  process.patType1CorrectedPFMetType01OnlyPF*
  #process.LJFilter*
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
    process.rootTupleMuons+
    process.rootTupleHPSTaus+
    process.rootTuplePhotons+
    process.rootTupleVertex+
    # MET objects for analysis                                                                                                                                   
    #process.rootTupleTCMET+
    # process.rootTupleCaloMET+
    #process.rootTupleCaloMETType1Cor+
    process.rootTuplePFMET+
    process.rootTuplePFMETType1Cor+
    process.rootTuplePFMETType01Cor+
    process.rootTuplePFMETType01XYCor+
    # MET objects for systematics                                                                                                                                
    process.rootTuplePFMETType01XYCorUnclusteredUp+
    process.rootTuplePFMETType01XYCorUnclusteredDown+
    process.rootTuplePFMETType01XYCorElectronEnUp+
    process.rootTuplePFMETType01XYCorElectronEnDown+
    process.rootTuplePFMETType01XYCorMuonEnUp+
    process.rootTuplePFMETType01XYCorTauEnUp+
    process.rootTuplePFMETType01XYCorTauEnDown+
    process.rootTuplePFMETType01XYCorJetResUp+
    process.rootTuplePFMETType01XYCorJetResDown+
    process.rootTuplePFMETType01XYCorJetEnUp+
    process.rootTuplePFMETType01XYCorJetEnDown+
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
