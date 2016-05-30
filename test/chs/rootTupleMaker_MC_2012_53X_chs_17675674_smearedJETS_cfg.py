#----------------------------------------------------------------------------------------------------
# Load PAT template
#----------------------------------------------------------------------------------------------------

isData=False
# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## Output Module Configuration (expects a path 'p')                                                                                                                                
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple.root'),
                               # save only events passing the full path                                                                                                            
                                      dropMetaData = cms.untracked.string('ALL'),
                                                      splitLevel = cms.untracked.int32(99),
                               SelectEvents = cms.untracked.PSet(      SelectEvents = cms.vstring('p')     ),
                               # save PAT Layer 1 output; you need a '*' to                                                                                                        
                               # unpack the list of commands 'patEventContent'                                                                                                     
                               outputCommands = cms.untracked.vstring('drop *' )
                               )

process.outpath = cms.EndPath(process.out)


# load the PAT config                                                                                                                                                              
process.load("PhysicsTools.PatAlgos.patSequences_cff")


from PhysicsTools.PatAlgos.tools.coreTools import *
###############################                                                                                                                                                    
####### Parameters ############                                                                                                                                                    
###############################                                                                                                                                                    

inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
if not isData:
    inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
    process.source = cms.Source ("PoolSource",
                                 fileNames=cms.untracked.vstring('file:/afs/cern.ch/work/j/jalmond/FE4C2F81-D0E1-E111-9080-0030487E0A2D.root'      ),
                                 skipEvents=cms.untracked.uint32(993)
                                 )

else :
    inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
    process.source.fileNames = [
        '/store/mc/Summer12/TTJets_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S7_START52_V5-v1/0000/FEBE99BB-3881-E111-B1F3-003048D42DC8.root'
        ]


# Load PF isolation for muons and electrons
#from PhysicsTools.PatAlgos.tools.pfTools import *
#usePFIso ( process )

# Options and Output Report
process.options.wantSummary = True

import os

############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.default.limit = 10
#################################################################

#----------------------------------------------------------------------------------------------------
# Load our RootTupleMakerV2 modules
#----------------------------------------------------------------------------------------------------

process.load('Leptoquarks.RootTupleMakerV2.Ntuple_cff')

# Output ROOT file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string( "file.root" )

)

process.source = cms.Source ("PoolSource",
                             fileNames=cms.untracked.vstring('file:/afs/cern.ch/work/j/jalmond/FE4C2F81-D0E1-E111-9080-0030487E0A2D.root'      ),
                             skipEvents=cms.untracked.uint32(993)
                             )

#----------------------------------------------------------------------------------------------------
# Set global settings (number of events, global tag, input files, etc)
#----------------------------------------------------------------------------------------------------

# GlobalTag
process.GlobalTag.globaltag = 'START53_V27::All'
#process.GlobalTag.globaltag = 'START53_V22::All'

# Events to process1

process.maxEvents.input = 1

process.readAK5PF    = cms.EDAnalyzer('JetCorrectorDBReader',  
        # below is the communication to the database 
        payloadName    = cms.untracked.string('AK5PF'),
                             #        # this is used ONLY for the name of the printed txt files. You can use any name that you like, 
        # but it is recommended to use the GT name that you retrieved the files from.
        globalTag      = cms.untracked.string('START53_V27::All'),
        printScreen    = cms.untracked.bool(False),
        createTextFile = cms.untracked.bool(True)
                                      )

pvSrc = 'offlinePrimaryVertices'

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( maxZ = cms.double(24.0),
                                     minNdof = cms.double(4.0) # this is >= 5                                                                                                      
                                     ),
    src=cms.InputTag(pvSrc)
    )

# manually adding type0 corrections in the sequence: has to be done before setting up PF2PAT                                                                                       

from PhysicsTools.PatAlgos.tools.pfTools import *


process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")

process.producePatPFMETCorrections.replace(process.pfCandMETcorr,
                                           process.type0PFMEtCorrection *
                                           process.patPFMETtype0Corr *
                                           process.pfCandMETcorr
                                           )

postfix="PFlow"


usePF2PAT(process, runPF2PAT=True, jetAlgo='AK5', runOnMC= not isData, postfix=postfix,
          jetCorrections=inputJetCorrLabel,
          pvCollection=cms.InputTag('goodOfflinePrimaryVertices'),
          typeIMetCorrections=True
          )


useGsfElectrons(process,postfix)



# top projections in PF2PAT:
getattr(process,"pfNoPileUp"+postfix).enable = True
getattr(process,"pfNoMuon"+postfix).enable = False #do old style cleaning
getattr(process,"pfNoElectron"+postfix).enable = False #do old style cleaning
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True



  

#----------------------------------------------------------------------------------------------------
# HEEP 4.0 (electron ID) still uses the 2011 definitions of rho for isolation corrections.
# 
# Recipe taken from here:
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaEARhoCorrection#Rho_for_2011_Effective_Areas
#----------------------------------------------------------------------------------------------------


# turn to false when running on data                                                                                                                                               
if isData :
    removeMCMatching( process, ['All'] )


# Do some configuration of the jet substructure things                                                                                                                             
for jetcoll in  (process.patJets,
                 ) :
    if not isData:
        jetcoll.addGenPartonMatch = True
        jetcoll.embedGenJetMatch = False
        jetcoll.getJetMCFlavour = True
    # Add the calo towers and PFCandidates.                                                                                                                                        
    # I'm being a little tricksy here, because I only                                                                                                                              
    # actually keep the products if the "writeFat" switch                                                                                                                          
    # is on. However, this allows for overlap checking                                                                                                                             
    # with the Refs so satisfies most use cases without                                                                                                                            
    # having to add to the object size                                                                                                                                             
    jetcoll.addBTagInfo = False
    jetcoll.embedCaloTowers = True
    jetcoll.embedPFCandidates = True

# Add CATopTag and b-tag info... piggy-backing on b-tag functionality                                                                                                              
process.patJets.addBTagInfo = True
#################################################                                                                                                                                  
#### Fix the PV collections for the future ######                                                                                                                                  
#################################################                                                                                                                                  
for module in [process.patJetCorrFactors,
               process.patJetCorrFactors
               ]:
    module.primaryVertices = "goodOfflinePrimaryVertices"



###############################                                                                                                                                                    
#### Selections Setup #########                                                                                                                                                    
###############################                                                                                                                                                    
process.patJets.addTagInfos = True
process.patJets.tagInfoSources = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfosAOD")
    )



process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJetsForHEEPIsolation = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForHEEPIsolation.Rho_EtaMax = cms.double(2.5)

#----------------------------------------------------------------------------------------------------
# Turn on trigger matching
# 
# Example taken from PAT SWGuide twiki
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTriggerMatchExercise#TrigProduce
#
# Don't include this line (recommended by instructions): removeCleaningFromTriggerMatching( process )
# It will break the matching.
#----------------------------------------------------------------------------------------------------

# load the PAT trigger Python tools
from PhysicsTools.PatAlgos.tools.trigTools import *

# switch on the trigger matching
switchOnTriggerMatching( process, triggerMatchers = [
       # electrons 
        'cleanElectronTriggerMatchHLTSingleElectron8',
        'cleanElectronTriggerMatchHLTSingleElectron17',
        'cleanElectronTriggerMatchHLTSingleElectronWP80',
        'cleanElectronTriggerMatchHLTDoubleElectron',
        # muons
        'cleanMuonTriggerMatchHLTSingleMuon',
        'cleanMuonTriggerMatchHLTSingleMuon5',
        'cleanMuonTriggerMatchHLTSingleMuon8',
        'cleanMuonTriggerMatchHLTSingleMuon12',
        'cleanMuonTriggerMatchHLTSingleMuon17',
        'cleanMuonTriggerMatchHLTSingleMuon24',
        'cleanMuonTriggerMatchHLTDoubleMuon',
        'cleanMuonTriggerMatchHLTSingleIsoMuon'
] )


from PhysicsTools.PatAlgos.tools.jetTools import *



#process.load("Leptoquarks.RootTupleMakerV2.ak5PFchsJets_cff")
#cms.Sequence() += process.ak5PFchsJetsSequence                                                                                                                   \
                                                                                                                                                                   

#addJetCollection(process,cms.InputTag('ak5PFchsJets'),
#    'AK5', 'PFchs',
#    doJetID      = True , # Perform jet ID algorithm and store ID info in the jet                                                                                 \
#                                                                                                                                                                   
#    doJTA        = True , # Perform jet track association and determine jet charge                                                                                \
#                                                                                                                                                                   
#    doBTagging   = True , # Perform b-tagging and store b-tagging info in the jet                                                                                 \
#                                                                                                                                                                   
#    doType1MET   = False, # Don't store Type1 PFMET information. This will be done by the runMEtUncertainties tool.                                               \
#                                                                                                                                                                   
#    jetIdLabel   = "ak5",# Which jet ID label should be used?                                                                                                     \
#                                                                                                                                                                   
#    jetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute']), # Which jet corrections should be used?                                               \
#                                                                                                                                                                   
#    jetCorrLabel = ('AK5PF', ['L1FastJet', 'L2Relative',  'L3Absolute']), # Which jet corrections should be used?                                                \
#                     
#                 genJetCollection = cms.InputTag("ak5GenJets") # Which GEN jets should be used?                                                                   \
#                                                                                                                                                                   
#)

#process.patJetCorrFactorsAK5PFchs.rho = cms.InputTag("kt6PFchsJets","rho")
#process.patJetCorrFactorsAK5PFchs.useRho = cms.bool(True)
#process.patJetCorrFactorsAK5PFchs.payload = cms.string('AK5PFchs')


jetSmearFileName='PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'
jetSmearHistogram='pfJetResolutionMCtoDataCorrLUT'

import RecoMET.METProducers.METSigParams_cfi as jetResolutions

process.smearedAnalysisPatJets = cms.EDProducer("SmearedPATJetProducer",
            src = cms.InputTag("selectedPatJetsPFlow"),
            dRmaxGenJetMatch = cms.string('TMath::Min(0.5, 0.1 + 0.3*TMath::Exp(-0.05*(genJetPt - 10.)))'),
            sigmaMaxGenJetMatch = cms.double(5.),
            inputFileName = cms.FileInPath(jetSmearFileName),
            lutName = cms.string(jetSmearHistogram),
            jetResolutions = jetResolutions.METSignificance_params,
            skipJetSelection = cms.string(
        'jecSetsAvailable & abs(energy - correctedP4("Uncorrected").energy) > (5.*min(energy, correctedP4("Uncorrected").energy))'
        ),
                                         skipRawJetPtThreshold = cms.double(10.), # GeV                                                                           \
                                                                                                                                                                   
                                         skipCorrJetPtThreshold = cms.double(1.e-2)
        )



## IVF and BCandidate producer for Vbb cross check analysis                                                                                                        
process.load('RecoVertex/AdaptiveVertexFinder/inclusiveVertexing_cff')


#----------------------------------------------------------------------------------------------------                                                             \
                                                                                                                                                                   
# Load our RootTupleMakerV2 modules                                                                                                                               \
                                                                                                                                                                   
#----------------------------------------------------------------------------------------------------                                                             \
                                                                                                                                                                   

process.load('Leptoquarks.RootTupleMakerV2.Ntuple_cff')


#----------------------------------------------------------------------------------------------------
# Ad PdFMET and TCMET
# See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#MET_Tools
#----------------------------------------------------------------------------------------------------

from PhysicsTools.PatAlgos.tools.metTools import *
#addPfMET(process, 'PF') # This adds Type1-corrected PFMET
addTcMET(process, 'TC') # This adds raw TC MET

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

process.selectedPatTaus.preselection = cms.string(' tauID("decayModeFinding") > 0.5 ')
process.selectedPatTaus.finalCut     = cms.string(' pt > 15.0 & abs(eta) < 2.5      ')


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
process.patElectrons.userData.userInts.src = cms.VInputTag('HEEPId')

#----------------------------------------------------------------------------------------------------
# Make analysisPatTaus and add them to the cleanPatCandidates sequence
#----------------------------------------------------------------------------------------------------

process.analysisPatTaus = process.selectedPatTaus.clone()
process.analysisPatTaus.preselection = cms.string(
    'tauID("decayModeFinding") > 0.5 &'
    ' tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 &'
    ' tauID("againstMuonLoose3") > 0.5 &'
    ' tauID("againstElectronLooseMVA3") > 0.5'
)
process.analysisPatTaus.finalCut = cms.string('pt > 20. & abs(eta) < 2.3')




process.selectedPatCandidates.replace ( process.selectedPatTaus, process.selectedPatTaus + process.analysisPatTaus )

#----------------------------------------------------------------------------------------------------
# Make analysisPatMuons and add them to the cleanPatCandidates sequence
#----------------------------------------------------------------------------------------------------

process.analysisPatMuons = process.selectedPatMuons.clone()
process.analysisPatMuons.finalCut = cms.string("isGlobalMuon & muonID('GlobalMuonPromptTight') & pt > 10")

process.selectedPatCandidates.replace ( process.selectedPatMuons, process.selectedPatMuons + process.analysisPatMuons )

#----------------------------------------------------------------------------------------------------
# Make analysisPatElectrons and add them to the cleanPatCandidates sequence
#----------------------------------------------------------------------------------------------------

process.analysisPatElectrons = process.selectedPatElectrons.clone()
process.analysisPatElectrons.finalCut = cms.string('pt > 10')

process.selectedPatCandidates.replace ( process.selectedPatElectrons, process.selectedPatElectrons + process.analysisPatElectrons )

#----------------------------------------------------------------------------------------------------
# Add MVA electron ID
#
# MVA electron ID details on this twiki:
# https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentification#MVA_based_Id_in_PAT
#
# Taken from the example:
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/EgammaAnalysis/ElectronTools/test/patTuple_electronId_cfg.py?revision=1.2&view=markup&pathrev=V00-00-09
#----------------------------------------------------------------------------------------------------

process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.mvaID = cms.Sequence(  process.mvaTrigV0 + process.mvaTrigNoIPV0 + process.mvaNonTrigV0 )
process.patElectrons.electronIDSources.mvaTrigV0     = cms.InputTag("mvaTrigV0")  
process.patElectrons.electronIDSources.mvaNonTrigV0  = cms.InputTag("mvaNonTrigV0") 
process.patElectrons.electronIDSources.mvaTrigNoIPV0 = cms.InputTag("mvaTrigNoIPV0")

process.patConversions = cms.EDProducer("PATConversionProducer",
    electronSource = cms.InputTag("selectedPatElectrons")  
)

#----------------------------------------------------------------------------------------------------
# Add the PFJets
# See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#Jet_Tools
# 
# Recommended JEC for PFJets:
# - L1FastJet : https://twiki.cern.ch/twiki/bin/view/CMS/JECAnalysesRecommendations#Corrections
# - L2Relative: https://twiki.cern.ch/twiki/bin/view/CMS/IntroToJEC#Mandatory_Jet_Energy_Corrections
# - L3Absolute: https://twiki.cern.ch/twiki/bin/view/CMS/IntroToJEC#Mandatory_Jet_Energy_Corrections
#----------------------------------------------------------------------------------------------------

from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.goodPatJets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                        filterParams = pfJetIDSelector.clone(),
                                        src = cms.InputTag("selectedPatJets")
                                        )



process.analysisPatJetsAK5PFchs = process.selectedPatJets.clone()
process.analysisPatJetsAK5PFchs.finalCut = cms.string("abs(eta)<2.5 & pt > 10")
process.patDefaultSequence.replace ( process.selectedPatJets, process.selectedPatJets + process.analysisPatJetsAK5PFchs )


switchJetCollection(process,cms.InputTag('ak5CaloJets'),
    doJetID      = True , # Perform jet ID algorithm and store ID info in the jet
    doJTA        = True , # Perform jet track association and determine jet charge
    doBTagging   = True , # Perform b-tagging and store b-tagging info in the jet
    doType1MET   = True , # Store Type1 PFMET information.  Label of resulting PFMET collection is: patMETsAK5Calo
    jetIdLabel   = "ak5",# Which jet ID label should be used?
    jetCorrLabel = ('AK5Calo', ['L1FastJet', 'L2Relative', 'L3Absolute']), # Which jet corrections should be used?
    genJetCollection = cms.InputTag("ak5GenJets") # Which GEN jets should be used?
)

#----------------------------------------------------------------------------------------------------
# Define the systematic shift correction
#----------------------------------------------------------------------------------------------------

#process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
#process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")

#process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc

#----------------------------------------------------------------------------------------------------
# Use the runMetUncertainties tool here
# See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#MET_Systematics_Tools
#----------------------------------------------------------------------------------------------------

print '=' * 60
print "Setting up PFMET from PAT"
print '=' * 60

getattr(process,'patPFMet'+postfix).addGenMET = cms.bool(True)
process.patPFMet.addGenMET = cms.bool(True)


process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")

process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")                                                                                                      
process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc                                                                               

from PhysicsTools.PatUtils.tools.metUncertaintyTools import runMEtUncertainties

runMEtUncertainties(
    process,
    jetCollection           = cms.InputTag('selectedPatJetsPFlow'), 
    doApplySysShiftCorr     = True,  # Apply correction for systematic x/y shift in MET
    doApplyType0corr        = True,  # Apply correction for pileup
    makeType1corrPFMEt      = True,  # Apply correction for jet energy scale
    makeType1p2corrPFMEt    = False, # DO NOT apply correction for unclustered energy (degrades MET resolution)
    makePFMEtByMVA          = False, # We don't use MVA PFMET
    doSmearJets             = True,  # Very important to smear the pfjets (MC ONLY)
    addToPatDefaultSequence = False,  # Add this to the PAT sequence
    electronCollection      = cms.InputTag('analysisPatElectrons'),
    tauCollection           = cms.InputTag('analysisPatTaus'),
    muonCollection          = cms.InputTag('analysisPatMuons'),
    sysShiftCorrParameter   = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc
)
getattr(process,'patPFJetMETtype1p2Corr'+postfix).skipEM = cms.bool(False)
getattr(process,'patPFJetMETtype1p2Corr'+postfix).skipMuons = cms.bool(False)


# Available pat::MET collections for analysis
# - process.patMETsTC                               : raw        TCMET   (NO  jet smearing)
# 
# - process.patMETsRawCalo                          : raw        CaloMET (NO  jet smearing)
# - process.patMETs                                 : Type1      CaloMET (NO  jet smearing)
# 
# - process.patMETsRawPF                            : raw        PFMET   (NO  jet smearing)
# - process.patType1CorrectedPFMet_Type1Only        : Type1      PFMET   (YES jet smearing)
# - process.patType1CorrectedPFMet_Type01Only       : Type0+1    PFMET   (YES jet smearing)
# - process.patType1CorrectedPFMet                  : Type0+1+XY PFMET   (YES jet smearing) <-- Recommended for analysis
# 
# Available pat::MET collections for systematic studies
# - process.patType1CorrectedPFMetElectronEnUp      : Type0+1+XY PFMET   (YES jet smearing), Electron energy shifted up 
# - process.patType1CorrectedPFMetElectronEnDown    : Type0+1+XY PFMET   (YES jet smearing), Electron energy shifted down
# - process.patType1CorrectedPFMetMuonEnUp          : Type0+1+XY PFMET   (YES jet smearing), Muon energy shifted up
# - process.patType1CorrectedPFMetMuonEnDown        : Type0+1+XY PFMET   (YES jet smearing), Muon energy shifted down 
# - process.patType1CorrectedPFMetTauEnUp           : Type0+1+XY PFMET   (YES jet smearing), Tau energy shifted up   
# - process.patType1CorrectedPFMetTauEnDown         : Type0+1+XY PFMET   (YES jet smearing), Tau energy shifted down 
# - process.patType1CorrectedPFMetJetResUp          : Type0+1+XY PFMET   (YES jet smearing), Jet resolution smeared up
# - process.patType1CorrectedPFMetJetResDown        : Type0+1+XY PFMET   (YES jet smearing), Jet resolution smeared down
# - process.patType1CorrectedPFMetJetEnUp           : Type0+1+XY PFMET   (YES jet smearing), Jet energy shifted up   
# - process.patType1CorrectedPFMetJetEnDown         : Type0+1+XY PFMET   (YES jet smearing), Jet energy shifted down  
# - process.patType1CorrectedPFMetUnclusteredEnUp   : Type0+1+XY PFMET   (YES jet smearing), Unclustered energy shifted up  
# - process.patType1CorrectedPFMetUnclusteredEnDown : Type0+1+XY PFMET   (YES jet smearing), Unclustered energy shifted down
# 
# Available shifted object collections:
# - process.shiftedPatElectronsEnUp                 : pat electrons, energy scale shifted up
# - process.shiftedPatElectronsEnDown               : pat electrons, energy scale shifted down
# - process.shiftedPatMuonsEnUp                     : pat muons    , energy scale shifted up
# - process.shiftedPatMuonsEnDown                   : pat muons    , energy scale shifted down
# - process.shiftedPatTausEnUp                      : pat taus     , energy scale shifted up
# - process.shiftedPatTausEnDown                    : pat taus     , energy scale shifted down
# - process.smearedPatJetsAK5PFchs                     : pat jets     , energy resolution smeared to match data  <-- Recommended for analysis
# - process.smearedPatJetsAK5PFchsresUp                : pat jets     , energy resolution smeared worse data
# - process.smearedPatJetsAK5PFchsresDown              : pat jets     , energy resolution sharpened better than data
# - process.shiftedPatJetsAK5PFchsenUpForCorrMEt       : pat jets     , energy scale shifted up  
# - process.shiftedPatJetsAK5PFchsenDownForCorrMEt     : pat jets     , energy scale shifted down 
# 
# Notes:
# - No Type0 corrections are available for CaloMET
# - No Type0 or Type1 corrections are available for TCMET
# - No Type2 corrections recommended for any MET, since they degrade the MET resolution
#
# Thanks to Jared Sturdy, Christian Veelken, and Guler Karapinar
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Add the raw CaloMET and raw PFMET
# Raw TCMET is already included: patMETsRawTC
#----------------------------------------------------------------------------------------------------

# CaloMET: raw

process.patMETsRawCalo = process.patMETsTC.clone()
process.patMETsRawCalo.metSource = cms.InputTag("met")

# PFMET: raw

process.patMETsRawPF = process.patMETsTC.clone()
process.patMETsRawPF.metSource = cms.InputTag("pfMet")

# PFMET: Type1, with jet smearing

process.patType1CorrectedPFMetType1Only = process.patType1CorrectedPFMet.clone()
process.patType1CorrectedPFMetType1Only.srcType1Corrections = cms.VInputTag(
    cms.InputTag("patPFJetMETtype1p2Corr","type1"), 
    # cms.InputTag("patPFMETtype0Corr"), 
    # cms.InputTag("pfMEtSysShiftCorr")
)

# PFMET: Type0+1, with jet smearing

process.patType1CorrectedPFMetType01Only = process.patType1CorrectedPFMet.clone()
process.patType1CorrectedPFMetType01Only.srcType1Corrections = cms.VInputTag(
    cms.InputTag("patPFJetMETtype1p2Corr","type1"), 
    cms.InputTag("patPFMETtype0Corr"), 
    # cms.InputTag("pfMEtSysShiftCorr")
)


#----------------------------------------------------------------------------------------------------
# This is MC, so analyze the smeared PFJets by default
#----------------------------------------------------------------------------------------------------

#process.rootTuplePFJets.InputTag = cms.InputTag('smearedAnalysisPatJetsAK5PFchs')
#process.rootTuplePFJets.InputTag = cms.InputTag('smearedAnalysisPatJetsAK5PFchs')       
#process.rootTuplePFJets.InputTag = cms.InputTag('analysisPatJetsAK5PFchs')
process.rootTuplePFJets.InputTag = cms.InputTag('smearedAnalysisPatJets')
#process.rootTuplePFJets.InputTagSmearedUp   = cms.InputTag('smearedAnalysisPatJetsAK5PFchsResUp')                                 
#process.rootTuplePFJets.InputTagSmearedDown = cms.InputTag('smearedAnalysisPatJetsAK5PFchsResDown')                                 
#process.rootTuplePFJets.InputTagScaledUp    = cms.InputTag('shiftedAnalysisPatJetsAK5PFchsEnUpForCorrMEt')                                 
#process.rootTuplePFJets.InputTagScaledDown  = cms.InputTag('shiftedAnalysisPatJetsAK5PFchsEnDownForCorrMEt')     
process.rootTuplePFJets.InputTagSmearedUp = cms.InputTag('smearedAnalysisPatJets')
process.rootTuplePFJets.InputTagSmearedDown = cms.InputTag('smearedAnalysisPatJets')
process.rootTuplePFJets.InputTagScaledUp    = cms.InputTag('smearedAnalysisPatJets')
process.rootTuplePFJets.InputTagScaledDown  =  cms.InputTag('smearedAnalysisPatJets')
#----------------------------------------------------------------------------------------------------
# Set Lepton-Gen Matching Parameters
#----------------------------------------------------------------------------------------------------

process.load("Leptoquarks.RootTupleMakerV2.leptonGenMatching_cfi")
process.patDefaultSequence.replace( process.electronMatch, process.elMatch )
process.patElectrons.genParticleMatch = cms.VInputTag( cms.InputTag("elMatch") )
process.patDefaultSequence.replace( process.muonMatch, process.muMatch )
process.patMuons.genParticleMatch = cms.VInputTag( cms.InputTag("muMatch") )
process.patDefaultSequence.replace( process.tauMatch, process.tauLepMatch )
process.patTaus.genParticleMatch = cms.VInputTag( cms.InputTag("tauLepMatch") )
process.patDefaultSequence.replace( process.tauGenJetMatch, process.tauJetMatch )
process.patTaus.genJetMatch = cms.InputTag("tauJetMatch")

#----------------------------------------------------------------------------------------------------
# Lepton + Jets filter
#----------------------------------------------------------------------------------------------------

process.load("Leptoquarks.LeptonJetFilter.leptonjetfilter_cfi")

#### Shared Muon/Electron/Tau Skim
process.LJFilter.tauLabel  = cms.InputTag("selectedPatTaus")                        
process.LJFilter.muLabel   = cms.InputTag("selectedPatMuons")
process.LJFilter.elecLabel = cms.InputTag("selectedPatElectrons")
process.LJFilter.jetLabel  = cms.InputTag("smearedAnalysisPatJetsAK5PFchs")
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
# -- WARNING :
# "customfilterEMuTauJet2012" configuration is hard-coded.
# (see: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Leptoquarks/LeptonJetFilter/src/LeptonJetFilter.cc )
# "customfilterEMuTauJet2012" is the desired mode of operation for the Lepton+Jets Filter in 2012.

#----------------------------------------------------------------------------------------------------
# PDF weights
#----------------------------------------------------------------------------------------------------

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
        'keep *_rootTupleCaloJets_*_*',
        'keep *_rootTupleElectrons_*_*',
        'keep *_rootTupleMuons_*_*',
        'keep *_rootTupleHPSTaus_*_*',
        'keep *_rootTuplePhotons_*_*',
        'keep *_rootTupleVertex_*_*',
        # MET objects for analysis
        'keep *_rootTupleTCMET_*_*',
        'keep *_rootTupleCaloMET_*_*',
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
        'keep *_rootTupleGenMETTrue_*_*',
        'keep *_rootTupleGenMETCalo_*_*'       
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



#process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.p = cms.Path(
    # gen particle skimmer modules#
    #process.readAK5PF*
    process.goodOfflinePrimaryVertices*
    process.primaryVertexFilter*
    process.inclusiveVertexing*
    process.ak5PFchsJetsSequence*

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
    process.mvaID*
    # HEEP rho for isolation correction
    process.kt6PFJetsForHEEPIsolation*
    # Good vertices
    process.goodVertices*
    # PFMET corrections
    process.type0PFMEtCorrection*
    process.patPFMETtype0Corr*
    process.producePFMETCorrections*
    process.pfMEtSysShiftCorrSequence*
    process.metUncertaintySequencePFlow*
    # MET filters (required):
    process.EcalDeadCellTriggerPrimitiveFilter*
    process.EcalDeadCellBoundaryEnergyFilter*
    process.HBHENoiseFilterResultProducer*
    process.trackingFailureFilter*
    process.eeBadScFilter*
    process.ecalLaserCorrFilter*

    getattr(process,"patPF2PATSequence"+postfix)*
    #process.patDefaultSequence*
    process.goodPatJets*
    #process.prunedGenParticles*

    process.smearedAnalysisPatJets*
    process.smearedAnalysisPatJets2*

    # Add the pileup MVA to the jets
    #process.puJetIdSqeuenceChs*
    # MET producers
    process.patMETsRawCalo*
    process.patMETsRawPF*
    process.patType1CorrectedPFMetType1Only*
    process.patType1CorrectedPFMetType01Only*
    # L+J Filter
    process.LJFilter*  
    # Run PAT conversions for electrons
    process.patConversions*
    # Re-run full HPS sequence to fully profit from the fix of high pT taus
    process.recoTauClassicHPSSequence*
    # RootTupleMakerV2
    (
    # Event information
    process.rootTupleEvent+
    process.rootTupleEventSelection+
    # Single objects
    process.rootTuplePFCandidates+
    process.rootTuplePFJets+
    process.rootTupleCaloJets+
    process.rootTupleElectrons+
    process.rootTupleMuons+
    process.rootTupleHPSTaus+
    process.rootTuplePhotons+
    process.rootTupleVertex+
    # MET objects for analysis
    process.rootTupleTCMET+
    process.rootTupleCaloMET+
    process.rootTupleCaloMETType1Cor+
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
    process.rootTuplePFMETType01XYCorMuonEnDown+
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
    # Put everything into the tree
    process.rootTupleTree
    #process.dump
)


#----------------------------------------------------------------------------------------------------
# Dump if necessary
#----------------------------------------------------------------------------------------------------

#process.dump = cms.OutputModule("PoolOutputModule",
#                                outputCommands = cms.untracked.vstring(
#                                'keep *',
#                                ),
#                                fileName = cms.untracked.string('dump.root')
#                                )
#process.DUMP    = cms.EndPath (process.dump)

# Delete predefined Endpath (needed for running with CRAB)
del process.out
del process.outpath


#----------------------------------------------------------------------------------------------------
# Run the path
#----------------------------------------------------------------------------------------------------

process.schedule = cms.Schedule(process.p)
