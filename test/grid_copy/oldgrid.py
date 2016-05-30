#----------------------------------------------------------------------------------------------------
# Load PAT template
#----------------------------------------------------------------------------------------------------

isData=False
# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

# Load PF isolation for muons and electrons
#from PhysicsTools.PatAlgos.tools.pfTools import *
#usePFIso ( process )

from PhysicsTools.PatAlgos.tools.coreTools import *

postfix = "PFlow"

inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])

if isData:
    inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
    
process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")
#process.producePatPFMETCorrections.replace(process.pfCandMETcorr,
#                                           process.type0PFMEtCorrection *
#                                           process.patPFMETtype0Corr *
#                                           process.pfCandMETcorr
#
#                                           )
from PhysicsTools.PatAlgos.tools.pfTools import *

#useGsfElectrons(process,postfix)

#usePFIso(process,postfix)    


usePF2PAT(process, runPF2PAT=True, jetAlgo='AK5', runOnMC=not isData, postfix=postfix, typeIMetCorrections=True,
          jetCorrections=inputJetCorrLabel,
          pvCollection=cms.InputTag('goodOfflinePrimaryVertices')
          )


process.pfIsolatedElectronsPFlow.doDeltaBetaCorrection = True
process.pfIsolatedMuonsPFlow.doDeltaBetaCorrection = False

process.patMuonsPFlow.isolationValues.user = cms.VInputTag("muPFIsoValueCharged03PFlow","muPFIsoValueNeutral03PFlow","muPFIsoValueGamma03PFlow","muPFIsoValuePU03PFlow","muPFIsoValueChargedAll03PFlow")

process.patElectronsPFlow.isolationValues.user = cms.VInputTag("elPFIsoValueCharged03PFIdPFlow","elPFIsoValueNeutral03PFIdPFlow","elPFIsoValueGamma03PFIdPFlow","elPFIsoValuePU03PFIdPFlow","elPFIsoValueChargedAll03PFIdPFlow")



#False == taus also in jet collection  (This is actually not needed since TP are applied later)          
process.pfNoTauPFlow.enable = False

## pile up corrections
from CommonTools.ParticleFlow.Tools.enablePileUpCorrection import enablePileUpCorrectionInPF2PAT
enablePileUpCorrectionInPF2PAT( process, postfix, sequence = "patPF2PATSequence"+postfix)


from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
for jetcoll in  (process.patJetsPFlow,
                 ) :
    
    if not isData:
        jetcoll.addGenPartonMatch = True
        jetcoll.embedGenJetMatch = False
        jetcoll.getJetMCFlavour = True
    jetcoll.addBTagInfo = False
    jetcoll.embedCaloTowers = True
    jetcoll.embedPFCandidates = True



for imod in [process.patMuonsPFlow,
             process.patElectronsPFlow] :
             
        imod.pvSrc = "goodOfflinePrimaryVertices"
        imod.embedTrack = True


 ## electron ID tool
process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
process.patElectronsPFlow.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
process.patElectronsPFlow.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
process.patDefaultSequencePFlow.replace( process.patElectronsPFlow, process.eidMVASequence * process.patElectronsPFlow )

process.patConversions = cms.EDProducer("PATConversionProducer",
    electronSource = cms.InputTag("selectedPatElectronsPFlow")
)


## Add this to main sequence
#process.p += getattr(process,"patPF2PATSequence"+postfix)

# top projections in PF2PAT:
getattr(process,"pfNoPileUp"+postfix).enable = True
getattr(process,"pfNoMuon"+postfix).enable = True
getattr(process,"pfNoElectron"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True
# verbose flags for the PF2PAT modules
getattr(process,"pfNoMuon"+postfix).verbose = False
# enable delta beta correction for muon selection in PF2PAT?
getattr(process,"pfIsolatedMuons"+postfix).doDeltaBetaCorrection = False


# change PV in Jet corrector
for module in [
               process.patJetCorrFactorsPFlow
               ]:
    module.primaryVertices = "goodOfflinePrimaryVertices"



#### object selection
process.selectedPatJetsPFlow.cut = cms.string("pt > 10")
process.patJetsPFlow.addTagInfos = True
process.patJetsPFlow.tagInfoSources = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfosAODPFlow")
    )


process.selectedPatElectronsPFlow.cut = cms.string('pt > 10.0 & abs(eta) < 2.5')
process.patElectronsPFlow.embedTrack = cms.bool(True)

process.selectedPatMuonsPFlow.cut = cms.string("pt > 10.0 & abs(eta) < 2.5")
process.patMuonsPFlow.embedTrack = cms.bool(True)
process.selectedPatTausPFlow.cut = cms.string("pt > 10.0 & abs(eta) < 3")
process.selectedPatTaus.cut = cms.string("pt > 10.0 & abs(eta) < 3")



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
                             skipEvents=cms.untracked.uint32(1108),
                             )

#----------------------------------------------------------------------------------------------------
# Set global settings (number of events, global tag, input files, etc)
#----------------------------------------------------------------------------------------------------

# GlobalTag
process.GlobalTag.globaltag = 'START53_V27::All'

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



#----------------------------------------------------------------------------------------------------
# HEEP 4.0 (electron ID) still uses the 2011 definitions of rho for isolation corrections.
# 
# Recipe taken from here:
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaEARhoCorrection#Rho_for_2011_Effective_Areas
#----------------------------------------------------------------------------------------------------

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
switchOnTriggerMatching( process, sequence = "patPF2PATSequence"+postfix , triggerMatchers = [
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

#----------------------------------------------------------------------------------------------------
# Add PFMET and TCMET
# See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#MET_Tools
#----------------------------------------------------------------------------------------------------

from PhysicsTools.PatAlgos.tools.metTools import *
#addPfMET(process, 'PF') # This adds Type1-corrected PFMET (This is added now using PFtoPAT)
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

process.cleanPatTaus.preselection = cms.string(' tauID("decayModeFinding") > 0.5 ')
process.cleanPatTaus.finalCut     = cms.string(' pt > 15.0 & abs(eta) < 2.5      ')


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
process.patElectronsPFlow.userData.userInts.src = cms.VInputTag('HEEPId')

#----------------------------------------------------------------------------------------------------
# Make analysisPatTaus and add them to the cleanPatCandidates sequence
#----------------------------------------------------------------------------------------------------

process.analysisPatTaus = process.cleanPatTaus.clone()
process.analysisPatTaus.preselection = cms.string(
    'tauID("decayModeFinding") > 0.5 &'
    ' tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 &'
    ' tauID("againstMuonLoose3") > 0.5 &'
    ' tauID("againstElectronLooseMVA3") > 0.5'
)
process.analysisPatTaus.finalCut = cms.string('pt > 20. & abs(eta) < 2.3')

process.cleanPatCandidates.replace ( process.cleanPatTaus, process.cleanPatTaus + process.analysisPatTaus )

#----------------------------------------------------------------------------------------------------
# Make analysisPatMuons and add them to the cleanPatCandidates sequence
#----------------------------------------------------------------------------------------------------

process.analysisPatMuons = process.selectedPatMuonsPFlow.clone()
process.analysisPatMuons.finalCut = cms.string("isGlobalMuon & muonID('GlobalMuonPromptTight') & pt > 10")

process.patDefaultSequencePFlow.replace(process.selectedPatMuonsPFlow , process.selectedPatMuonsPFlow +  process.analysisPatMuons )

#----------------------------------------------------------------------------------------------------
# Make analysisPatElectrons and add them to the cleanPatCandidates sequence
#----------------------------------------------------------------------------------------------------

process.analysisPatElectrons = process.selectedPatElectronsPFlow.clone()
process.analysisPatElectrons.finalCut = cms.string('pt > 10')

process.patDefaultSequencePFlow.replace(process.selectedPatElectronsPFlow , process.selectedPatElectronsPFlow + process.analysisPatElectrons)



process.analysisPatJetsAK5PFchs = process.selectedPatJetsPFlow.clone()
process.analysisPatJetsAK5PFchs.finalCut = cms.string("abs(eta)<2.5 & pt > 10")
process.patDefaultSequencePFlow.replace ( process.selectedPatJetsPFlow, process.selectedPatJetsPFlow + process.analysisPatJetsAK5PFchs)

#----------------------------------------------------------------------------------------------------

switchJetCollection(process,cms.InputTag('ak5CaloJets'),
    doJetID      = True , # Perform jet ID algorithm and store ID info in the jet
    doJTA        = True , # Perform jet track association and determine jet charge
    doBTagging   = True , # Perform b-tagging and store b-tagging info in the jet
    doType1MET   = True , # Store Type1 PFMET information.  Label of resulting PFMET collection is: patMETsAK5Calo
    jetIdLabel   = "ak5",# Which jet ID label should be used?
    jetCorrLabel = ('AK5Calo', ['L1FastJet', 'L2Relative', 'L3Absolute']), # Which jet corrections should be used?
    genJetCollection = cms.InputTag("ak5GenJets") # Which GEN jets should be used?
)


### setting up MET
getattr(process,'patPFMet'+postfix).addGenMET = cms.bool(not isData)
process.patPFMet.addGenMET = cms.bool(not isData)


#----------------------------------------------------------------------------------------------------
# Define the systematic shift correction
#----------------------------------------------------------------------------------------------------

process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")

process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc

#----------------------------------------------------------------------------------------------------
# Use the runMetUncertainties tool here
# See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#MET_Systematics_Tools
#----------------------------------------------------------------------------------------------------

from PhysicsTools.PatUtils.tools.metUncertaintyTools import runMEtUncertainties

runMEtUncertainties(
    process,
    jetCollection           = cms.InputTag('analysisPatJetsAK5PFchs'), 
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

#these flags are false for '+postfix' mets by default, but true for non-postfix ones!                                                                              
getattr(process,'patPFJetMETtype1p2Corr'+postfix).skipEM = cms.bool(False)
getattr(process,'patPFJetMETtype1p2Corr'+postfix).skipMuons = cms.bool(False)



jetSmearFileName='PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'
jetSmearHistogram='pfJetResolutionMCtoDataCorrLUT'

import RecoMET.METProducers.METSigParams_cfi as jetResolutions

process.smearedAnalysisPatJets = cms.EDProducer("SmearedPATJetProducer",
            src = cms.InputTag("analysisPatJetsAK5PFchs"),
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

process.rootTupleMuons.InputTag = cms.InputTag('selectedPatMuonsPFlow')
process.rootTupleElectrons.InputTag = cms.InputTag('selectedPatElectronsPFlow')

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
process.LJFilter.tauLabel  = cms.InputTag("cleanPatTaus")                        
process.LJFilter.muLabel   = cms.InputTag("cleanPatMuons")
process.LJFilter.elecLabel = cms.InputTag("cleanPatElectrons")
process.LJFilter.jetLabel  = cms.InputTag("analysisPatJetsAK5PFchs")
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
    #process.ak5PFchsJetsSequence*
    process.goodOfflinePrimaryVertices * 
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
    process.type0PFMEtCorrection*
    process.patPFMETtype0Corr*
    process.producePFMETCorrections*
    process.pfMEtSysShiftCorrSequence*
    # MET filters (required):
    process.EcalDeadCellTriggerPrimitiveFilter*
    process.EcalDeadCellBoundaryEnergyFilter*
    process.HBHENoiseFilterResultProducer*
    process.trackingFailureFilter*
    process.eeBadScFilter*
    process.ecalLaserCorrFilter*
    # Now the regular PAT default sequence
    getattr(process, "patPF2PATSequence" + postfix) * 
    process.patDefaultSequence*
    process.metUncertaintySequence * 
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
    process.smearedAnalysisPatJets*
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
