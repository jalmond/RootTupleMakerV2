# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

isData=False
###############################
####### Parameters ############
###############################
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

options.register ('useData',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Run this on real data")

options.register ('hltProcess',
                  'HLT',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "HLT process name to use.")
#options.register ('triggerResults',
#                  'HLT',
#                  VarParsing.multiplicity.singleton,
#                  VarParsing.varType.string,
#                  "trigger results to be used.")

options.register ('writeFat',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Output tracks and PF candidates (and GenParticles for MC)")

options.register ('writeSimpleInputs',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Write four-vector and ID of PF candidates")

options.register ('writeGenParticles',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Output GenParticles collection")

options.register ('forceCheckClosestZVertex',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Force the check of the closest z vertex")


options.register ('useSusyFilter',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Use the SUSY event filter")


options.register ('useExtraJetColls',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Write extra jet collections for substructure studies")

options.parseArguments()


if not options.useData :
    inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])

    process.source.fileNames = [
        'file:/afs/cern.ch/work/j/jalmond/FE4C2F81-D0E1-E111-9080-0030487E0A2D.root'
        ]    

else :
    inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
    process.source.fileNames = [
        '/store/mc/Summer12/TTJets_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S7_START52_V5-v1/0000/FEBE99BB-3881-E111-B1F3-003048D42DC8.root'
        ]

    

#process.source.eventsToProcess = cms.untracked.VEventRange( ['1:86747'] )

#process.source.skipEvents = cms.untracked.uint32(17268) 

print options

print 'Running jet corrections: '
print inputJetCorrLabel

import sys


###############################
####### Global Setup ##########
###############################

# 4.2.x or 52x configuration
fileTag = "52x"
if options.useData :
    #process.GlobalTag.globaltag = cms.string( 'GR_R_52_V7::All' )
    process.GlobalTag.globaltag = cms.string( 'GR_P_V42_AN3::All' )
else :
    process.GlobalTag.globaltag = cms.string( 'START53_V22::All' ) 
    #process.GlobalTag.globaltag = cms.string( 'START53_V15::All' )
    #process.GlobalTag.globaltag = cms.string( 'START53_V19::All' )


# require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter = cms.untracked.bool(True),
                                    debugOn = cms.untracked.bool(False),
                                    numtrack = cms.untracked.uint32(10),
                                    thresh = cms.untracked.double(0.25)
                                    )
# HB + HE noise filtering
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
# Modify defaults setting to avoid an over-efficiency in the presence of OFT PU
process.HBHENoiseFilter.minIsolatedNoiseSumE = cms.double(999999.)
process.HBHENoiseFilter.minNumIsolatedNoiseChannels = cms.int32(999999)
process.HBHENoiseFilter.minIsolatedNoiseSumEt = cms.double(999999.)



# switch on PAT trigger 
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger( process, hltProcess=options.hltProcess)




###############################
####### DAF PV's     ##########
###############################

pvSrc = 'offlinePrimaryVertices'

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag("goodOfflinePrimaryVertices"),
                                           minimumNDOF = cms.uint32(5) , # this is > 4
                                           maxAbsZ = cms.double(24), 
                                           maxd0 = cms.double(2) 
                                           )




from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( maxZ = cms.double(24.0),
                                     minNdof = cms.double(5.0) # this is >= 5
                                     ),
    src=cms.InputTag(pvSrc)
    )


###############################
########## Gen Setup ##########
###############################

process.load("RecoJets.Configuration.GenJetParticles_cff")
#from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets


#process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

# add the flavor history
#process.load("PhysicsTools.HepMCCandAlgos.flavorHistoryPaths_cfi")


# prune gen particles
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.prunedGenParticles = cms.EDProducer("GenParticlePruner",
                                            src = cms.InputTag("genParticles"),
                                            select = cms.vstring(
                                                "drop  *"
                                                ,"keep status = 3" #keeps  particles from the hard matrix element
                                                ,"keep (abs(pdgId) >= 11 & abs(pdgId) <= 16) & status = 1" #keeps e/mu and nus with status 1
                                                ,"keep (abs(pdgId)  = 15) & status = 3" #keeps taus
                                                )
                                            )


## process.prunedGenParticles = cms.EDProducer("GenParticlePruner",
##                                             src = cms.InputTag("genParticles"),
##                                             select = cms.vstring(
##                                                 "drop  *"
##                                                 ,"keep++ (abs(pdgId) =6) "
##                                                 )
##                                             )

###############################
#### Jet RECO includes ########
###############################

#from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.CaloJetParameters_cfi import *
#from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
#from RecoJets.JetProducers.CATopJetParameters_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *


###############################
########## PF Setup ###########
###############################

# Default PF2PAT with AK5 jets. Make sure to turn ON the L1fastjet stuff. 
from PhysicsTools.PatAlgos.tools.pfTools import *


# manually adding type0 corrections in the sequence: has to be done before setting up PF2PAT
process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")
process.producePatPFMETCorrections.replace(process.pfCandMETcorr,
                                           process.type0PFMEtCorrection *
                                           process.patPFMETtype0Corr *
                                           process.pfCandMETcorr
                                           )

inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
postfix="PFlow"
if isData:
    inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])


usePF2PAT(process, runPF2PAT=True, jetAlgo='AK5', runOnMC=not options.useData, postfix=postfix, 
          jetCorrections=inputJetCorrLabel,
          pvCollection=cms.InputTag('goodOfflinePrimaryVertices'),
          typeIMetCorrections=True
          )


if not options.forceCheckClosestZVertex :
    process.pfPileUpPFlow.checkClosestZVertex = False

getattr(process,'patType1CorrectedPFMet'+postfix).srcType1Corrections = cms.VInputTag(
     cms.InputTag("patPFJetMETtype1p2Corr"+postfix,"type1"),
     cms.InputTag("patPFMETtype0Corr"+postfix),
     )
getattr(process,'patType1p2CorrectedPFMet'+postfix).srcType1Corrections = cms.VInputTag(
     cms.InputTag("patPFJetMETtype1p2Corr"+postfix,"type1"),
     cms.InputTag("patPFMETtype0Corr"+postfix),
     )



# Turn on the delta-beta corrections
process.pfIsolatedElectronsPFlow.doDeltaBetaCorrection = True
process.pfIsolatedMuonsPFlow.doDeltaBetaCorrection = True

# Set up "loose" leptons. 
process.pfIsolatedMuonsPFlow.isolationCut = cms.double(0.5) 
process.pfIsolatedElectronsPFlow.isolationCut = cms.double(0.5)


# turn to false when running on data
if options.useData :
    removeMCMatching( process, ['All'] )


###############################
###### Electron ID ############
###############################

process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi') 
process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
#Electron ID
process.patElectronsPFlow.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
process.patElectronsPFlow.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0") 
process.patPF2PATSequencePFlow.replace( process.patElectronsPFlow, process.eidMVASequence * process.patElectronsPFlow )


###############################
###### Bare KT 0.6 jets #######
###############################

from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets



# Use the good primary vertices everywhere. 
for imod in [process.patMuonsPFlow,
             process.patElectronsPFlow,
             process.patMuons,
             process.patElectrons] :
    imod.pvSrc = "goodOfflinePrimaryVertices"
    imod.embedTrack = True
    

###############################
### TagInfo and Matching Setup#
###############################

# Do some configuration of the jet substructure things
for jetcoll in  (process.patJetsPFlow,
                 ) :
    if options.useData == False :
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
process.patJetsPFlow.addBTagInfo = True




#################################################
#### Fix the PV collections for the future ######
#################################################
for module in [process.patJetCorrFactors,
               process.patJetCorrFactorsPFlow
               ]:
    module.primaryVertices = "goodOfflinePrimaryVertices"

    


###############################
#### Selections Setup #########
###############################

# AK5 Jets 
process.selectedPatJetsPFlow.cut = cms.string("pt > 10  & abs(rapidity) < 2.4  & chargedEmEnergyFraction < 0.99 & neutralEmEnergyFraction < 0.99 & neutralHadronEnergyFraction < 0.99 & chargedHadronEnergyFraction > 0 & chargedMultiplicity > 0 & numberOfDaughters >1 ")
process.patJetsPFlow.addTagInfos = True
process.patJetsPFlow.tagInfoSources = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfosAODPFlow")
    )
process.patJetsPFlow.userData.userFunctions = cms.vstring( "? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
                                                      "tagInfoSecondaryVertex('secondaryVertex').secondaryVertex(0).p4().mass() : 0")
process.patJetsPFlow.userData.userFunctionLabels = cms.vstring('secvtxMass')



# electrons
process.selectedPatElectrons.cut = cms.string('pt > 10.0 & abs(eta) < 2.5')
process.patElectrons.embedTrack = cms.bool(True)
process.selectedPatElectronsPFlow.cut = cms.string('pt > 10.0 & abs(eta) < 2.5')
process.patElectronsPFlow.embedTrack = cms.bool(True)
# muons
process.patMuons.usePV= cms.bool(True)
process.selectedPatMuons.cut = cms.string('pt > 10.0 & abs(eta) < 2.5')
process.patMuons.embedTrack = cms.bool(True)
process.selectedPatMuonsPFlow.cut = cms.string("pt > 10.0 & abs(eta) < 2.5")
process.patMuonsPFlow.embedTrack = cms.bool(True)


# Apply jet ID to all of the jets upstream. We aren't going to screw around
# with this, most likely. So, we don't really to waste time with it
# at the analysis level. 
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.goodPatJetsPFlow = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                        filterParams = pfJetIDSelector.clone(),
                                        src = cms.InputTag("selectedPatJetsPFlow")
                                        )




## IVF and BCandidate producer for Vbb cross check analysis
process.load('RecoVertex/AdaptiveVertexFinder/inclusiveVertexing_cff')


# let it run

process.patseq = cms.Sequence(
#    process.trigger*     #prova
    process.scrapingVeto*
    process.HBHENoiseFilter*
  
    
    #process.offlinePrimaryVerticesDAF*    
    process.goodOfflinePrimaryVertices*
    process.primaryVertexFilter*
    process.inclusiveVertexing*
    process.genParticlesForJetsNoNu*
    getattr(process,"patPF2PATSequence"+postfix)*
    process.goodPatJetsPFlow*
#    process.flavorHistorySeq*
    process.prunedGenParticles
#    process.miniPFLeptonSequence
    )



if options.useData == True :
    process.patseq.remove( process.genParticlesForJetsNoNu )
    process.patseq.remove( process.genJetParticles )
    process.patseq.remove( process.prunedGenParticles )


if options.writeSimpleInputs :
    process.patseq *= cms.Sequence(process.pfInputs)

if options.useSusyFilter :
    process.patseq.remove( process.HBHENoiseFilter )
    process.load( 'PhysicsTools.HepMCCandAlgos.modelfilter_cfi' )
    process.modelSelector.parameterMins = [500.,    0.] # mstop, mLSP
    process.modelSelector.parameterMaxs = [7000., 200.] # mstop, mLSP
    process.p0 = cms.Path(
        process.modelSelector *
        process.patseq
        )



else :
    process.p0 = cms.Path(
        process.patseq
        )





process.out.SelectEvents.SelectEvents = cms.vstring('p0')

# rename output file
if options.useData :
    if options.writeFat :
        process.out.fileName = cms.untracked.string('ttbsm_' + fileTag + '_data_fat.root')
    else :
        process.out.fileName = cms.untracked.string('ttbsm_' + fileTag + '_data.root')
else :
    if options.writeFat :
        process.out.fileName = cms.untracked.string('ttbsm_' + fileTag + '_mc_fat.root')
    else :
        #process.out.fileName = cms.untracked.string('ttbsm_' + fileTag + '_mc.root')
        process.out.fileName = cms.untracked.string('Step4.root')


# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)


# process all the events
process.maxEvents.input = 100
process.options.wantSummary = False
process.out.dropMetaData = cms.untracked.string("DROPPED")


process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")



process.out.outputCommands = [
    'drop *_cleanPat*_*_*',
    'keep *_selectedPat*_*_*',
    'keep *_goodPat*_*_*',
    'keep patJets_selectedPat*_*_*',
    'keep *_selectedPatJets_*_*',    
    'keep *_patMETs*_*_*',
#    'keep *_offlinePrimaryVertices*_*_*',
    'drop *_kt6PFJets*_*_*',
    'keep *_goodOfflinePrimaryVertices*_*_*',    
    'drop patPFParticles_*_*_*',
#    'drop patTaus_*_*_*',
#    'keep recoPFJets_caPruned*_*_*',
#    'keep recoPFJets_ca*Filtered*_*_*',
#    'keep recoPFJets_caTopTag*_*_*',
    'keep patTriggerObjects_patTriggerPFlow_*_*',
    'keep patTriggerFilters_patTriggerPFlow_*_*',
    'keep patTriggerPaths_patTriggerPFlow_*_*',
    'keep patTriggerEvent_patTriggerEventPFlow_*_*',
    'drop *_cleanPatPhotonsTriggerMatch*_*_*',
    'drop *_cleanPatElectronsTriggerMatch*_*_*',
    'keep *_cleanPatMuonsTriggerMatch*_*_*',
    'drop *_cleanPatTausTriggerMatch*_*_*',
    'drop *_cleanPatJetsTriggerMatch*_*_*',
    'drop *_patMETsTriggerMatch*_*_*',
    'keep double_*_*_PAT',
    'keep *_TriggerResults_*_*',
    'keep edmTriggerResults_*_*_*',
    'keep *_hltTriggerSummaryAOD_*_*',
    'keep *_prunedGenParticles_*_*',
    'drop recoPFCandidates_selectedPatJets*_*_*',
    'drop recoPFCandidates_selectedPatJetsPFlow_*_*',
    'keep CaloTowers_selectedPatJets*_*_*',
    'keep CaloTowers_goodPatJets*_*_*',
    'drop recoBasicJets_*_*_*',
    'keep *_*Lite_*_*',
    'drop recoGenJets_selectedPatJets*_*_*',
    'drop recoGenJets_goodPatJets*_*_*',
    'keep *_*_rho_*',
    'drop *_*PFlowLoose*_*_*',
    'drop patElectrons_*PFlowLoose*_*_*',
    'drop patMuons_*PFlowLoose*_*_*',
#   'keep patTaus_*PFlowLoose*_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_allConversions_*_*',
    'keep *_genParticles_*_*',
    'drop *_ak*Jets*_*_*',
    'drop *_kt*Jets*_*_*',
    'drop *_iterativeCone*_*_*',
    'drop *_pfJetsPFlow_*_*',
    'drop *_selectedPatJetsForMET*_*_*',  #se droppo questi 67 MB per 500 eventi
    'drop *patTaus_*_*_*',
    'drop GenEventInfoProduct_*_*_*',
    'drop recoBaseTagInfosOwned_*_*_*',
    'drop recoConversion_*_*_*'           #se droppo questi arrivo a 31 MB
    #'drop *_goodPatJets*_*_*'            #se droppo questi arrivo a 30 MB
    #'keep recoTracks_generalTracks_*_*'
    ]

if options.useData :
    process.out.outputCommands += ['drop *_MEtoEDMConverter_*_*',
                                   'keep LumiSummary_lumiProducer_*_*'
                                   ]
else :
    process.out.outputCommands += [
        #     'keep recoGenJets_ca8GenJetsNoNu_*_*',
        #   'keep recoGenJets_ak5GenJetsNoNu_*_*',
        #   'keep recoGenJets_ak7GenJetsNoNu_*_*',
        #   'keep recoGenJets_ak8GenJetsNoNu_*_*',
        #   'keep recoGenJets_caFilteredGenJetsNoNu_*_*',
        #   'keep recoGenJets_caPrunedGen_*_*',
        #   'keep *_caTopTagGen_*_*',
                                   'keep GenRunInfoProduct_generator_*_*',
                                   'keep GenEventInfoProduct_generator_*_*',
  #                                 'keep *_flavorHistoryFilter_*_*',
                                   'keep PileupSummaryInfos_*_*_*',
                                      'keep recoGenJets_selectedPatJetsPFlow_*_*',
                                   ]

if options.writeFat :

    process.out.outputCommands += [
        'keep *_pfNoElectron*_*_*',
        'keep recoTracks_generalTracks_*_*',
        'drop recoPFCandidates_selectedPatJets*_*_*',
        'keep recoBaseTagInfosOwned_selectedPatJets*_*_*',
        'keep CaloTowers_selectedPatJets*_*_*'
        ]
if options.writeFat or options.writeGenParticles :
    if options.useData == False :
        process.out.outputCommands += [
            'keep *_genParticles_*_*'
            ]


if options.writeSimpleInputs :
    process.out.outputCommands += [
        'keep *_pfInputs_*_*'
        ]


#open('junk.py','w').write(process.dumpPython())
