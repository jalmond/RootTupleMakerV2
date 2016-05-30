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

isData=True
###############################
####### Parameters ############
###############################

inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
if not isData:
    inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
    process.source.fileNames = [
        'root://xrootd.unl.edu//store/data/Run2012A/SingleMu/AOD/22Jan2013-v1/20000/002F5062-346F-E211-BF00-1CC1DE04DF20.root'
        ]



else :
    inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
process.source.fileNames = [
  'root://xrootd.unl.edu//store/data/Run2012A/SingleMu/AOD/22Jan2013-v1/20000/002F5062-346F-E211-BF00-1CC1DE04DF20.root'
]


process.maxEvents.input = 1000


print 'Running jet corrections: '
print inputJetCorrLabel

import sys


###############################
####### Global Setup ##########
###############################

fileTag = "53x"
if isData:
    process.GlobalTag.globaltag = cms.string( 'GR_P_V42_AN3::All' )
else :
    process.GlobalTag.globaltag = cms.string( 'START53_V27::All' ) 


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


# Default PF2PAT with AK5 jets. Make sure to turn ON the L1fastjet stuff. 
from PhysicsTools.PatAlgos.tools.pfTools import *


# manually adding type0 corrections in the sequence: has to be done before setting up PF2PAT
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


# turn to false when running on data
if isData :
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
process.patJetsPFlow.addTagInfos = True
process.patJetsPFlow.tagInfoSources = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfosAODPFlow")
    )


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



from PhysicsTools.PatAlgos.tools.jetTools import *



process.load("Leptoquarks.RootTupleMakerV2.ak5PFchsJets_cff")
#cms.Sequence() += process.ak5PFchsJetsSequence                                                                                                                                    

addJetCollection(process,cms.InputTag('ak5PFchsJets'),
    'AK5', 'PFchs',
    doJetID      = True , # Perform jet ID algorithm and store ID info in the jet                                                                                                  
    doJTA        = True , # Perform jet track association and determine jet charge                                                                                                 
    doBTagging   = True , # Perform b-tagging and store b-tagging info in the jet                                                                                                  
    doType1MET   = False, # Don't store Type1 PFMET information. This will be done by the runMEtUncertainties tool.                                                                
    jetIdLabel   = "ak5",# Which jet ID label should be used?                                                                                                                      
    jetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute']), # Which jet corrections should be used?                                                                
#    jetCorrLabel = ('AK5PF', ['L1FastJet', 'L2Relative',  'L3Absolute']), # Which jet corrections should be used?                                                                 
                 genJetCollection = cms.InputTag("ak5GenJets") # Which GEN jets should be used?                                                                                                 
)

process.patJetCorrFactorsAK5PFchs.rho = cms.InputTag("kt6PFchsJets","rho")
process.patJetCorrFactorsAK5PFchs.useRho = cms.bool(True)
process.patJetCorrFactorsAK5PFchs.payload = cms.string('AK5PFchs')

# top projections in PF2PAT:
getattr(process,"pfNoPileUp"+postfix).enable = True
getattr(process,"pfNoMuon"+postfix).enable = False #do old style cleaning
getattr(process,"pfNoElectron"+postfix).enable = False #do old style cleaning
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True



## IVF and BCandidate producer for Vbb cross check analysis
process.load('RecoVertex/AdaptiveVertexFinder/inclusiveVertexing_cff')


#----------------------------------------------------------------------------------------------------                                                                              
# Load our RootTupleMakerV2 modules                                                                                                                                                
#----------------------------------------------------------------------------------------------------                                                                              

process.load('Leptoquarks.RootTupleMakerV2.Ntuple_cff')
process.load("Leptoquarks.RootTupleMakerV2.metFilters_cfi")


#----------------------------------------------------------------------------------------------------                                                                              
# This is MC, so analyze the smeared PFJets by default                                                                                                                             
#----------------------------------------------------------------------------------------------------                                                                              
process.rootTuplePFJets.InputTag = cms.InputTag('selectedPatJetsPFlow')
process.rootTuplePFJets.InputTagSmearedUp = cms.InputTag('selectedPatJetsPFlow')
process.rootTuplePFJets.InputTagSmearedDown = cms.InputTag('selectedPatJetsPFlow')
process.rootTuplePFJets.InputTagScaledUp    = cms.InputTag('selectedPatJetsPFlow')
process.rootTuplePFJets.InputTagScaledDown  =  cms.InputTag('selectedPatJetsPFlow')

# let it run


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


process.p = cms.Path(
#    process.trigger*     #prova
    #process.offlinePrimaryVerticesDAF*    
    #process.genTausFromWs*
    #process.genMuonsFromWs*
    #process.genElectronsFromWs*
    #process.genTausFromZs*
    #process.genMuonsFromZs*
    #process.genElectronsFromZs*
    #process.goodOfflinePrimaryVertices*
    process.goodOfflinePrimaryVertices*
    process.primaryVertexFilter*
    process.inclusiveVertexing*
    process.ak5PFchsJetsSequence*
    getattr(process,"patPF2PATSequence"+postfix)*
    process.goodPatJetsPFlow*
    process.patDefaultSequence*
    (
      process.rootTuplePFJets
     )
    #process.rootTupleTree
    #process.dump                                                                                                                                                                  
)



# Delete predefined Endpath (needed for running with CRAB)                                                                                                                         
del process.out
del process.outpath


#----------------------------------------------------------------------------------------------------                                                                              
# Run the path                                                                                                                                                                     
#----------------------------------------------------------------------------------------------------                                                                              

process.schedule = cms.Schedule(process.p)
