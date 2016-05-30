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

isData=False
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
process.maxEvents.input = 1


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


###############################
########## Gen Setup ##########
###############################

process.load("RecoJets.Configuration.GenJetParticles_cff")
#from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets

#from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.CaloJetParameters_cfi import *
#from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
#from RecoJets.JetProducers.CATopJetParameters_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *


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

# Default PF2PAT with AK5 jets. Make sure to turn ON the L1fastjet stuff. 
from PhysicsTools.PatAlgos.tools.pfTools import *


# manually adding type0 corrections in the sequence: has to be done before setting up PF2PAT
process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")
process.producePatPFMETCorrections.replace(process.pfCandMETcorr,
                                           process.type0PFMEtCorrection *
                                           process.patPFMETtype0Corr *
                                           process.pfCandMETcorr
                                           )

postfix=""


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
process.pfIsolatedElectrons.doDeltaBetaCorrection = True
process.pfIsolatedMuons.doDeltaBetaCorrection = True


# turn to false when running on data
if isData :
    removeMCMatching( process, ['All'] )


###############################
###### Electron ID ############
###############################




###############################
###### Bare KT 0.6 jets #######
###############################

from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets


# Use the good primary vertices everywhere. 
    

###############################
### TagInfo and Matching Setup#
###############################

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


process.selectedPatJets.cut = cms.string("pt > 10  & abs(rapidity) < 2.4  & chargedEmEnergyFraction < 0.99 & neutralEmEnergyFraction < 0.99 & neutralHadronEnergyFraction < 0.99 & chargedHadronEnergyFraction > 0 & chargedMultiplicity > 0 & numberOfDaughters >1 ")
process.patJets.addTagInfos = True
process.patJets.tagInfoSources = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfosAOD")
    )
process.patJets.userData.userFunctions = cms.vstring( "? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
                                                      "tagInfoSecondaryVertex('secondaryVertex').secondaryVertex(0).p4().mass() : 0")
process.patJets.userData.userFunctionLabels = cms.vstring('secvtxMass')





# Apply jet ID to all of the jets upstream. We aren't going to screw around
# with this, most likely. So, we don't really to waste time with it
# at the analysis level. 
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.goodPatJets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                        filterParams = pfJetIDSelector.clone(),
                                        src = cms.InputTag("selectedPatJets")
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


jetSmearFileName='PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'
jetSmearHistogram='pfJetResolutionMCtoDataCorrLUT'

import RecoMET.METProducers.METSigParams_cfi as jetResolutions

process.smearedAnalysisPatJets = cms.EDProducer("SmearedPATJetProducer",
            src = cms.InputTag("goodPatJets"),
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

process.smearedAnalysisPatJets2 = cms.EDProducer("SmearedPATJetProducer",
            src = cms.InputTag("cleanPatJetsAK5PFchs"),
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
process.rootTuplePFJets.InputTag = cms.InputTag('smearedAnalysisPatJets2')
process.rootTuplePFJets.InputTagSmearedUp = cms.InputTag('smearedAnalysisPatJets2')
process.rootTuplePFJets.InputTagSmearedDown = cms.InputTag('smearedAnalysisPatJets2')
process.rootTuplePFJets.InputTagScaledUp    = cms.InputTag('smearedAnalysisPatJets2')
process.rootTuplePFJets.InputTagScaledDown  =  cms.InputTag('smearedAnalysisPatJets2')

# let it run
process.load("Leptoquarks.RootTupleMakerV2.leptonGenMatching_cfi")
process.patDefaultSequence.replace( process.electronMatch, process.elMatch )


process.load("Leptoquarks.LeptonJetFilter.leptonjetfilter_cfi")

#### Shared Muon/Electron/Tau Skim                                                                                                                                                 
process.LJFilter.tauLabel  = cms.InputTag("cleanPatTaus")
process.LJFilter.muLabel   = cms.InputTag("cleanPatMuons")
process.LJFilter.elecLabel = cms.InputTag("cleanPatElectrons")
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
    process.genParticlesForJetsNoNu*
    process.ak5PFchsJetsSequence*
    getattr(process,"patPF2PATSequence"+postfix)*
    process.goodPatJets*
    process.prunedGenParticles*
    process.patDefaultSequence*
    process.smearedAnalysisPatJets*
    process.smearedAnalysisPatJets2*
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
