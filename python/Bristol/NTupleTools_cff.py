# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

process.load('Leptoquarks.RootTupleMakerV2.Ntuple_cff')

postfix = "PFlow"

###############################
####### Parameters ############
###############################
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

FILETAG = '53X'
GLOBALTAG_DATA = 'START53_V27::All'
GLOBALTAG_MC = 'START53_V27::All'
TEST_DATA_FILE = 'file:/afs/cern.ch/work/j/jalmond/FE4C2F81-D0E1-E111-9080-0030487E0A2D.root'
TEST_MC_FILE = 'file:/afs/cern.ch/work/j/jalmond/FE4C2F81-D0E1-E111-9080-0030487E0A2D.root'

removeTausFromJetCollection = False

options.register ('useData',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Run this on real data")


options.register ('CMSSW',
                  '53X',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "CMSSW version used: 53X (default), 52X or 44X")

options.register ('useGSFelectrons',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Use GSF instead of PF electrons in PAT")

options.register ('writeIsolatedPFLeptons',
                  True,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Output isolated PF leptons")

options.register ('forceCheckClosestZVertex',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Force the check of the closest z vertex")
options.register ('printEventContent',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Outputs the event content at the end of the path")


options.register ('maxLooseLeptonRelIso',
                  5.,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.float,
                  "Maximum (PF)relIso value for leptons to be stored.")


options.parseArguments()

print options


GLOBALTAG_DATA = 'START53_V27::All'
GLOBALTAG_MC = 'START53_V27::All'
FILETAG = '53X'
TEST_DATA_FILE = 'file:/afs/cern.ch/work/j/jalmond/FE4C2F81-D0E1-E111-9080-0030487E0A2D.root'
TEST_MC_FILE = 'file:/afs/cern.ch/work/j/jalmond/FE4C2F81-D0E1-E111-9080-0030487E0A2D.root'

maxLooseLeptonRelIso = options.maxLooseLeptonRelIso

if not options.useData :
    process.source.fileNames = [
            TEST_MC_FILE
            ]
else:
    process.source.fileNames = [
            TEST_DATA_FILE
            ]


###############################
####### Global Setup ##########
###############################
if options.useData :
        process.GlobalTag.globaltag = cms.string(GLOBALTAG_DATA)
else :
        process.GlobalTag.globaltag = cms.string(GLOBALTAG_MC)

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
    
process.goodOfflinePrimaryVertices = cms.EDFilter(
  "PrimaryVertexObjectFilter",
  filterParams=cms.PSet(
                        minNdof=cms.double(4.),
                        maxZ=cms.double(24.),
                        maxRho=cms.double(2.)
                        ),
  filter=cms.bool(True),
  src=cms.InputTag('offlinePrimaryVertices')
)

from Leptoquarks.RootTupleMakerV2.PF2PAT_Setup_cff import *
setup_PF2PAT(process, cms, options, postfix=postfix, removeTausFromJetCollection=removeTausFromJetCollection)
setup_looseLeptons(process, cms, options, postfix=postfix, maxLooseLeptonRelIso=maxLooseLeptonRelIso)


from Leptoquarks.RootTupleMakerV2.MET_Setup_cff import *
setup_MET(process, cms, options, postfix=postfix)

from Leptoquarks.RootTupleMakerV2.Jets_Setup_cff import *

#setup_jets(process, cms, options, postfix=postfix)
inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
process.load("RecoJets.Configuration.GenJetParticles_cff")

from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.CaloJetParameters_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *


from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
for jetcoll in  (process.patJetsPFlow,
                 ) :
    
    jetcoll.addGenPartonMatch = True
    jetcoll.embedGenJetMatch = False
    jetcoll.getJetMCFlavour = True
    jetcoll.addBTagInfo = False
    jetcoll.embedCaloTowers = True
    jetcoll.embedPFCandidates = True

process.patJetsPFlow.addBTagInfo = True

for module in [process.patJetCorrFactors,
               process.patJetCorrFactorsPFlow
               ]:
    module.primaryVertices = "goodOfflinePrimaryVertices"

process.patJetsPFlow.addTagInfos = True
process.patJetsPFlow.tagInfoSources = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfosAODPFlow")
    )
addJetCollection(process,cms.InputTag('ak5PFchsJets'),
    'AK5', 'PFchs',
    doJetID      = True , # Perform jet ID algorithm and store ID info in the jet                                                                                 \
                                                                                                                                                                   
    doJTA        = True , # Perform jet track association and determine jet charge                                                                                \
                                                                                                                                                                   
    doBTagging   = True , # Perform b-tagging and store b-tagging info in the jet                                                                                 \
                                                                                                                                                                   
    doType1MET   = False, # Don't store Type1 PFMET information. This will be done by the runMEtUncertainties tool.                                               \
                                                                                                                                                                   
    jetIdLabel   = "ak5",# Which jet ID label should be used?                                                                                                     \
                                                                                                                                                                   
    jetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute']), # Which jet corrections should be used?                                               \
                                                                                                                                                                   
#    jetCorrLabel = ('AK5PF', ['L1FastJet', 'L2Relative',  'L3Absolute']), # Which jet corrections should be used?                                                \
                                                                                                                                                                   
                 genJetCollection = cms.InputTag("ak5GenJets") # Which GEN jets should be used?                                                                   \
                                                                                                                                                                   
)

process.patJetCorrFactorsAK5PFchs.rho = cms.InputTag("kt6PFchsJets","rho")
process.patJetCorrFactorsAK5PFchs.useRho = cms.bool(True)
process.patJetCorrFactorsAK5PFchs.payload = cms.string('AK5PFchs')

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
                                         skipRawJetPtThreshold = cms.double(10.), # GeV                                                                           \
                                                                                                                                                                   
                                                                                                                                                                  \
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




#electron ID
print "Leptoquarks.RootTupleMakerV2.ElectronID_cff"
from Leptoquarks.RootTupleMakerV2.ElectronID_cff import *
setup_electronID(process, cms)

print "Leptoquarks.RootTupleMakerV2.ObjectSelection_cff"

from Leptoquarks.RootTupleMakerV2.ObjectSelection_cff import *
selectObjects(process, cms)



from PhysicsTools.PatAlgos.tools.trigTools import *

# switch on the trigger matching                                                                                                                                                                                                  
#switchOnTriggerMatching( process, triggerMatchers = [
#        'cleanElectronTriggerMatchHLTSingleElectron8',
#        'cleanElectronTriggerMatchHLTSingleElectron17',
#        'cleanElectronTriggerMatchHLTSingleElectronWP80',
#        'cleanElectronTriggerMatchHLTDoubleElectron',
        # muons                                                                                                                                                                                                                   
#        'cleanMuonTriggerMatchHLTSingleMuon',
#        'cleanMuonTriggerMatchHLTSingleMuon5',
#        'cleanMuonTriggerMatchHLTSingleMuon8',
#        'cleanMuonTriggerMatchHLTSingleMuon12',
#        'cleanMuonTriggerMatchHLTSingleMuon17',
#        'cleanMuonTriggerMatchHLTSingleMuon24',
#        'cleanMuonTriggerMatchHLTDoubleMuon',
#        'cleanMuonTriggerMatchHLTSingleIsoMuon'
#] )

process.load("Leptoquarks.RootTupleMakerV2.ak5PFchsJets_cff")

process.patseq = cms.Sequence(
    process.goodOfflinePrimaryVertices * 
    process.eidCiCSequence*
    process.genParticlesForJetsNoNu * 
    process.ak5PFchsJetsSequence*
    getattr(process, "patPF2PATSequence" + postfix) * 
    process.patDefaultSequence * 
    process.metUncertaintySequence*
    process.smearedAnalysisPatJets*
    process.smearedAnalysisPatJets2
    )

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


if options.useData:
    process.patseq.remove(process.genParticlesForJetsNoNu)



#process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")


process.TFileService = cms.Service("TFileService",
                           fileName=cms.string('ntuple.root')
                           )

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

del process.outpath


