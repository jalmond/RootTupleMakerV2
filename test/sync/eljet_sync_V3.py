import sys

import FWCore.ParameterSet.Config as cms

# setup 'standard' options
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('standard')
options.register('runOnMC', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "decide if run on MC or data")
# parsing command line arguments
if( hasattr(sys, "argv") ):
  #options.parseArguments()
  if(len(sys.argv) > 1):
    print "Parsing command line arguments:"
  for args in sys.argv :
    arg = args.split(',')
    for val in arg:
      val = val.split('=')
      if(len(val)==2):
        print "Setting *", val[0], "* to:", val[1]
        setattr(options,val[0], val[1])


process = cms.Process( 'PAT' )


### ======================================================================== ###
###                                                                          ###
###                                 Constants                                ###
###                            (user job steering)                           ###
###                                                                          ###
### ======================================================================== ###


### Data or MC?
runOnMC = options.runOnMC
print runOnMC

### Switch on/off selection steps

# Step 0a
useTrigger      = False
# Step 0b
useGoodVertex   = True
# Step 1
useGoodMuon     = True
# Step 2
useMuonVeto     = True
# Step 3
useElectronVeto = True
# Step 4a
use1Jet         = True
# Step 4b
use2Jets        = True
# Step 4c (choice depends on trigger)
use3JetsTight   = True
# Step 5
use4Jets        = True

### Trigger matching?
addTriggerMatching = False

### Reference selection

from Leptoquarks.RootTupleMakerV2.patRefSel_refMuJets import *
# Muons general
#muonsUsePV     = False
#muonEmbedTrack = True
# Muons
#muonCut       = ''
#signalMuonCut = ''
#muonVertexMaxDZ = 0.5
# Electrons
#electronCut = ''
# Jets
#jetCut          = ''
#veryLooseJetCut = 'pt > 20.' # transverse momentum (all jets)
#looseJetCut     = 'pt > 35.' # transverse momentum (3rd jet, optional for 'use3JetsLoose = True')
#tightJetCut     = 'pt > 45.' # transverse momentum (leading jets)

# Trigger and trigger object
#triggerSelectionData       = ''
#triggerObjectSelectionData = ''
#triggerSelectionMC       = ''
#triggerObjectSelectionMC = ''

### Particle flow

postfix = 'PF'

# subtract charged hadronic pile-up particles (from wrong PVs)
# effects also JECs
usePFnoPU       = True # before any top projection
usePfIsoLessCHS = True # switch to new PF isolation with L1Fastjet CHS

# other switches for PF top projections (default: all 'True')
useNoMuon     = True # before electron top projection
useNoElectron = True # before jet top projection
useNoJet      = True # before tau top projection
useNoTau      = True # before MET top projection

# cuts used in top projections
# vertices
#pfVertices  = 'goodOfflinePrimaryVertices'
#pfD0Cut     = 0.2
#pfDzCut     = 0.5
# muons
#pfMuonSelectionCut = 'pt > 5.'
useMuonCutBasePF = False # use minimal (veto) muon selection cut on top of 'pfMuonSelectionCut'
#pfMuonIsoConeR03 = False
#pfMuonCombIsoCut = 0.2
# electrons
#pfElectronSelectionCut  = 'pt > 5. && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits < 2'
useElectronCutBasePF  = False # use minimal (veto) electron selection cut on top of 'pfElectronSelectionCut'
#pfElectronIsoConeR03 = True
#pfElectronCombIsoCut  = 0.2

### JEC levels

# levels to be accessible from the jets
# jets are corrected to L3Absolute (MC), L2L3Residual (data) automatically, if enabled here
# and remain uncorrected, if none of these levels is enabled here
useL1FastJet    = True  # needs useL1Offset being off, error otherwise
useL1Offset     = False # needs useL1FastJet being off, error otherwise
useL2Relative   = True
useL3Absolute   = True
useL2L3Residual = True  # takes effect only on data
useL5Flavor     = False
useL7Parton     = False

typeIMetCorrections = True

### Input

# list of input files
inputFiles = ['file:/afs/cern.ch/work/j/jalmond/FE4C2F81-D0E1-E111-9080-0030487E0A2D.root'] # overwritten, if "useRelVals" is 'True'

# maximum number of events
maxEvents = -1 # reduce for testiang

### Conditions



### Output

# output file
outputFile = 'patRefSel_elJets.root'

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

process.load( "Leptoquarks.RootTupleMakerV2.patRefSel_basics_cff" )
process.MessageLogger.cerr.FwkReport.reportEvery = fwkReportEvery
process.options.wantSummary = wantSummary

process.GlobalTag.globaltag = 'START53_V27::All'
  


###
### Input configuration
###


process.load( "Leptoquarks.RootTupleMakerV2.patRefSel_inputModule_cfi" )
process.source.fileNames = inputFiles
#process.source.skipEvents=cms.untracked.uint32(1108)
process.maxEvents.input  = maxEvents

#process.source.eventsToProcess = cms.untracked.VEventRange( '1:58892:17663898', '1:58892:17664085', '1:58893:17664133', '1:58893:17664160', '1:58909:17669045', '1:58911:17669630', '1:58937:17677468', '1:94363:28303059')

###
### Output configuration
###

process.load( "Leptoquarks.RootTupleMakerV2.patRefSel_outputModule_cff" )
# output file name
process.out.fileName = outputFile
# event content
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out.outputCommands += patEventContent
# clear event selection
process.out.SelectEvents.SelectEvents = []


###
### Cleaning and trigger selection configuration
###

### Trigger selection
if runOnMC:
  triggerSelection = triggerSelectionMC
else:
  if useRelVals:
    triggerSelection = triggerSelectionDataRelVals
  else:
    triggerSelection = triggerSelectionData
#from Leptoquarks.RootTupleMakerV2.patRefSel_triggerSelection_cff import triggerResults
#process.step0a = triggerResults.clone(
#  triggerConditions = [ triggerSelection ]
#)

### Good vertex selection
process.load( "Leptoquarks.RootTupleMakerV2.patRefSel_goodVertex_cfi" )
process.step0b = process.goodOfflinePrimaryVertices.clone( filter = True )

### Event cleaning
process.load( 'Leptoquarks.RootTupleMakerV2.patRefSel_eventCleaning_cff' )
process.trackingFailureFilter.VertexSource = cms.InputTag( pfVertices )
process.step0c = process.eventCleaning

process.step0c += process.eventCleaningMC


###
### PAT/PF2PAT configuration
###

process.load( "PhysicsTools.PatAlgos.patSequences_cff" )

### Check JECs

jecLevels = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
if not runOnMC:
  jecLevels = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])

### Switch configuration

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

## pile up corrections
from CommonTools.ParticleFlow.Tools.enablePileUpCorrection import *
#enablePileUpCorrectionInPF2PAT( process, postfix, sequence = "patPF2PATSequence"+postfix)
enablePileUpCorrection( process, postfix, sequence = "patPF2PATSequence"+postfix)

################################################
################################################   ADD TOP PROJECTION CUTS TO PF Objects
################################################

process.pfIsolatedElectronsPF.doDeltaBetaCorrection = True
# Turn on the delta-beta corrections
process.pfSelectedMuonsPF.cut = 'abs(eta)<2.5 && pt>10.'
process.pfIsolatedMuonsPF.doDeltaBetaCorrection = True
process.pfIsolatedMuonsPF.deltaBetaFactor = -0.5
process.pfIsolatedMuonsPF.isolationCut = 0.20

process.load( "EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi" )
process.patElectrons.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0"   )
process.patElectrons.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")

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
                                  process.mvaTrigV0 + process.pfIdentifiedElectronsPF + process.pfSelectedElectronsPF + process.elPFIsoValueEA03 )

process.patElectronsPF.isolationValues.pfPhotons = 'elPFIsoValueGamma03PFId'+postfix
process.patElectronsPF.isolationValues.pfNeutralHadrons = 'elPFIsoValueNeutral03PFId'+postfix
process.patElectronsPF.isolationValues.pfChargedHadrons = 'elPFIsoValueCharged03PFId'+postfix
process.patElectronsPF.isolationValues.pfChargedAll= 'elPFIsoValueChargedAll03PFId'+postfix
process.patElectronsPF.isolationValues.pfPUChargedHadrons = 'elPFIsoValuePU03PFId'+postfix
process.patElectronsPF.isolationValues.user = cms.VInputTag(cms.InputTag("elPFIsoValueEA03"))
process.patElectronsPF.electronIDSources = cms.PSet( mvaTrigV0 = cms.InputTag("mvaTrigV0")) ### This is updated later


# Set up "loose" leptons. 
#process.pfIsolatedMuonsPF.isolationCut = cms.double(0.5) 
#process.pfIsolatedElectronsPF.isolationCut = cms.double(0.5)


if pfElectronIsoConeR03:

  getattr( process, 'patElectrons' + postfix ).isolationValues.pfNeutralHadrons   = cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
  getattr( process, 'patElectrons' + postfix ).isolationValues.pfChargedAll       = cms.InputTag( 'elPFIsoValueChargedAll03PFId' + postfix )
  getattr( process, 'patElectrons' + postfix ).isolationValues.pfPUChargedHadrons = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
  getattr( process, 'patElectrons' + postfix ).isolationValues.pfPhotons          = cms.InputTag( 'elPFIsoValueGamma03PFId' + postfix )
  getattr( process, 'patElectrons' + postfix ).isolationValues.pfChargedHadrons   = cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix )


# Use the good primary vertices everywhere. 
for imod in [process.patMuonsPF,
             process.patElectronsPF,
             process.patMuons,
             process.patElectrons] :
    imod.pvSrc = "goodOfflinePrimaryVertices"
    imod.embedTrack = True
    imod.usePV = False

###### Added after ga1
#process.patJetsPF.addTagInfos = True
#    process.patJetsPF.tagInfoSources = cms.VInputTag(
#        cms.InputTag("secondaryVertexTagInfosAODPFlow")
#        )



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


process.load('Leptoquarks.RootTupleMakerV2.Ntuple_cff')



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
if not runOnMC:
  runOnData( process
           , names = [ 'PFAll' ]
           , postfix = postfix
           )
removeSpecificPATObjects( process
                        , names = [ 'Photons', 'Taus' ]
                        , postfix = postfix
                        ) # includes 'removeCleaning'

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

### Muons

intermediatePatElectrons.src = cms.InputTag( 'selectedPatElectrons' + postfix )
setattr( process, 'intermediatePatElectrons' + postfix, intermediatePatElectrons )

goodPatElectrons.src = cms.InputTag( 'selectedPatElectrons' + postfix )
setattr( process, 'goodPatElectrons' + postfix, goodPatElectrons)
step_3elconv.src = cms.InputTag( 'goodPatElectrons' + postfix )  
setattr( process, 'step_3elconv' + postfix, step_3elconv) 

step1.src = cms.InputTag( 'intermediatePatElectrons' + postfix )
setattr( process, 'step1' + postfix, step1 )
step2el.src = cms.InputTag( 'selectedPatMuons' + postfix )
setattr( process, 'step2el' + postfix, step2el )


### Jets
veryLoosePatJets.src = cms.InputTag( 'selectedPatJets' + postfix )
veryLoosePatJets.cut = veryLooseJetCut
setattr( process, 'veryLoosePatJets' + postfix, veryLoosePatJets )

BTaggedPatJets.src = cms.InputTag( 'veryLoosePatJets' + postfix )
BTaggedPatJets.cut = bjetCut
setattr( process, 'BTaggedPatJets' + postfix,BTaggedPatJets)

loosePatJets.src = cms.InputTag( 'veryLoosePatJets' + postfix )
loosePatJets.cut = looseJetCut
setattr( process, 'loosePatJets' + postfix, loosePatJets )

tightPatJets.src = cms.InputTag( 'loosePatJets' + postfix )
tightPatJets.cut = tightJetCut
setattr( process, 'tightPatJets' + postfix, tightPatJets )

vtightPatJets.src = cms.InputTag( 'tightPatJets' + postfix )
vtightPatJets.cut = vtightJetCut
setattr( process, 'vtightPatJets' + postfix, vtightPatJets )

step4a.src = cms.InputTag( 'vtightPatJets' + postfix )
setattr( process, 'step4a' + postfix, step4a )
step4b.src = cms.InputTag( 'tightPatJets' + postfix )
setattr( process, 'step4b' + postfix, step4b )
step4c.src = cms.InputTag( 'loosePatJets' + postfix )
setattr( process, 'step4c' + postfix, step4c )

step5.src = cms.InputTag( 'veryLoosePatJets' + postfix )
setattr( process, 'step5' + postfix, step5  )

step6.src  = cms.InputTag( 'BTaggedPatJets' + postfix )
setattr( process, 'step6' + postfix, step6  )

### Electrons
step3el.src = cms.InputTag( 'selectedPatElectrons' + postfix )
setattr( process, 'step3el' + postfix, step3el )

process.out.outputCommands.append( 'keep *_goodPatElectrons*_*_*' )
process.out.outputCommands.append( 'keep *_veryLoosePatJets*_*_*' )
process.out.outputCommands.append( 'keep *_loosePatJets*_*_*' )
process.out.outputCommands.append( 'keep *_tightPatJets*_*_*' )
process.out.outputCommands.append( 'keep *_vtightPatJets*_*_*' )
process.out.outputCommands.append( 'keep *_pfMuonsFromVertex*_*_*' )

###
### Selection configuration
###

### Muons

getattr( process, 'patMuons' + postfix ).usePV      = muonsUsePV
getattr( process, 'patMuons' + postfix ).embedTrack = muonEmbedTrack
getattr( process, 'selectedPatMuons' + postfix ).cut = muonCut

### Jets
getattr( process, 'selectedPatJets'  + postfix ).cut = jetCut
getattr( process, 'veryLoosePatJets' + postfix ).cut = veryLooseJetCut
getattr( process, 'loosePatJets'     + postfix ).cut = looseJetCut
getattr( process, 'tightPatJets'     + postfix ).cut = tightJetCut
getattr( process, 'vtightPatJets'     + postfix ).cut = vtightJetCut

### Electrons

getattr( process, 'patElectrons' + postfix ).electronIDSources = electronIDSources
getattr( process, 'intermediatePatElectrons' + postfix ).cut = signalElectronCut
getattr( process, 'goodPatElectrons' + postfix ).cut = signalElectronCutTight
getattr( process, 'selectedPatElectrons' + postfix ).cut = electronCut

process.patConversions = cms.EDProducer("PATConversionProducer",
    electronSource = cms.InputTag('selectedPatElectrons' + postfix)
)


###
### Scheduling
###

# MVA electron ID

process.eidMVASequence = cms.Sequence(
  process.mvaTrigV0
+ process.mvaNonTrigV0
)

# The additional sequence

patAddOnSequence = cms.Sequence(
  getattr( process, 'intermediatePatElectrons' + postfix )
* getattr( process, 'goodPatElectrons'     + postfix )    
* getattr( process, 'veryLoosePatJets'     + postfix )
* getattr( process, 'BTaggedPatJets'     + postfix )
* getattr( process, 'loosePatJets'         + postfix )
* getattr( process, 'tightPatJets'         + postfix )
* getattr( process, 'vtightPatJets'         + postfix )
)
setattr( process, 'patAddOnSequence' + postfix, patAddOnSequence )

# The paths




jetSmearFileName='PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'
jetSmearHistogram='pfJetResolutionMCtoDataCorrLUT'

import RecoMET.METProducers.METSigParams_cfi as jetResolutions

process.smearedAnalysisPatJets = cms.EDProducer("SmearedPATJetProducer",
            src = cms.InputTag('veryLoosePatJets'     + postfix),
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


process.rootTuplePFJets.InputTag = cms.InputTag('smearedAnalysisPatJets')
process.rootTupleElectrons.InputTag = cms.InputTag('selectedPatElectrons' + postfix)
process.rootTupleMuons.InputTag = cms.InputTag('patMuons'+ postfix)
process.rootTupleMuons.InputTagEnUp  = cms.InputTag('patMuons' + postfix)
process.rootTupleMuons.InputTagEnDown  = cms.InputTag('patMuons' + postfix)
process.rootTupleElectrons.InputTagEnUp  =  cms.InputTag('selectedPatElectrons' + postfix)
process.rootTupleElectrons.InputTagEnDown =  cms.InputTag('selectedPatElectrons' + postfix)

process.run_ntuple= cms.Sequence(
      process.rootTuplePFJets+
      process.rootTupleElectrons+
      process.rootTupleMuons
)


process.p = cms.Path(
  process.goodOfflinePrimaryVertices*
  process.eidMVASequence*

  getattr( process, 'patPF2PATSequence' + postfix )*
  getattr( process, 'patAddOnSequence' + postfix )*
  process.patConversions*
  process.patDefaultSequence*
  process.smearedAnalysisPatJets*
  (
   process.rootTuplePFJets+
   process.rootTupleElectrons+
   process.rootTupleMuons
   )
)

if useTrigger:
  process.p += process.step0a
#process.p += process.goodOfflinePrimaryVertices
if useGoodVertex:
  process.p += process.step0b
process.p += process.step0c
#process.p += process.eidMVASequence
#process.p += getattr( process, 'patPF2PATSequence' + postfix )
#process.p += getattr( process, 'patAddOnSequence' + postfix )
#process.p += process.smearedAnalysisPatJets
#process.p += process.run_ntuple
if useGoodMuon:
  process.p += getattr( process, 'step1' + postfix )
if useMuonVeto:
  process.p += getattr( process, 'step2el' + postfix )
if useElectronVeto:
  process.p += getattr( process, 'step3el' + postfix )
process.p += getattr( process, 'step_3elconv' + postfix )
if use1Jet:
  process.p += getattr( process, 'step4a' + postfix )
if use2Jets:
  process.p += getattr( process, 'step4b' + postfix )
if use3JetsTight:
  process.p += getattr( process, 'step4c' + postfix )
if use4Jets:
  process.p += getattr( process, 'step5' + postfix )
  process.p += getattr( process, 'step6' + postfix )

process.out.SelectEvents.SelectEvents.append( 'p' )

