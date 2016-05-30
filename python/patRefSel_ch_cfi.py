import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.patSequences_cff import *
from Leptoquarks.RootTupleMakerV2.patRefSel_common_cfi import *

### Muons

intermediatePatMuons = selectedPatMuons.clone(
  src = cms.InputTag( 'selectedPatMuons' )
, cut = '' # signalMuonCut
)

goodPatMuons = cms.EDProducer(
  "MuonSelectorVertex"
, muonSource   = cms.InputTag( 'intermediatePatMuons' )
, vertexSource = cms.InputTag( 'offlinePrimaryVertices' )
, maxDZ        = cms.double( 999. ) # muonVertexMaxDZ
)

intermediatePatElectrons = selectedPatElectrons.clone(
  src = cms.InputTag( 'selectedPatElectrons' )
, cut = '' # signalElectronCut                                                                                                                                                     
)

goodPatElectrons =selectedPatElectrons.clone(
  src = cms.InputTag( 'selectedPatElectrons' )
, cut = '' # signalElectronCut                                                                                                                                               
)



step1 = cms.EDFilter(
  "PATCandViewCountFilter"
, src = cms.InputTag( 'goodPatMuons' )
, minNumber = cms.uint32( 1 )
, maxNumber = cms.uint32( 1 )
)

step2 = countPatMuons.clone(
  maxNumber = 1 # includes the signal muon
)

step3el = countPatElectrons.clone(
  maxNumber = 1 # includes the signal muon                                                                                                                                     
)
BTaggedPatJets = selectedPatJets.clone(
  src = 'veryLoosePatJets'
, cut = '' # looseJetCut                                                                                                                                                      
)

loosePatJets = selectedPatJets.clone(
  src = 'veryLoosePatJets'
, cut = '' # looseJetCut
)
tightPatJets = selectedPatJets.clone(
  src = 'loosePatJets'
, cut = '' # tightJetCut
)


vtightPatJets = selectedPatJets.clone(
  src = 'tightPatJets'
, cut = '' # tightJetCut                                                                                                                                                      
)

step_3elconv = cms.EDFilter(
  "PATCandViewCountFilter"
, src = cms.InputTag( 'goodPatElectrons' )
, minNumber = cms.uint32( 1 )
, maxNumber = cms.uint32( 1 )
)

step4a = cms.EDFilter(
  "PATCandViewCountFilter"
, src = cms.InputTag( 'tightPatJets' )
, minNumber = cms.uint32( 1 )
, maxNumber = cms.uint32( 999999 )
)
step4b = step4a.clone(
  minNumber = 2
)
step4c = step4a.clone(
  minNumber = 3
)

step5  = step4a.clone(
  src       = 'veryLoosePatJets'
, minNumber = 4
)
step6  = step5.clone(
  src       = 'veryLoosePatJets'
, minNumber = 1
)

veryLoosePatJets = selectedPatJets.clone(
  src = 'selectedPatJets'
, cut = '' # veryLooseJetCut
)

step3 = countPatElectrons.clone( maxNumber = 0 )

step2el = countPatMuons.clone( maxNumber = 0 )
