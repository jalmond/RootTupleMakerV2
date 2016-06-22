import FWCore.ParameterSet.Config as cms

rootTuplePFJets = cms.EDProducer("RootTupleMakerV2_PFJets",
    InputTag = cms.InputTag('selectedPatJetsAK5PF'),
    # InputTagL1Offset    = cms.InputTag('selectedPatJetsAK5PFL1Offset'),
    InputTagSmearedUp   = cms.InputTag('smearedPatJetsAK5PFResUp'),                                 
    InputTagSmearedDown = cms.InputTag('smearedPatJetsAK5PFResDown'),                                 
    InputTagScaledUp    = cms.InputTag('shiftedPatJetsAK5PFEnUpForCorrMEt'),                                 
    InputTagScaledDown  = cms.InputTag('shiftedPatJetsAK5PFEnDownForCorrMEt'),                                 
    Prefix = cms.string('PFJet'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(30),
    JECUncertainty = cms.string('AK5PF'),
    ReadJECuncertainty = cms.bool(True),
    VertexInputTag = cms.InputTag('offlinePrimaryVertices')
)
