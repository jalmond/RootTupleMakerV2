import FWCore.ParameterSet.Config as cms

rootTuplePFJets = cms.EDProducer("RootTupleMakerV2_PFJets",
    InputTag = cms.InputTag('selectedPatJetsAK5PFchs'),
    # InputTagL1Offset    = cms.InputTag('selectedPatJetsAK5PFL1Offset'),
    InputTagSmearedUp   = cms.InputTag('smearedPatJetsAK5PFchsResUp'),                                 
    InputTagSmearedDown = cms.InputTag('smearedPatJetsAK5PFchsResDown'),                                 
    InputTagScaledUp    = cms.InputTag('shiftedPatJetsAK5PFchsEnUpForCorrMEt'),                                 
    InputTagScaledDown  = cms.InputTag('shiftedPatJetsAK5PFchsEnDownForCorrMEt'),                                 
    Prefix = cms.string('PFJet'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(30),
    JECUncertainty = cms.string('AK5PFchs'),
    ReadJECuncertainty = cms.bool(True),
    VertexInputTag = cms.InputTag('offlinePrimaryVertices')
)
