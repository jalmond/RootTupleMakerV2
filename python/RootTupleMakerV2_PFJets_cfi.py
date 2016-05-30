import FWCore.ParameterSet.Config as cms

rootTuplePFJets = cms.EDProducer("RootTupleMakerV2_PFJets",
    InputTag = cms.InputTag('selectedPatJetsAK5PF'),
    # InputTagL1Offset    = cms.InputTag('selectedPatJetsAK5PFL1Offset'),
    FastJetForJECInputTag       =cms.InputTag('kt6PFchsJets','rho'),
                                 #cms.InputTag('fixedGridRhoAll'),
                                 #
    InputTagSmearedUp   = cms.InputTag('smearedPatJetsAK5PFresUp'),                                 
    InputTagSmearedDown = cms.InputTag('smearedPatJetsAK5PFresDown'),                                 
    InputTagScaledUp    = cms.InputTag('shiftedPatJetsAK5PFenUpForCorrMEt'),                                 
    InputTagScaledDown  = cms.InputTag('shiftedPatJetsAK5PFenDownForCorrMEt'),                                 
    JECUncertaintyFile      = cms.InputTag('Leptoquarks/RootTupleMakerV2/data/Summer13_V5_DATA_UncertaintySources_AK5PF.txt'),
    Prefix = cms.string('PFJet'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(30),
    JECUncertainty = cms.string('AK5PF'),
    ReadJECuncertainty = cms.bool(True),
    VertexInputTag = cms.InputTag('offlinePrimaryVertices')
)
