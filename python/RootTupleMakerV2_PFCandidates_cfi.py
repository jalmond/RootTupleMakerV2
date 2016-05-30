import FWCore.ParameterSet.Config as cms

rootTuplePFCandidates = cms.EDProducer("RootTupleMakerV2_PFCandidates",
                                       ReducedPFCandidateInputTag = cms.InputTag('pfCandsNotInJetPF'),
                                       ElectronInputTag = cms.InputTag('selectedPatElectronsPF'),
                                       MuonInputTag = cms.InputTag('selectedPatMuonsPF'),                                    
                                       Prefix = cms.string('PFCand'),
                                       Suffix = cms.string(''),
                                       MaxSize = cms.uint32(10),
                                       DRmatch = cms.double(0.1)
                                       )
