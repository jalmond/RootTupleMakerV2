import FWCore.ParameterSet.Config as cms

rootTupleMuons = cms.EDProducer("RootTupleMakerV2_Muons",
    InputTag = cms.InputTag('cleanPatMuons'),
    TriggerEventInputTag = cms.InputTag ('patTriggerEvent'),                                
    Prefix = cms.string('Muon'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10),
    MuonIso = cms.double(0.05),
    MuonID = cms.string('GlobalMuonPromptTight'),
    BeamSpotCorr = cms.bool(True),
    UseCocktailRefits = cms.bool(True),
    VertexInputTag = cms.InputTag('offlinePrimaryVertices'),
    SingleMuonTriggerMatch5 = cms.string ("cleanMuonTriggerMatchHLTSingleMuon5"),
    SingleMuonTriggerMatch8 = cms.string ("cleanMuonTriggerMatchHLTSingleMuon8"),
    SingleMuonTriggerMatch12 = cms.string ("cleanMuonTriggerMatchHLTSingleMuon12"),
    SingleMuonTriggerMatch17= cms.string ("cleanMuonTriggerMatchHLTSingleMuon17"),
    SingleMuonTriggerMatch24 = cms.string ("cleanMuonTriggerMatchHLTSingleMuon24"),
    SingleMuonTriggerMatch = cms.string ("cleanMuonTriggerMatchHLTSingleMuon"),
    DoubleMuonTriggerMatch = cms.string ("cleanMuonTriggerMatchHLTDoubleMuon"),
    SingleIsoMuonTriggerMatch = cms.string ("cleanMuonTriggerMatchHLTSingleIsoMuon"),
    EMuTriggerMatch8     = cms.string ("cleanElectronTriggerMatchHLTSingleElectronemu8"),
    EMuTriggerMatch17     = cms.string ("cleanElectronTriggerMatchHLTSingleElectronemu17")
)

cleanMuonTriggerMatchHLTSingleMuon = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'cleanPatMuons' )
, matched = cms.InputTag( 'patTrigger' )          
, matchedCuts = cms.string( 'type( "TriggerMuon" ) && path ( "HLT_Mu40_eta2p1_v*" )' )
, maxDeltaR = cms.double( 0.1 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

cleanMuonTriggerMatchHLTSingleMuon5 = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'cleanPatMuons' )
, matched = cms.InputTag( 'patTrigger' )
, matchedCuts = cms.string( 'type( "TriggerMuon" ) && path ( "HLT_Mu5_v*" )' )
, maxDeltaR = cms.double( 0.1 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object                                                                     
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)                                      
)

cleanMuonTriggerMatchHLTSingleMuon8 = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'cleanPatMuons' )
, matched = cms.InputTag( 'patTrigger' )
, matchedCuts = cms.string( 'type( "TriggerMuon" ) && path ( "HLT_Mu8_v*" )' )
, maxDeltaR = cms.double( 0.1 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object                          
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above) 
)

cleanMuonTriggerMatchHLTSingleMuon12 = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'cleanPatMuons' )
, matched = cms.InputTag( 'patTrigger' )
, matchedCuts = cms.string( 'type( "TriggerMuon" ) && path ( "HLT_Mu12_v*" )' )
, maxDeltaR = cms.double( 0.1 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object                                                                     
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)                                      
)


cleanMuonTriggerMatchHLTSingleMuon17 = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'cleanPatMuons' )
, matched = cms.InputTag( 'patTrigger' )
, matchedCuts = cms.string( 'type( "TriggerMuon" ) && path ( "HLT_Mu17_v*" )' )
, maxDeltaR = cms.double( 0.1 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object                                                                     
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)                                      
)


cleanMuonTriggerMatchHLTSingleMuon24 = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'cleanPatMuons' )
, matched = cms.InputTag( 'patTrigger' )
, matchedCuts = cms.string( 'type( "TriggerMuon" ) && path ( "HLT_Mu24_v*" )' )
, maxDeltaR = cms.double( 0.1 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object                                                                     
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)                                      
)


cleanMuonTriggerMatchHLTSingleIsoMuon = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'cleanPatMuons' )
, matched = cms.InputTag( 'patTrigger' )          
, matchedCuts = cms.string( 'type( "TriggerMuon" ) && path ( "HLT_IsoMu24_eta2p1_v*" )' )
, maxDeltaR = cms.double( 0.1 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)


cleanMuonTriggerMatchHLTDoubleMuon = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'cleanPatMuons' )
, matched = cms.InputTag( 'patTrigger' )
, matchedCuts = cms.string( 'type( "TriggerMuon" ) && path ( "HLT_Mu17_TkMu8_v*" )' )
, maxDeltaR = cms.double( 0.1 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object                                                                                          
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)                                                           
)


cleanElectronTriggerMatchHLTSingleElectronemu8 = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'cleanPatMuons' )
, matched = cms.InputTag( 'patTrigger' )
, matchedCuts = cms.string( 'type( "TriggerMuon" ) && path( "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")' )
, maxDeltaR = cms.double( 0.1 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object                                                                           
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)                          

)


cleanElectronTriggerMatchHLTSingleElectronemu17 = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'cleanPatMuons' )
, matched = cms.InputTag( 'patTrigger' )
  , matchedCuts = cms.string( 'type( "TriggerMuon" ) && path( "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")')
, maxDeltaR = cms.double( 0.1 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object                                                  
                                                                                                                                                                       
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)                           
                                                                                                                                                                       
)
