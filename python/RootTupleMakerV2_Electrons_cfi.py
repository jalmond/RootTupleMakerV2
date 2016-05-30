import FWCore.ParameterSet.Config as cms

#pelectrons=selectedPatElectrons

rootTupleElectrons = cms.EDProducer("RootTupleMakerV2_Electrons",
    TracksInputTag = cms.InputTag('generalTracks'),
    DCSInputTag = cms.InputTag('scalersRawToDigi'),
    Prefix = cms.string('Electron'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10),
    ElectronIso = cms.double(0.1),
    MuonPt = cms.double(10.),
    MuonIso = cms.double(0.05),
    MuonID = cms.string('GlobalMuonPromptTight'),
    VertexInputTag = cms.InputTag('offlinePrimaryVertices'),
    BeamSpotInputTag = cms.InputTag('offlineBeamSpot'),
    ConversionsInputTag = cms.InputTag('allConversions'),
    TriggerEventInputTag = cms.InputTag('patTriggerEvent'),                                    
    LikelihoodInputTag = cms.InputTag('egammaIDLikelihood') ,
    RhoInputTag = cms.InputTag('kt6PFJets','rho'),
    SingleEleTriggerMatch8     = cms.string ("cleanElectronTriggerMatchHLTSingleElectron8"),
    SingleEleTriggerMatch17     = cms.string ("cleanElectronTriggerMatchHLTSingleElectron17"),
    SingleEleTriggerMatchWP80 = cms.string ("cleanElectronTriggerMatchHLTSingleElectronWP80"),
    DoubleEleTriggerMatch     = cms.string ("cleanElectronTriggerMatchHLTDoubleElectron"),
    EMuTriggerMatch8     = cms.string ("cleanElectronTriggerMatchHLTSingleElectronemu8"),
    EMuTriggerMatch17     = cms.string ("cleanElectronTriggerMatchHLTSingleElectronemu17")
                                          
)
rootTupleElectronsLoose = rootTupleElectrons.clone(
    Prefix = cms.string('LooseElectron'),
    storePFIsolation = cms.bool(True)
)


cleanElectronTriggerMatchHLTSingleElectron17 = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'selectedPatElectronsPF' )
, matched = cms.InputTag( 'patTrigger' )          
, matchedCuts = cms.string( 'type( "TriggerElectron" ) && path ( "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*" )' )
, maxDeltaR = cms.double( 0.2 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

cleanElectronTriggerMatchHLTSingleElectronWP80 = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'selectedPatElectronsPF' )
, matched = cms.InputTag( 'patTrigger' )          
, matchedCuts = cms.string( 'type( "TriggerElectron" ) && path( "HLT_Ele27_WP80_v*" )' )
, maxDeltaR = cms.double( 0.2 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

cleanElectronTriggerMatchHLTDoubleElectron = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'selectedPatElectronsPF' )
, matched = cms.InputTag( 'patTrigger' )          
, matchedCuts = cms.string( 'type( "TriggerElectron" ) && path( "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*" )' )
, maxDeltaR = cms.double( 0.2 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

cleanElectronTriggerMatchHLTSingleElectron8 = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'selectedPatElectronsPF' )
, matched = cms.InputTag( 'patTrigger' )
, matchedCuts = cms.string( 'type( "TriggerElectron" ) && path( "HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")' )
, maxDeltaR = cms.double( 0.2 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object                                                                            
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)                                          
)


cleanElectronTriggerMatchHLTSingleElectronemu8 = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'selectedPatElectronsPF' )
, matched = cms.InputTag( 'patTrigger' )
, matchedCuts = cms.string( 'type( "TriggerElectron" ) && path( "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")' )
, maxDeltaR = cms.double( 0.2 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)                                
                                                                                                                                                                       
)


cleanElectronTriggerMatchHLTSingleElectronemu17 = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'selectedPatElectronsPF' )
, matched = cms.InputTag( 'patTrigger' )
  , matchedCuts = cms.string( 'type( "TriggerElectron" ) && path( "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")')
, maxDeltaR = cms.double( 0.2 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object                                                                                  
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)                                               \   

)




# Extra trigger matching (for QCD estimate).  Leave commented for now.' )
# , matched = cms.InputTag( 'patTrigger' )          
# , matchedCuts = cms.string( 'type( "TriggerPhoton" ) && path( "HLT_Photon*_CaloIdVL_v*" )' )
# , maxDeltaR = cms.double( 0.5 )
# , resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
# , resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
# )
# 
# cleanElectronTriggerMatchHLTPhoton135 = cms.EDProducer(
#   "PATTriggerMatcherDRLessByR"
# , src     = cms.InputTag( 'selectedPatElectrons' )
# , matched = cms.InputTag( 'patTrigger' )          
# , matchedCuts = cms.string( 'type( "TriggerPhoton" ) && path( "HLT_Photon135_v*" )' )
# , maxDeltaR = cms.double( 0.5 )
# , resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
# , resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
# )
# 
# cleanElectronTriggerMatchHLTPhoton150 = cms.EDProducer(
#   "PATTriggerMatcherDRLessByR"
# , src     = cms.InputTag( 'selectedPatElectrons' )
# , matched = cms.InputTag( 'patTrigger' )          
# , matchedCuts = cms.string( 'type( "TriggerPhoton" ) && path( "HLT_Photon150_v*" )' )
# , maxDeltaR = cms.double( 0.5 )
# , resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
# , resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
# )

#  LocalWords:  scalersRawToDigi
