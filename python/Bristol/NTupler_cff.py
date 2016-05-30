def setup_ntupler(process, cms, options, includeCA08Jets = False):
    print '=' * 60
    print "Setting up NTupler"
    print '=' * 60
    ######################################################################################################
    ################################## nTuple Configuration ##############################################
    ######################################################################################################


    #calo jets

    #PF2PAT jets
    process.rootTuplePFJets.InputTag = cms.InputTag('smearedAnalysisPatJets2')
    
    #GSF Electrons
    #process.rootTupleElectrons.InputTag = cms.InputTag('selectedPatElectrons')

    #isolated PF Electrons
    #process.rootTuplePFElectrons.InputTag = cms.InputTag('selectedPatElectronsPFlow')

    #muons
    #process.rootTupleMuons.InputTag = cms.InputTag('selectedPatMuons')

    #standard PF muons
    #process.rootTuplePFMuons.InputTag = cms.InputTag('selectedPatMuonsPF')

    print '=' * 60
    print "Setting up NTupler"
    print '=' * 60
    process.rootTupleTree = cms.EDAnalyzer("RootTupleMakerV2_Tree",
        outputCommands=cms.untracked.vstring(
           'drop *',
           # Event information                                                                                                                                                                                              
           'keep *_rootTupleEvent_*_*',
           'keep *_rootTupleEventSelection_*_*',
           # Single objects                                                                                                                                                                                          
           'keep *_rootTuplePFCandidates_*_*',
           'keep *_rootTuplePFJets_*_*',
           'keep *_rootTupleCaloJets_*_*',
           'keep *_rootTupleElectrons_*_*',
           'keep *_rootTuplePFElectrons_*_*',
           'keep *_rootTupleMuons_*_*',
           'keep *_rootTuplePFMuons_*_*',
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
    
    print "NTupler_cff.py: process.rootNTuples"

    process.rootNTuples = cms.Sequence((
 # Event information                                                                                                                                                                                           
    #process.rootTupleEvent+
    #process.rootTupleEventSelection+
    # Single objects                                                                                                                                                                                              
    #process.rootTuplePFCandidates+
    process.rootTuplePFJets
   # process.rootTupleCaloJets+
    #process.rootTupleElectrons+
    #process.rootTuplePFElectrons+
    #process.rootTupleMuons+
    #process.rootTuplePFMuons+
    #process.rootTupleHPSTaus+
    #process.rootTuplePhotons+
    #process.rootTupleVertex+
    # MET objects for analysis                                                                                                                                                                                    
    #process.rootTupleTCMET+
    #process.rootTupleCaloMET+
    #process.rootTupleCaloMETType1Cor+
    #process.rootTuplePFMET+
    #process.rootTuplePFMETType1Cor+
    #process.rootTuplePFMETType01Cor+
    #process.rootTuplePFMETType01XYCor+
    # MET objects for systematics                                                                                                                                                                                 
    #process.rootTuplePFMETType01XYCorUnclusteredUp+
    #process.rootTuplePFMETType01XYCorUnclusteredDown+
    #process.rootTuplePFMETType01XYCorElectronEnUp+
    #process.rootTuplePFMETType01XYCorElectronEnDown+
    #process.rootTuplePFMETType01XYCorMuonEnUp+
    #process.rootTuplePFMETType01XYCorMuonEnDown+
    #process.rootTuplePFMETType01XYCorTauEnUp+
    #process.rootTuplePFMETType01XYCorTauEnDown+
    #process.rootTuplePFMETType01XYCorJetResUp+
    #process.rootTuplePFMETType01XYCorJetResDown+
    #process.rootTuplePFMETType01XYCorJetEnUp+
    #process.rootTuplePFMETType01XYCorJetEnDown+
    ## Trigger objects                                                                                                                                                                                             
    #process.rootTupleTrigger+
    #process.rootTupleTriggerObjects+
    # GEN objects                                                                                                                                                                                                 
    #process.rootTupleGenEventInfo+
    #process.rootTupleGenParticles+
    #process.rootTupleGenJets+
    #process.rootTupleGenElectronsFromWs+
    #process.rootTupleGenElectronsFromZs+
    #process.rootTupleGenMuonsFromWs+
    #process.rootTupleGenMuonsFromZs+
    #process.rootTupleGenTausFromWs+
    #process.rootTupleGenTausFromZs+
    #process.rootTupleGenMETTrue+
    #process.rootTupleGenMETCalo
    )*
                                       process.rootTupleTree)
    print '=' * 60
    print "Finished seeting up NTupler"
    print '=' * 60
