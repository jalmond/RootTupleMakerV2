def selectObjects(process, cms):
    print '=' * 60
    print "Setting up selected objects"
    print '=' * 60
    ###############################
    #### Selections Setup #########
    ###############################
    
    ### AK5 Jets
    #calo jets
    process.selectedPatJets.cut = cms.string("pt > 20 & abs(eta) < 2.5")
    #PF jets
    process.selectedPatJetsPFlow.cut = cms.string("pt > 20 & abs(eta) < 2.5")

    process.patJetsPFlow.addTagInfos = True
    process.patJetsPFlow.tagInfoSources = cms.VInputTag(
        cms.InputTag("secondaryVertexTagInfosAODPFlow")
        )
    
    
    # electrons
    process.selectedPatElectrons.cut = cms.string('pt > 10.0 & abs(eta) < 2.5')
    process.patElectrons.embedTrack = cms.bool(True)
    process.selectedPatElectronsPFlow.cut = cms.string('pt > 10.0 & abs(eta) < 2.5')
    process.patElectronsPFlow.embedTrack = cms.bool(True)
    # muons
    process.selectedPatMuons.cut = cms.string('pt > 10.0 & abs(eta) < 2.5')
    process.patMuons.embedTrack = cms.bool(True)
    process.selectedPatMuonsPFlow.cut = cms.string("pt > 10.0 & abs(eta) < 2.5")
    process.patMuonsPFlow.embedTrack = cms.bool(True)
    # taus
    process.selectedPatTausPFlow.cut = cms.string("pt > 10.0 & abs(eta) < 3")
    process.selectedPatTaus.cut = cms.string("pt > 10.0 & abs(eta) < 3")
    process.patTausPFlow.isoDeposits = cms.PSet()
    process.patTaus.isoDeposits = cms.PSet()
    # photons
    process.patPhotonsPFlow.isoDeposits = cms.PSet()
    process.patPhotons.isoDeposits = cms.PSet()
    
    
    # Apply jet ID to all of the jets upstream. We aren't going to screw around
    # with this, most likely. So, we don't really to waste time with it
    # at the analysis level. 
    #default is version=PURE09, quality=loose ID
    #https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Documentation

    process.selectedPatJetsForMETtype1p2Corr.src = cms.InputTag('selectedPatJetsPFlow')
    process.selectedPatJetsForMETtype2Corr.src = cms.InputTag('selectedPatJetsPFlow')
