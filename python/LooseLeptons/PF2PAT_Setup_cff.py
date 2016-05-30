vertexCollection = 'goodOfflinePrimaryVertices'
from PhysicsTools.PatAlgos.tools.pfTools import *

def setup_looseLeptons(process, cms, postfix="PF", maxLooseLeptonRelIso=999):
    print '=' * 60
    print "Setting up loose leptons"
    print '=' * 60
    
    # In order to have a coherent semileptonic channel also, add
    # some "loose" leptons to do QCD estimates.
    # Use the good primary vertices everywhere.

    #### create new pfIsolatedMuonsLoosePF objects: will not chanhe input like for electrons since TopProjection cuts on muon pt is ok for analysis
    process.pfIsolatedMuonsLoosePF = process.pfIsolatedMuonsPF.clone(
        isolationCut=cms.double(0.5)
        )
    
    #### create new patMuon collection with looser cuts for analysis
    process.patMuonsLoosePF = process.patMuonsPF.clone(
        pfMuonSource=cms.InputTag("pfIsolatedMuonsLoose" + postfix),
        genParticleMatch=cms.InputTag("muonMatchLoose" + postfix)
        )

    ### run pattools: AdaptMuons here to setup new patmuon collection (This is better than calling function to make sure matching is done correctly)
    process.patMuonsLoosePF.pfMuonSource=cms.InputTag("pfIsolatedMuonsLoose" + postfix)
    process.patMuonsLoosePF.genParticleMatch=cms.InputTag("muonMatchLoose" + postfix)
    process.patMuonsLoosePF.useParticleFlow = True
    process.patMuonsLoosePF.userIsolation   = cms.PSet()
    process.patMuonsLoosePF.isoDeposits = cms.PSet(
        pfChargedHadrons = cms.InputTag("muPFIsoDepositCharged" + postfix),
        pfChargedAll = cms.InputTag("muPFIsoDepositChargedAll" + postfix),
        pfPUChargedHadrons = cms.InputTag("muPFIsoDepositPU" + postfix),
        pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutral" + postfix),
        pfPhotons = cms.InputTag("muPFIsoDepositGamma" + postfix)
        )

    process.patMuonsLoosePF.isolationValues = cms.PSet(
        pfChargedHadrons = cms.InputTag("muPFIsoValueCharged04"+ postfix),
        pfChargedAll = cms.InputTag("muPFIsoValueChargedAll04"+ postfix),
        pfPUChargedHadrons = cms.InputTag("muPFIsoValuePU04" + postfix),
        pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral04" + postfix),
        pfPhotons = cms.InputTag("muPFIsoValueGamma04" + postfix)
        )
    
    ### Create new gen matching collection: This is only needed for muons since electrons do not match with pf electrons
    process.muonMatchLoosePF = process.muonMatchPF.clone(
        src=cms.InputTag("pfIsolatedMuonsLoose" + postfix)
        )
    ### Create selectedPatMuonsLoosePF objects for use in analysis
    process.selectedPatMuonsLoosePF = process.selectedPatMuonsPF.clone(
        src=cms.InputTag("patMuonsLoose" + postfix)
        )

    #### Add pfSelectedElectronsLoosePF to loosen pt cut to 10 GeV for loose leptons
    process.pfSelectedElectronsLoosePF = process.pfSelectedElectronsPF.clone()
    process.pfSelectedElectronsLoosePF.src = cms.InputTag("pfElectronsFromVertex"+postfix)
    process.pfSelectedElectronsLoosePF.cut = 'abs(eta)<2.5 && pt>10. && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits<2'
    
    #process.elPFIsoValueEA03Loose = process.elPFIsoValueEA03.clone()
    #process.elPFIsoValueEA03.pfElectrons = cms.InputTag('pfSelectedElectronsLoose' + postfix)



    #### Add pfIsolatedElectronsLoosePF to loosen isolation cut to 0.5
    process.pfIsolatedElectronsLoosePF = process.pfIsolatedElectronsPF.clone(
        src=cms.InputTag("pfSelectedElectronsLoose"+postfix),
        isolationCut=cms.double(0.5)
        )
    #process.pfIsolatedElectronsLoosePF.doDeltaBetaCorrection = True
    #process.pfIsolatedElectronsLoosePF.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFId"+postfix))
    #process.pfIsolatedElectronsLoosePF.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFId"+postfix),
    #                                                                  cms.InputTag("elPFIsoValueGamma03PFId"+postfix))
    #process.pfIsolatedElectronsLoosePF.deltaBetaIsolationValueMap = 'elPFIsoValueEA03' # new Effective Area edm::ValueMap                                                              
    #process.pfIsolatedElectronsLoosePF.doDeltaBetaCorrection = True                    # not really a 'deltaBeta' correction, but it serves                                            
    #process.pfIsolatedElectronsLoosePF.deltaBetaFactor = -1.0
    #process.pfIsolatedElectronsLoosePF.isolationCut = 0.5


    
    #### Define new patElectronsLoosePF objetcs with looser pfelectron selection
    process.patElectronsLoosePF = process.patElectronsPF.clone(
        pfElectronSource=cms.InputTag("pfIsolatedElectronsLoose" + postfix)
        )
    #adaptPFElectrons(process, process.patElectronsLoosePF, postfix)
    ### run pftool: adaptElectron ourselves
    process.patElectronsLoosePF.useParticleFlow = True
    process.patElectronsLoosePF.pfElectronSource = cms.InputTag("pfIsolatedElectronsLoose" + postfix)
    process.patElectronsLoosePF.userIsolation   = cms.PSet()
    process.patElectronsLoosePF.isoDeposits = cms.PSet(
        pfChargedHadrons = cms.InputTag("elPFIsoDepositCharged" + postfix),
        pfChargedAll = cms.InputTag("elPFIsoDepositChargedAll" + postfix),
        pfPUChargedHadrons = cms.InputTag("elPFIsoDepositPU" + postfix),
        pfNeutralHadrons = cms.InputTag("elPFIsoDepositNeutral" + postfix),
        pfPhotons = cms.InputTag("elPFIsoDepositGamma" + postfix)
        )
    #### change dR cone to 0.3 vs standard cone of 0.4
    process.patElectronsLoosePF.isolationValues = cms.PSet(
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFId"+ postfix),
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFId"+ postfix),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFId" + postfix),
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFId" + postfix),
        pfPhotons = cms.InputTag("elPFIsoValueGamma03PFId" + postfix)
        )
    
    ### comment removal of patelectron object since this is done already
    #removeIfInSequence(process,  "patElectronIsolation",  "patDefaultSequence", postfix)
    
    ### make new selectedPatElectronsLoosePF objects for use in analyis (veto + QCD)
    process.selectedPatElectronsLoosePF = process.selectedPatElectronsPF.clone(
        src=cms.InputTag("patElectronsLoose" + postfix)
        )
    
    #### define loose lepton sequence here to run in submit script
    process.looseLeptonSequence = cms.Sequence(
        process.pfIsolatedMuonsLoosePF + 
        process.muonMatchLoosePF+ 
        process.patMuonsLoosePF + 
        process.selectedPatMuonsLoosePF + 
        process.pfSelectedElectronsLoosePF + 
        #process.elPFIsoValueEA03Loose+
        process.pfIsolatedElectronsLoosePF + 
        process.patElectronsLoosePF + 
        process.selectedPatElectronsLoosePF
        )

    for imod in [process.patMuonsPF,
                 process.patMuonsLoosePF,
                 process.patElectronsPF,
                 process.patElectronsLoosePF] :

        imod.pvSrc = vertexCollection
        imod.embedTrack = True
