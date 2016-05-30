import FWCore.ParameterSet.Config as cms

#----------------------------------------------------------------------------------------------------
# PFMET collections
#----------------------------------------------------------------------------------------------------

# Raw PFMET, no jet smearing

rootTuplePFMET = cms.EDProducer("RootTupleMakerV2_MET",
    #InputTag = cms.InputTag('RawMET'),
    InputTag = cms.InputTag('patPFMetPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string(''),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

# Type1 PFMET, with jet smearing

rootTuplePFMETType1Cor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patType1CorrectedPFMetType1OnlyPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1Cor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

# Type0+1 PFMET, with jet smearing

rootTuplePFMETType01Cor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patType1CorrectedPFMetType01OnlyPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01Cor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

# Type0+1+XY PFMET, with jet smearing

rootTuplePFMETType01XYCor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patType1CorrectedPFMetType01XYOnlyPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

#----------------------------------------------------------------------------------------------------
# PFMET systematics collections                                                                
#----------------------------------------------------------------------------------------------------

# Shift unclustered energy up

rootTuplePFMETType01XYCorUnclusteredUp = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patType1CorrectedPFMetType01XYOnlyPFUnclustup'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorUnclusteredUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

# Shift unclustered energy down

rootTuplePFMETType01XYCorUnclusteredDown = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patType1CorrectedPFMetType01XYOnlyPFUnclustdown'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorUnclusteredDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

# Shift unclustered electron energy up

rootTuplePFMETType01XYCorElectronEnUp   = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patType1CorrectedPFMetType01OnlyPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorElectronEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

# Shift unclustered electron energy down

rootTuplePFMETType01XYCorElectronEnDown = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patType1CorrectedPFMetType1OnlyPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorElectronEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

# Shift unclustered muon energy up

rootTuplePFMETType01XYCorMuonEnUp       = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patType1CorrectedPFMetType01OnlyPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorMuonEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

# Shift unclustered muon energy down

rootTuplePFMETType01XYCorMuonEnDown     = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patType1CorrectedPFMetType01XYOnlyPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorMuonEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

# Shift unclustered tau energy up

rootTuplePFMETType01XYCorTauEnUp        = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patType1CorrectedPFMetType01XYOnlyPFUnclustdown'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorTauEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

# Shift unclustered tau energy down

rootTuplePFMETType01XYCorTauEnDown      = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patType1CorrectedPFMetType01XYOnlyPFUnclustup'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorTauEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

# Shift unclustered jet energy resolution shifted up

rootTuplePFMETType01XYCorJetResUp       = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patPFMetPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorJetResUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

# Shift unclustered jet energy resolution shifted down

rootTuplePFMETType01XYCorJetResDown     = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patPFMetPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorJetResDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

# Shift unclustered jet energy shifted up

rootTuplePFMETType01XYCorJetEnUp        = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patPFMetPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorJetEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

# Shift unclustered jet energy shifted down

rootTuplePFMETType01XYCorJetEnDown      = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patPFMetPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorJetEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

