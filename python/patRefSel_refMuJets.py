
# This file contains the Top PAG reference selection for mu + jets analysis.
#

### Muon configuration

# PAT muons
muonsUsePV     = False # use beam spot rather than PV, which is necessary for 'dB' cut
muonEmbedTrack = True  # embedded track needed for muon ID cuts

### Jet configuration

# Jet algorithm
jetAlgo = 'AK5'

# JEC sets
jecSetBase = jetAlgo


### ------------------------- Reference selection -------------------------- ###


# PF2PAT settings
from Leptoquarks.RootTupleMakerV2.patRefSel_PF2PAT import *

# HLT selection
triggerSelectionDataRelVals = 'HLT_IsoMu17_eta2p1_TriCentralPFJet30_v*' # 2012A RelVals
triggerSelectionData        = 'HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v*'
triggerSelectionMC          = 'HLT_*' # not recommended

### Muon selection

# Minimal selection for all muons, also basis for signal and veto muons

muonCut  =     'isPFMuon'                                                                      # general reconstruction property
muonCut += ' && (isGlobalMuon || isTrackerMuon)'                                               # general reconstruction property
muonCut += ' && pt > 10.'                                                                      # transverse momentum
muonCut += ' && abs(eta) < 2.5'                                                                # pseudo-rapisity range
muonCut += ' && (chargedHadronIso+neutralHadronIso+photonIso-0.5*puChargedHadronIso)/pt < 0.2' # relative isolation w/ Delta beta corrections (factor 0.5)

LQmuonCut = 'pt > 10.'
LQmuonCut += ' && abs(eta) < 2.5'  


# Signal muon selection on top of 'muonCut'
signalMuonCut  =     'isPFMuon'                                                                       # general reconstruction property
signalMuonCut += ' && isGlobalMuon'                                                                   # general reconstruction property
signalMuonCut += ' && pt > 26.'                                                                       # transverse momentum
signalMuonCut += ' && abs(eta) < 2.1'                                                                 # pseudo-rapisity range
signalMuonCut += ' && globalTrack.normalizedChi2 < 10.'                                               # muon ID: 'isGlobalMuonPromptTight'
signalMuonCut += ' && track.hitPattern.trackerLayersWithMeasurement > 5'                              # muon ID: 'isGlobalMuonPromptTight'
signalMuonCut += ' && globalTrack.hitPattern.numberOfValidMuonHits > 0'                               # muon ID: 'isGlobalMuonPromptTight'
signalMuonCut += ' && abs(dB) < 0.2'                                                                  # 2-dim impact parameter with respect to beam spot (s. "PAT muon configuration" above)
signalMuonCut += ' && innerTrack.hitPattern.numberOfValidPixelHits > 0'                               # tracker reconstruction
signalMuonCut += ' && numberOfMatchedStations > 1'                                                    # muon chamber reconstruction
signalMuonCut += ' && (chargedHadronIso+neutralHadronIso+photonIso-0.5*puChargedHadronIso)/pt < 0.12' # relative isolation w/ Delta beta corrections (factor 0.5)

muonVertexMaxDZ = 0.5 # DeltaZ between muon vertex and PV

### Jet selection

# Signal jet selection
jetCut  =     'abs(eta) < 2.5'                                        # pseudo-rapisity range
jetCut += ' && numberOfDaughters > 1'                                 # PF jet ID:
jetCut += ' && neutralHadronEnergyFraction < 0.99'                    # PF jet ID:
jetCut += ' && neutralEmEnergyFraction < 0.99'                        # PF jet ID:
jetCut += ' && (chargedEmEnergyFraction < 0.99 || abs(eta) >= 2.4)'   # PF jet ID:
jetCut += ' && (chargedHadronEnergyFraction > 0. || abs(eta) >= 2.4)' # PF jet ID:
jetCut += ' && (chargedMultiplicity > 0 || abs(eta) >= 2.4)'          # PF jet ID:
# varying jet pt threshods
veryLooseJetCut = 'pt > 20.' # transverse momentum (all jets)

LQJetCut = 'pt > 10.' 
LQJetCut += ' && abs(eta) < 2.5' 

looseJetCut     = 'pt > 35.' # transverse momentum (3rd jet, optional)
tightJetCut     = 'pt > 45.' # transverse momentum (leading jets)
vtightJetCut    = 'pt > 55.' # transverse momentum (leading jets)

### Electron selection

# Minimal selection for all electrons, also basis for signal and veto muons
electronCut  =     'ecalDrivenMomentum.pt > 20.'                                                                              # transverse energy
electronCut += ' && abs(ecalDrivenMomentum.eta) < 2.5'                                                                        # pseudo-rapisity range
electronCut += ' && electronID("mvaTrigV0") > 0.5'                                                          # MVA electrons ID
electronCut += ' &&  (chargedHadronIso+max(0.,neutralHadronIso+photonIso-1.0*userIsolation("User1Iso")))/ecalDrivenMomentum.pt < 0.15 '

LQelectronCut  =     'ecalDrivenMomentum.pt > 10.'                                                                              
LQelectronCut +=' && abs(ecalDrivenMomentum.eta) < 2.5'


signalElectronCut  = 'abs(ecalDrivenMomentum.eta) < 2.5'  
signalElectronCut  += '&& ecalDrivenMomentum.pt > 30.' 
signalElectronCut  += '&& ( abs(superCluster.eta) < 1.4442  ||'
signalElectronCut  += 'abs(superCluster.eta) > 1.5660 )'
signalElectronCut  += '&& abs(dB)  <  0.02 '   
signalElectronCut+= ' && electronID("mvaTrigV0") > 0.5'
signalElectronCut+= ' &&gsfTrack.trackerExpectedHitsInner.numberOfHits = 0   '  
signalElectronCut += ' &&  (chargedHadronIso+max(0.,neutralHadronIso+photonIso-1.0*userIsolation("User1Iso")))/ecalDrivenMomentum.pt < 0.10 '

signalElectronCutTight  = 'abs(ecalDrivenMomentum.eta) < 2.5'
signalElectronCutTight  += '&& ecalDrivenMomentum.pt > 30.'
signalElectronCutTight  += '&& ( abs(superCluster.eta) < 1.4442  ||'
signalElectronCutTight  += 'abs(superCluster.eta) > 1.5660 )'
signalElectronCutTight  += '&& abs(dB)  <  0.02 '
signalElectronCutTight  += '&&passConversionVeto()'
signalElectronCutTight+= ' && electronID("mvaTrigV0") > 0.5'
signalElectronCutTight+= ' &&gsfTrack.trackerExpectedHitsInner.numberOfHits = 0   '
signalElectronCutTight += ' &&  (chargedHadronIso+max(0.,neutralHadronIso+photonIso-1.0*userIsolation("User1Iso")))/ecalDrivenMomentum.pt < 0.10 '


bjetCut = 'bDiscriminator("combinedSecondaryVertexBJetTags") > 0.679'


### ------------------------------------------------------------------------ ###

# Trigger object selection
triggerObjectSelectionDataRelVals = 'type("TriggerMuon") && ( path("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v*") )' # 2012A RelVals
triggerObjectSelectionData        = 'type("TriggerMuon") && ( path("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v*") )'
triggerObjectSelectionMC          = 'type("TriggerMuon") && ( path("HLT_*") )' # not recommended
