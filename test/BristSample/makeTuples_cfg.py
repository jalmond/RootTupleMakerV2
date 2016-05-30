from Leptoquarks.RootTupleMakerV2.NTupleTools_cff import *  

from Leptoquarks.RootTupleMakerV2.NTupler_cff import *
setup_ntupler(process, cms, options)

process.GlobalTag.globaltag = 'START53_V27::All'


process.source = cms.Source ("PoolSource",
                             fileNames=cms.untracked.vstring('file:/afs/cern.ch/work/j/jalmond/FE4C2F81-D0E1-E111-9080-0030487E0A2D.root'      ),
                             skipEvents=cms.untracked.uint32(993)                                                                                                                                                                                                )

process.maxEvents.input = 1


process.options.wantSummary = True

process.load("Leptoquarks.LeptonJetFilter.leptonjetfilter_cfi")
#### Shared Muon/Electron/Tau Skim                                                                                                                                                                                
process.LJFilter.tauLabel  = cms.InputTag("selectedPatTaus")
process.LJFilter.muLabel   = cms.InputTag("selectedPatMuons")
process.LJFilter.elecLabel = cms.InputTag("selectedPatElectrons")
process.LJFilter.jetLabel  = cms.InputTag("selectedPatJets")
process.LJFilter.muonsMin = 0
process.LJFilter.muPT     = 10.0
process.LJFilter.electronsMin = 0
process.LJFilter.elecPT       = 15.0
process.LJFilter.tausMin = 0
process.LJFilter.tauPT   = 15.0
process.LJFilter.jetsMin = 0
process.LJFilter.jetPT   = 15.0
process.LJFilter.counteitherleptontype = True
process.LJFilter.customfilterEMuTauJet2012 = False
# -- WARNING :                                      


process.makingNTuples = cms.Path(
    process.pdfWeights * 
    process.patseq *
    getattr(process, "producePatPFMETCorrections" + postfix) * 
    getattr(process, "patMETs" + postfix) * 
    #process.printEventContent * 
    process.rootNTuples
    )

#process.makingNTuples.remove(getattr(process, "producePatPFMETCorrections" + postfix))
#process.makingNTuples.remove(getattr(process, "patMETs" + postfix))


process.out.SelectEvents.SelectEvents = cms.vstring('makingNTuples')
