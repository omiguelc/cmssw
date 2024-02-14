import FWCore.ParameterSet.Config as cms

hltEG108EtL1SeededFilter = cms.EDFilter("HLTEgammaEtFilter",
    etcutEB = cms.double(108.0),
    etcutEE = cms.double(9999999.0),
    inputTag = cms.InputTag("hltEgammaCandidatesWrapperL1Seeded"),
    l1EGCand = cms.InputTag("hltEgammaCandidatesL1Seeded"),
    ncandcut = cms.int32(1),
    saveTags = cms.bool(True)
)
