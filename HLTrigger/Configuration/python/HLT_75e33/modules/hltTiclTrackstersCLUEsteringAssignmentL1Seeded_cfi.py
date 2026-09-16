import FWCore.ParameterSet.Config as cms

# As hltTiclTrackstersCLUEsteringAssignment, for the L1-seeded regions:
hltTiclTrackstersCLUEsteringAssignmentL1Seeded = cms.EDProducer("TrackstersCLUEsteringProducer@alpaka",
    layerClusters = cms.VInputTag("hltHgcalSoALayerClustersProducerL1Seeded",
                                  "hltHgcalSoALayerClustersProducerHSciL1Seeded",
                                  "hltHgcalSoALayerClustersProducerHSiL1Seeded"),
    sigmaT = cms.vdouble(0.003, 0.012, 0.006),
    filtered_mask = cms.InputTag("hltFilteredLayerClustersCLUE3DHighL1Seeded","CLUE3DHigh")
)
