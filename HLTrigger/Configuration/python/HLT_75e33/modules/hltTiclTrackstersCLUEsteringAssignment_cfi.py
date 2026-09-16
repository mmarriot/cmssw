import FWCore.ParameterSet.Config as cms

# Layer-cluster to trackster assignment produced on device by CLUEstering.
hltTiclTrackstersCLUEsteringAssignment = cms.EDProducer("TrackstersCLUEsteringProducer@alpaka",
    layerClusters = cms.VInputTag("hltHgcalSoALayerClustersProducer",
                                  "hltHgcalSoALayerClustersProducerHSci",
                                  "hltHgcalSoALayerClustersProducerHSi"),
    sigmaT = cms.vdouble(0.003, 0.012, 0.006),
    filtered_mask = cms.InputTag("hltFilteredLayerClustersCLUE3DHigh","CLUE3DHigh")
)
