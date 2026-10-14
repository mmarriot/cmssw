#ifndef RecoLocalCalo_HGCalRecProducers_plugins_alpaka_HGCalLayerClustersSoAAlgoWrapper_h
#define RecoLocalCalo_HGCalRecProducers_plugins_alpaka_HGCalLayerClustersSoAAlgoWrapper_h

#include <cstdint>

#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "DataFormats/CaloRecHit/interface/alpaka/CaloClusterDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class HGCalLayerClustersSoAAlgoWrapper {
  public:
    // Fills the nLayerClusters layer clusters (energy, number of cells, seed and position) from the rechits and their
    // clustering (the layer cluster and the seed flag of each rechit)
    void run(Queue& queue,
             uint32_t nLayerClusters,
             float thresholdW0,
             float positionDeltaRho2,
             bool isScintillator,
             HGCalSoARecHitsDeviceCollection::ConstView recHits,
             HGCalSoARecHitsExtraDeviceCollection::ConstView recHitsClustering,
             reco::CaloClusterDeviceCollection::View layerClusters) const;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoLocalCalo_HGCalRecProducers_plugins_alpaka_HGCalLayerClustersSoAAlgoWrapper_h
