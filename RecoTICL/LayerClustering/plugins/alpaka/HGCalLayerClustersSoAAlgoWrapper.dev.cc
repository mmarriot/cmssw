#include <algorithm>
#include <cmath>
#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "HGCalLayerClustersSoAAlgoWrapper.h"
#include "ConstantsForClusters.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;
  using namespace hgcal::constants;

  namespace {

    // Scintillator (BH) rechits store (eta, phi) in (dim1, dim2). Recover the
    // cartesian transverse position from them and the cell z.
    ALPAKA_FN_ACC ALPAKA_FN_INLINE void etaPhiZToXY(
        const float eta, const float phi, const float z, float& x, float& y) {
      const float r = z / std::sinh(eta);
      x = r * std::cos(phi);
      y = r * std::sin(phi);
    }

    // Per rechit -> per layer cluster: the sum of the energies, the number of cells and the DetId of the seed
    struct KernelLayerClusterEnergy {
      ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                    HGCalSoARecHitsDeviceCollection::ConstView recHits,
                                    HGCalSoARecHitsExtraDeviceCollection::ConstView recHitsClustering,
                                    reco::CaloClusterDeviceCollection::View layerClusters) const {
        for (auto i : uniform_elements(acc, recHits.metadata().size())) {
          const int c = recHitsClustering[i].clusterIndex();
          if (c == kInvalidCluster) {
            continue;
          }
          alpaka::atomicAdd(acc, &layerClusters.energy()[c].energy(), recHits[i].energy());
          alpaka::atomicAdd(acc, &layerClusters.position()[c].cells(), 1);
          if (recHitsClustering[i].isSeed() == 1) {
            layerClusters.indexes()[c].seedID() = recHits[i].detid();
          }
        }
      }
    };

    // Per rechit -> per layer cluster: the sum of the energies and the index of the most energetic rechit
    struct KernelLayerClusterMostEnergeticHit {
      ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                    HGCalSoARecHitsDeviceCollection::ConstView recHits,
                                    HGCalSoARecHitsExtraDeviceCollection::ConstView recHitsClustering,
                                    float* energySum,
                                    int* maxEnergyHit) const {
        for (int32_t i : uniform_elements(acc, recHits.metadata().size())) {
          const int c = recHitsClustering[i].clusterIndex();
          if (c == kInvalidCluster) {
            continue;
          }
          alpaka::atomicAdd(acc, &energySum[c], recHits[i].energy());

          // Replace the most energetic rechit of the cluster with this one while this one has more energy. If
          // another thread replaced it in the meantime, the compare-and-swap fails and returns the new one.
          int current = maxEnergyHit[c];
          float currentEnergy = (current == kInvalidIndex) ? 0.f : recHits[current].energy();
          while (recHits[i].energy() > currentEnergy) {
            const int previous = alpaka::atomicCas(acc, &maxEnergyHit[c], current, i);
            if (previous == i) {
              break;
            }
            current = previous;
            currentEnergy = (current == kInvalidIndex) ? 0.f : recHits[current].energy();
          }
        }
      }
    };

    // Per rechit -> per layer cluster: the weighted sums of the transverse positions and the sum of the weights.
    // The scintillator rechits are weighted by their energy, the silicon rechits within 1.3 cm (1+6 cells, for all
    // thicknesses) of the most energetic one by max(thresholdW0 + log(energy / cluster energy), 0).
    struct KernelLayerClusterPositionSums {
      ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                    float thresholdW0,
                                    float positionDeltaRho2,
                                    bool isScintillator,
                                    HGCalSoARecHitsDeviceCollection::ConstView recHits,
                                    HGCalSoARecHitsExtraDeviceCollection::ConstView recHitsClustering,
                                    reco::CaloClusterDeviceCollection::View layerClusters,
                                    const float* energySum,
                                    float* positionWeightSum,
                                    const int* maxEnergyHit) const {
        for (auto i : uniform_elements(acc, recHits.metadata().size())) {
          const int c = recHitsClustering[i].clusterIndex();
          if (c == kInvalidCluster) {
            continue;
          }

          if (isScintillator) {
            const float energy = recHits[i].energy();
            float x;
            float y;
            etaPhiZToXY(recHits[i].dim1(), recHits[i].dim2(), recHits[i].dim3(), x, y);
            alpaka::atomicAdd(acc, &layerClusters.position()[c].x(), x * energy);
            alpaka::atomicAdd(acc, &layerClusters.position()[c].y(), y * energy);
            alpaka::atomicAdd(acc, &positionWeightSum[c], energy);
            continue;
          }

          const int maxHit = maxEnergyHit[c];
          const float d1 = recHits[i].dim1() - recHits[maxHit].dim1();
          const float d2 = recHits[i].dim2() - recHits[maxHit].dim2();
          if (std::fmaf(d1, d1, d2 * d2) > positionDeltaRho2) {
            continue;
          }
          const float weight = std::max(thresholdW0 + std::log(recHits[i].energy() / energySum[c]), 0.f);
          alpaka::atomicAdd(acc, &layerClusters.position()[c].x(), recHits[i].dim1() * weight);
          alpaka::atomicAdd(acc, &layerClusters.position()[c].y(), recHits[i].dim2() * weight);
          alpaka::atomicAdd(acc, &positionWeightSum[c], weight);
        }
      }
    };

    // Per layer cluster: the transverse position from the weighted sums, or the one of its most energetic rechit if
    // no rechit has a weight, and z and the layer of its most energetic rechit
    struct KernelLayerClusterPosition {
      ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                    bool isScintillator,
                                    HGCalSoARecHitsDeviceCollection::ConstView recHits,
                                    reco::CaloClusterDeviceCollection::View layerClusters,
                                    const float* positionWeightSum,
                                    const int* maxEnergyHit) const {
        for (auto c : uniform_elements(acc, layerClusters.position().metadata().size())) {
          const int maxHit = maxEnergyHit[c];
          auto position = layerClusters.position()[c];
          if (positionWeightSum[c] > 0.f) {
            const float invWeightSum = 1.f / positionWeightSum[c];
            position.x() *= invWeightSum;
            position.y() *= invWeightSum;
          } else if (isScintillator) {
            float x;
            float y;
            etaPhiZToXY(recHits[maxHit].dim1(), recHits[maxHit].dim2(), recHits[maxHit].dim3(), x, y);
            position.x() = x;
            position.y() = y;
          } else {
            position.x() = recHits[maxHit].dim1();
            position.y() = recHits[maxHit].dim2();
          }
          position.z() = recHits[maxHit].dim3();
          position.layer() = recHits[maxHit].layer();
        }
      }
    };

  }  // namespace

  void HGCalLayerClustersSoAAlgoWrapper::run(Queue& queue,
                                             uint32_t nLayerClusters,
                                             float thresholdW0,
                                             float positionDeltaRho2,
                                             bool isScintillator,
                                             HGCalSoARecHitsDeviceCollection::ConstView recHits,
                                             HGCalSoARecHitsExtraDeviceCollection::ConstView recHitsClustering,
                                             reco::CaloClusterDeviceCollection::View layerClusters) const {
    if (nLayerClusters == 0) {
      return;
    }

    auto x = make_device_view<float>(queue, layerClusters.position().x());
    alpaka::memset(queue, x, 0x0);
    auto y = make_device_view<float>(queue, layerClusters.position().y());
    alpaka::memset(queue, y, 0x0);
    auto cells = make_device_view<int>(queue, layerClusters.position().cells());
    alpaka::memset(queue, cells, 0x0);
    auto energy = make_device_view<float>(queue, layerClusters.energy().energy());
    alpaka::memset(queue, energy, 0x0);
    // per layer cluster: the sum of the energies, the sum of the position weights and the most energetic rechit
    auto energySum = make_device_buffer<float[]>(queue, nLayerClusters);
    alpaka::memset(queue, energySum, 0x0);
    auto positionWeightSum = make_device_buffer<float[]>(queue, nLayerClusters);
    alpaka::memset(queue, positionWeightSum, 0x0);
    auto maxEnergyHit = make_device_buffer<int[]>(queue, nLayerClusters);
    alpaka::memset(queue, maxEnergyHit, kInvalidIndexByte);

    constexpr uint32_t blockSize = 64;
    const auto workDivRecHits = make_workdiv<Acc1D>(divide_up_by(recHits.metadata().size(), blockSize), blockSize);
    const auto workDivLayerClusters = make_workdiv<Acc1D>(divide_up_by(nLayerClusters, blockSize), blockSize);
    alpaka::exec<Acc1D>(queue, workDivRecHits, KernelLayerClusterEnergy{}, recHits, recHitsClustering, layerClusters);
    alpaka::exec<Acc1D>(queue,
                        workDivRecHits,
                        KernelLayerClusterMostEnergeticHit{},
                        recHits,
                        recHitsClustering,
                        energySum.data(),
                        maxEnergyHit.data());
    alpaka::exec<Acc1D>(queue,
                        workDivRecHits,
                        KernelLayerClusterPositionSums{},
                        thresholdW0,
                        positionDeltaRho2,
                        isScintillator,
                        recHits,
                        recHitsClustering,
                        layerClusters,
                        energySum.data(),
                        positionWeightSum.data(),
                        maxEnergyHit.data());
    alpaka::exec<Acc1D>(queue,
                        workDivLayerClusters,
                        KernelLayerClusterPosition{},
                        isScintillator,
                        recHits,
                        layerClusters,
                        positionWeightSum.data(),
                        maxEnergyHit.data());
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
