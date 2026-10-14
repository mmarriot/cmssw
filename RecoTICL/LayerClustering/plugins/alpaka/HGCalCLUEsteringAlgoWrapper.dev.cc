#include <cmath>
#include <cstdint>
#include <span>

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "HGCalCLUEsteringAlgoWrapper.h"

#include "CLUEstering/core/detail/defines.hpp"
#include "CLUEstering/core/Clusterer.hpp"
#include "CLUEstering/data_structures/PointsDevice.hpp"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;

  namespace {

    struct KernelSetSeeds {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int32_t* seeds,
                                    HGCalSoARecHitsExtraDeviceCollection::View outputs,
                                    uint32_t nseeds) const {
        for (auto k : uniform_elements(acc, nseeds)) {
          outputs[seeds[k]].isSeed() = 1;
        }
      }
    };

  }  // namespace

  void HGCalCLUEsteringAlgoWrapper::run(Queue& queue,
                                        uint32_t nRecHits,
                                        float dc,
                                        float kappa,
                                        float outlierDeltaFactor,
                                        bool isScintillator,
                                        std::span<const uint32_t> batchItemSizes,
                                        HGCalSoARecHitsDeviceCollection::ConstView recHits,
                                        HGCalSoARecHitsExtraDeviceCollection::View recHitsClustering) const {
    // Nothing to do for an empty event: return 0 clusters
    if (nRecHits == 0) {
      auto nClusters = make_device_view<unsigned int>(queue, recHitsClustering.numberOfClustersScalar());
      alpaka::memset(queue, nClusters, 0x0);
      return;
    }

    const auto items = 256u;

    auto isSeedView = make_device_view(queue, recHitsClustering.isSeed().data(), nRecHits);
    alpaka::fill(queue, isSeedView, static_cast<uint8_t>(0));

    clue::ConstPointsDevice<2, float> d_points(queue,
                                               static_cast<int32_t>(nRecHits),
                                               recHits.dim1().data(),
                                               recHits.dim2().data(),
                                               recHits.energy().data(),
                                               recHitsClustering.clusterIndex().data());
    d_points.set_density_uncertainty(recHits.sigmaNoise());
    d_points.set_tags(recHits.detid());

    clue::Clusterer<2> algo(queue, dc, kappa, dc * outlierDeltaFactor);
    if (isScintillator) {
      // Scintillator (BH) cells cluster in (eta, phi)
      constexpr float kTwoPi = 2.f * static_cast<float>(M_PI);
      algo.setWrappedCoordinates(0, 1);
      algo.make_clusters(queue,
                         d_points,
                         batchItemSizes,
                         clue::PeriodicEuclideanMetric<2, float>{0.f, kTwoPi},
                         clue::FlatKernel<float>{0.5f});
    } else {
      algo.make_clusters(queue, d_points, batchItemSizes);
    }
    alpaka::wait(queue);

    auto h_nClusters = make_host_buffer<unsigned int>(queue);
    *h_nClusters.data() = static_cast<unsigned int>(d_points.n_clusters());
    auto d_nClusters = make_device_view<unsigned int>(queue, recHitsClustering.numberOfClustersScalar());
    alpaka::memcpy(queue, d_nClusters, h_nClusters);

    auto seeds = algo.getSeeds();
    const auto nseeds = static_cast<uint32_t>(seeds.size());
    if (nseeds > 0) {
      const auto seedGroups = divide_up_by(nseeds, items);
      const auto seedWorkDiv = make_workdiv<Acc1D>(seedGroups, items);
      alpaka::exec<Acc1D>(queue, seedWorkDiv, KernelSetSeeds{}, seeds.data(), recHitsClustering, nseeds);
    }
    alpaka::wait(queue);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
