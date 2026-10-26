#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/HGCalReco/interface/Common.h"
#include "DataFormats/Math/interface/normalizedPhi.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/prefixScan.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "RecoTICL/Linking/interface/CornettoGrouping.h"

#include "CornettoGroupsAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;
  using ticl::cornetto::Feature;
  using ticl::cornetto::PairTest;

  namespace {

    // The tiles of TracksterLinkingbyCornetto: one TICLLayerTileT<TileConstants> per endcap, index 1 for z > 0, here
    // with the bins of the two endcaps one after the other. The bin functions are the ones of TICLLayerTileT.
    using Tile = ticl::TileConstants;
    constexpr uint32_t kBins = 2 * Tile::nBins;

    ALPAKA_FN_ACC int etaBinOf(float eta) {
      constexpr float etaRange = Tile::maxEta - Tile::minEta;
      float r = Tile::nEtaBins / etaRange;
      int bin = (std::abs(eta) - Tile::minEta) * r;
      return std::clamp(bin, 0, Tile::nEtaBins - 1);
    }

    ALPAKA_FN_ACC int phiBinOf(float phi) {
      auto normPhi = normalizedPhi(phi);
      float r = Tile::nPhiBins * M_1_PI * 0.5f;
      int bin = (normPhi + M_PI) * r;
      return bin;
    }

    // The bin of a trackster. A phi of exactly pi gives the phi bin nPhiBins, which TICLLayerTileT::fill rejects with
    // an exception; here it is the first phi bin, which is how the search reads that bin.
    ALPAKA_FN_ACC uint32_t tileBin(const Feature& f) {
      return (f.z > 0.f ? Tile::nBins : 0) + phiBinOf(f.phi) % Tile::nPhiBins + etaBinOf(f.eta) * Tile::nPhiBins;
    }

    // TICLLayerTileT::searchBoxEtaPhi for |eta| limits, which are positive: the eta bins and the phi bins, with the
    // upper phi bin above nPhiBins when the window crosses the phi boundary.
    ALPAKA_FN_ACC std::array<int, 4> searchBox(float etaMin, float etaMax, float phiMin, float phiMax) {
      if (etaMax - etaMin < 0)
        return {{0, 0, 0, 0}};
      int phiBinMin = phiBinOf(phiMin);
      int phiBinMax = phiBinOf(phiMax);
      if (phiBinMax < phiBinMin)
        phiBinMax += Tile::nPhiBins;
      return {{etaBinOf(etaMin), etaBinOf(etaMax), phiBinMin, phiBinMax}};
    }

    // The number of tracksters in each bin: as TICLLayerTileT, the tiles hold the valid tracksters only.
    struct KernelCountTiles {
      ALPAKA_FN_ACC void operator()(Acc1D const& acc, const Feature* features, uint32_t n, uint32_t* counts) const {
        for (auto i : uniform_elements(acc, n))
          if (features[i].valid)
            alpaka::atomicAdd(acc, &counts[tileBin(features[i])], 1u, alpaka::hierarchy::Blocks{});
      }
    };

    // The tracksters of each bin, from offsets[bin] to offsets[bin + 1], in no particular order.
    struct KernelFillTiles {
      ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                    const Feature* features,
                                    uint32_t n,
                                    const uint32_t* offsets,
                                    uint32_t* cursor,
                                    uint32_t* content) const {
        for (auto i : uniform_elements(acc, n)) {
          if (features[i].valid) {
            const uint32_t bin = tileBin(features[i]);
            content[offsets[bin] + alpaka::atomicAdd(acc, &cursor[bin], 1u, alpaka::hierarchy::Blocks{})] = i;
          }
        }
      }
    };

    // Step 1, TracksterLinkingbyCornetto::findParentsIn with the cone of the pair test: the parent of trackster i is
    // the highest-ranked trackster above it in the search window of its endcap that is compatible with it as anchor,
    // else i itself. The result does not depend on the order of the tracksters in a bin.
    struct KernelParents {
      ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                    PairTest pairTest,
                                    const Feature* features,
                                    uint32_t n,
                                    const uint32_t* offsets,
                                    const uint32_t* content,
                                    uint32_t* parents) const {
        for (auto i : uniform_elements(acc, n)) {
          const Feature& fi = features[i];
          uint32_t best = i;
          if (fi.valid && fi.rad > 0.f) {
            const float absEta = std::abs(fi.eta);
            const auto window = ticl::cornetto::searchWindow(pairTest.cone, absEta, fi.rad, pairTest.etaWindow);
            // local copies of the tile limits: std::clamp takes references, which the device cannot take to them
            const float minEta = Tile::minEta;
            const float maxEta = Tile::maxEta;
            const float etaMin = std::clamp(absEta - window.dEta, minEta, maxEta);
            const float etaMax = std::clamp(absEta + window.dEta, minEta, maxEta);
            const uint32_t endcap = fi.z > 0.f ? Tile::nBins : 0;
            const auto box = searchBox(etaMin, etaMax, fi.phi - window.dPhi, fi.phi + window.dPhi);
            for (int etaBin = box[0]; etaBin <= box[1]; ++etaBin) {
              for (int phiBin = box[2]; phiBin <= box[3]; ++phiBin) {
                const uint32_t bin = endcap + phiBin % Tile::nPhiBins + etaBin * Tile::nPhiBins;
                for (uint32_t k = offsets[bin]; k < offsets[bin + 1]; ++k) {
                  const uint32_t j = content[k];
                  const bool candidate =
                      j != i && ticl::cornetto::ranksAbove(features[j].energy, j, fi.energy, i) &&
                      (best == i || ticl::cornetto::ranksAbove(features[j].energy, j, features[best].energy, best));
                  if (candidate && pairTest.compatible(features[j], fi))
                    best = j;
                }
              }
            }
          }
          parents[i] = best;
        }
      }
    };

    // Steps 2 and 3, cornetto::roots and cornetto::groups: the root of trackster i along its parents, and its label,
    // the root when the pair (root, i) is compatible, else i itself.
    struct KernelLabels {
      ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                    PairTest pairTest,
                                    const Feature* features,
                                    uint32_t n,
                                    const uint32_t* parents,
                                    uint32_t* labels) const {
        for (auto i : uniform_elements(acc, n)) {
          uint32_t root = i;
          while (parents[root] != root)
            root = parents[root];
          labels[i] = (root == i || pairTest.compatible(features[root], features[i])) ? root : i;
        }
      }
    };

  }  // namespace

  void CornettoGroupsAlgo::run(
      Queue& queue, const PairTest& pairTest, const Feature* features, uint32_t n, uint32_t* labels) const {
    auto counts = make_device_buffer<uint32_t[]>(queue, kBins);
    auto cursor = make_device_buffer<uint32_t[]>(queue, kBins);
    auto offsets = make_device_buffer<uint32_t[]>(queue, kBins + 1);
    auto content = make_device_buffer<uint32_t[]>(queue, n);
    auto parents = make_device_buffer<uint32_t[]>(queue, n);
    alpaka::memset(queue, counts, 0x0);
    alpaka::memset(queue, cursor, 0x0);
    alpaka::memset(queue, offsets, 0x0);

    constexpr uint32_t blockSize = 256;
    const auto workDiv = make_workdiv<Acc1D>(divide_up_by(n, blockSize), blockSize);
    alpaka::exec<Acc1D>(queue, workDiv, KernelCountTiles{}, features, n, counts.data());
    // inclusive scan of the counts, so offsets[bin + 1] is where the tracksters of bin end
    iterativePrefixScan<Acc1D>(counts.data(), offsets.data() + 1, kBins, queue);
    alpaka::exec<Acc1D>(queue, workDiv, KernelFillTiles{}, features, n, offsets.data(), cursor.data(), content.data());
    alpaka::exec<Acc1D>(
        queue, workDiv, KernelParents{}, pairTest, features, n, offsets.data(), content.data(), parents.data());
    alpaka::exec<Acc1D>(queue, workDiv, KernelLabels{}, pairTest, features, n, parents.data(), labels);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
