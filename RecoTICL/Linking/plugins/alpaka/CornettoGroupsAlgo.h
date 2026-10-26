#ifndef RecoTICL_Linking_plugins_alpaka_CornettoGroupsAlgo_h
#define RecoTICL_Linking_plugins_alpaka_CornettoGroupsAlgo_h

#include <cstdint>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoTICL/Linking/interface/CornettoPairTest.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // The first pass of the Cornetto grouping (TracksterLinkingbyCornetto::firstPassGroups) on the device, one thread per
  // trackster: the parent of each trackster in the tiles of the endcaps, its root, and its group label, the root when
  // the pair (root, trackster) is compatible, else the trackster itself.
  class CornettoGroupsAlgo {
  public:
    void run(Queue& queue,
             const ticl::cornetto::PairTest& pairTest,
             const ticl::cornetto::Feature* features,
             uint32_t n,
             uint32_t* labels) const;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoTICL_Linking_plugins_alpaka_CornettoGroupsAlgo_h
