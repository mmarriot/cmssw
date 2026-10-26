#ifndef RecoTICL_Linking_CornettoFeatures_h
#define RecoTICL_Linking_CornettoFeatures_h

// The per-trackster quantities that the Cornetto linking reads, from the input tracksters: TracksterLinkingbyCornetto
// on the host and CornettoGroupsProducer for the device compute them with computeFeatures.

#include <vector>

#include "DataFormats/Common/interface/MultiSpan.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "RecoTICL/Linking/interface/CornettoPairTest.h"

namespace ticl::cornetto {

  struct Features {
    std::vector<math::XYZVectorF> bary;
    std::vector<float> energy;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<float> rad;   // |bary|
    std::vector<float> em;    // EM probability, negative without a PID
    std::vector<float> time;  // local to the barycenter [ns]
    std::vector<float> timeError;
    std::vector<bool> valid;  // finite barycenter and energy
  };

  Features computeFeatures(const edm::MultiSpan<Trackster>& tracksters);

  // The quantities of trackster i.
  inline Feature feature(const Features& f, unsigned int i) {
    return {f.bary[i].x(),
            f.bary[i].y(),
            f.bary[i].z(),
            f.rad[i],
            f.eta[i],
            f.phi[i],
            f.energy[i],
            f.em[i],
            f.time[i],
            f.timeError[i],
            f.valid[i]};
  }

}  // namespace ticl::cornetto

#endif
