#ifndef RecoTICL_Linking_CornettoPairTest_h
#define RecoTICL_Linking_CornettoPairTest_h

// The Cornetto pair test (see RecoTICL/Linking/plugins/TracksterLinkingbyCornetto.h), on the host and on the device:
// TracksterLinkingbyCornetto and CornettoGroupsProducer run these functions.

#include <cmath>

#include <alpaka/core/Common.hpp>

#include "DataFormats/Math/interface/deltaPhi.h"
#include "RecoTICL/Linking/interface/CornettoGeometry.h"

namespace edm {
  class ParameterSet;
  class ParameterSetDescription;
}  // namespace edm

namespace ticl::cornetto {

  // The quantities of one trackster that the pair test and the parent search read.
  struct Feature {
    float x, y, z;  // barycenter
    float rad;      // |barycenter|
    float eta, phi;
    float energy;
    float em;    // EM probability, negative without a PID
    float time;  // local to the barycenter [ns]
    float timeError;
    bool valid;  // finite barycenter and energy
  };

  struct PairTest {
    float etaWindow;
    Cone cone;
    Cone emCone;  // cone of an EM anchor, inside cone
    bool useTimeGate;
    float timeNSigma2;
    float timeFloor2;
    float typeVetoProbability;

    // A NaN time error or EM probability disables the time gate or the type veto for that pair.
    ALPAKA_FN_HOST_ACC bool timeCompatible(
        float tAnchor, float errAnchor, float tOther, float errOther, float path) const {
      if (!useTimeGate || !(errAnchor > 0.f) || !(errOther > 0.f))
        return true;
      constexpr float kInvSpeedOfLight = 1.f / 29.9792458f;  // ns/cm
      const float dt = tOther - tAnchor - path * kInvSpeedOfLight;
      return dt * dt <= timeNSigma2 * (errAnchor * errAnchor + errOther * errOther + timeFloor2);
    }

    // Whether b can join anchor a, with the cones of the test (the merge cones in the second pass).
    ALPAKA_FN_HOST_ACC bool compatibleIn(const Cone& testCone,
                                         const Cone& testEMCone,
                                         const Feature& a,
                                         const Feature& b) const {
      if (!(a.z * b.z > 0.f))
        return false;  // same endcap only
      if (!(std::abs(b.eta - a.eta) <= etaWindow) || !(std::abs(reco::deltaPhi(b.phi, a.phi)) <= etaWindow))
        return false;
      const bool emAnchor = typeVetoProbability > 0.f && a.em > typeVetoProbability;
      float s;
      if (!separationInCone(emAnchor ? testEMCone : testCone, a.x, a.y, a.z, a.rad, b.x, b.y, b.z, s))
        return false;
      if (!timeCompatible(a.time, a.timeError, b.time, b.timeError, s))
        return false;
      if (typeVetoProbability > 0.f && a.em >= 0.f && b.em >= 0.f) {
        const float hi = typeVetoProbability;
        const float lo = 1.f - typeVetoProbability;
        if ((a.em > hi && b.em < lo) || (a.em < lo && b.em > hi))
          return false;
      }
      return true;
    }

    ALPAKA_FN_HOST_ACC bool compatible(const Feature& a, const Feature& b) const {
      return compatibleIn(cone, emCone, a, b);
    }
  };

  // Whether a cone and the cone of its EM anchors are valid: the EM cone inside the cone.
  inline bool coneValid(const Cone& c, const Cone& em) {
    return c.maxLongitudinal > 0.f && c.maxBackward >= 0.f && c.radius0 >= 0.f && c.slope >= 0.f && em.radius0 >= 0.f &&
           em.radius0 <= c.radius0 && em.slope >= 0.f && em.slope <= c.slope;
  }

  // The pair test from its parameters (fillPairTestDescription); throws on invalid parameters.
  PairTest makePairTest(const edm::ParameterSet& conf);
  void fillPairTestDescription(edm::ParameterSetDescription& desc);

}  // namespace ticl::cornetto

#endif
