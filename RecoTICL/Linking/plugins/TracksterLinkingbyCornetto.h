#ifndef RecoTICL_Linking_TracksterLinkingbyCornetto_h
#define RecoTICL_Linking_TracksterLinkingbyCornetto_h

// Trackster linking by a parent tree.
//
// Pair test with anchor a and other trackster b (RecoTICL/Linking/interface/CornettoPairTest.h, which the device
// version runs too, see TracksterLinkingbyCornettoDevice.h): the
// axis points from the origin to the barycenter of a. With D = bary_b - bary_a, s = D . axis and
// dT = |D - s axis|, the pair is compatible when -maxBackwardDistance <= s <= maxLongitudinalDistance,
// dT <= transverseRadius0 + transverseSlope * |s|, |deta| and |dphi| <= etaWindow, both are in the same
// endcap, and:
// - EM anchor: when typeVetoProbability > 0 and the EM probability of a is above typeVetoProbability, the
//   transverse limit is emTransverseRadius0 + emTransverseSlope * |s| instead.
// - time: when timeCompatibilityNSigma > 0 and both time errors are positive, t_b - t_a - s / c agrees
//   with 0 within timeCompatibilityNSigma times the combined error, with timeResolutionFloor added in
//   quadrature. The trackster times must be local to the barycenter.
// - type: when typeVetoProbability > 0 and both carry a PID, the pair is refused if one EM probability is
//   above typeVetoProbability and the other below 1 - typeVetoProbability.
//
// Grouping, each step independent per trackster:
// 1. The parent of b is the highest-energy trackster a with more energy than b (equal energy: smaller
//    index) for which the pair (a, b) is compatible. A trackster without a parent is a root. A trackster
//    with a non-finite barycenter or energy is a root and is never a parent.
// 2. Each trackster follows its parents to its root.
// 3. A trackster stays with its root only if the pair (root, trackster) is compatible; else it is a group
//    of its own.
//
// Optional second pass on the groups (mergeFragments, mergeDropped): the same parent tree on the groups, with the
// pair test of the merge cone. A group has the raw energy, the energy-weighted barycenter, the energy-weighted EM
// probability of its members with a PID, and the error-weighted time of its members moved to its barycenter. Only a
// group that passes the emission rule is a parent. A group that passes it can have a parent when mergeFragments is
// true, a group that does not pass it when mergeDropped is true. With mergeHadronicOnly and a positive
// typeVetoProbability, a pair with an EM group (EM probability above typeVetoProbability) is refused.
//
// Emission rule: a group is emitted when its raw energy is at least minEmittedEnergy and its raw pt is at least
// minEmittedPt, or its raw energy is at least minEmittedEnergyAnyPt (when positive). A group with an EM probability
// above typeVetoProbability uses minEmittedEnergyEM and minEmittedPtEM instead. The pt uses the direction of the
// energy-weighted barycenter of the members.

#include <array>
#include <vector>

#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "RecoTICL/Linking/interface/CornettoFeatures.h"
#include "RecoTICL/Linking/interface/CornettoGeometry.h"
#include "RecoTICL/Linking/interface/CornettoGrouping.h"
#include "RecoTICL/Linking/interface/CornettoPairTest.h"
#include "RecoTICL/Linking/interface/TracksterLinkingAlgoBase.h"

namespace ticl {

  class TracksterLinkingbyCornetto : public TracksterLinkingAlgoBase {
  public:
    TracksterLinkingbyCornetto(const edm::ParameterSet& conf,
                               edm::ConsumesCollector iC,
                               cms::Ort::ONNXRuntime const* onnxRuntime = nullptr);
    ~TracksterLinkingbyCornetto() override = default;

    void linkTracksters(const Inputs& input,
                        std::vector<Trackster>& resultTracksters,
                        std::vector<std::vector<unsigned int>>& linkedResultTracksters,
                        std::vector<std::vector<unsigned int>>& linkedTracksterIdToInputTracksterId) override;

    void initialize(const HGCalDDDConstants* hgcons,
                    const ticlgeom::Tools rhtools,
                    const edm::ESHandle<MagneticField> bfieldH,
                    const edm::ESHandle<Propagator> propH) override {}

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc) {
      cornetto::fillPairTestDescription(iDesc);
      iDesc.add<float>("minEmittedEnergy", 1.f)
          ->setComment("Min raw energy [GeV] of an emitted group. Applies after linking.");
      iDesc.add<float>("minEmittedPt", 0.5f)
          ->setComment("Min raw pt [GeV] of an emitted group. Applies after linking.");
      iDesc.add<float>("minEmittedEnergyEM", 1.f)
          ->setComment("Min raw energy [GeV] of an emitted EM group (EM probability above typeVetoProbability).");
      iDesc.add<float>("minEmittedPtEM", 0.5f)->setComment("Min raw pt [GeV] of an emitted EM group.");
      iDesc.add<float>("minEmittedEnergyAnyPt", 0.f)
          ->setComment(
              "A group at or above this raw energy [GeV] is emitted at any pt. Values <= 0 disable it. The minimum "
              "energy still applies.");
      iDesc.add<bool>("mergeFragments", false)
          ->setComment("Second pass: an emitted group can join an emitted group of higher energy.");
      iDesc.add<bool>("mergeDropped", false)
          ->setComment("Second pass: a group that is not emitted can join an emitted group.");
      iDesc.add<bool>("mergeHadronicOnly", true)->setComment("Second pass: refuse a pair with an EM group.");
      iDesc.add<float>("mergeMinParentEnergy", 0.f)->setComment("Second pass: min raw energy [GeV] of a parent.");
      iDesc.add<float>("mergeMaxLongitudinalDistance", 60.f)
          ->setComment("Second pass: max separation [cm] downstream of the parent. Must be positive.");
      iDesc.add<float>("mergeMaxBackwardDistance", 20.f)
          ->setComment("Second pass: max separation [cm] upstream of the parent. Must not be negative.");
      iDesc.add<float>("mergeTransverseRadius0", 4.f)
          ->setComment("Second pass: max transverse distance at zero separation [cm]. Must not be negative.");
      iDesc.add<float>("mergeTransverseSlope", 0.05f)
          ->setComment(
              "Second pass: growth of the max transverse distance per cm of separation. Must not be negative.");
      iDesc.add<float>("mergeEMTransverseRadius0", 1.5f)
          ->setComment("Second pass: max transverse distance at zero separation [cm] for an EM parent.");
      iDesc.add<float>("mergeEMTransverseSlope", 0.05f)
          ->setComment("Second pass: growth of the max transverse distance per cm of separation for an EM parent.");
      TracksterLinkingAlgoBase::fillPSetDescription(iDesc);
    }

  protected:
    using Features = cornetto::Features;

    // The groups of the first pass: steps 1 to 3 of the grouping.
    virtual std::vector<std::vector<unsigned int>> firstPassGroups(const Features& f) const;

  private:
    using Tiles = std::array<TICLLayerTile, 2>;  // one per endcap, index 1 for z > 0

    // Features of the groups, from the features f of their members.
    Features groupFeatures(const Features& f, const std::vector<std::vector<unsigned int>>& groups) const;
    bool compatible(const Features& f, unsigned int a, unsigned int b) const {
      return pairTest_.compatible(cornetto::feature(f, a), cornetto::feature(f, b));
    }
    bool compatibleIn(const cornetto::Cone& cone,
                      const cornetto::Cone& emCone,
                      const Features& f,
                      unsigned int a,
                      unsigned int b) const {
      return pairTest_.compatibleIn(cone, emCone, cornetto::feature(f, a), cornetto::feature(f, b));
    }
    bool mergeCompatible(const Features& g, unsigned int a, unsigned int b) const;
    bool emits(float energy, float pt, float em) const;
    // The parent of each node: the highest-ranked node a above it, with allowed(a, node) true, that passes
    // test(a, node).
    template <typename Allowed, typename Test>
    std::vector<unsigned int> findParentsIn(
        const Features& f, const Tiles& tiles, const cornetto::Cone& cone, Allowed&& allowed, Test&& test) const;
    std::vector<unsigned int> findParents(const Features& f, const Tiles& tiles) const;
    // The groups after the second pass.
    std::vector<std::vector<unsigned int>> mergeGroups(const Features& f,
                                                       std::vector<std::vector<unsigned int>> groups) const;

    const cornetto::PairTest pairTest_;
    const float minEmittedEnergy_;
    const float minEmittedPt_;
    const float minEmittedEnergyEM_;
    const float minEmittedPtEM_;
    const float minEmittedEnergyAnyPt_;
    const bool mergeFragments_;
    const bool mergeDropped_;
    const bool mergeHadronicOnly_;
    const float mergeMinParentEnergy_;
    const cornetto::Cone mergeCone_;
    const cornetto::Cone mergeEMCone_;  // merge cone of an EM parent, inside mergeCone_
  };

}  // namespace ticl

#endif
