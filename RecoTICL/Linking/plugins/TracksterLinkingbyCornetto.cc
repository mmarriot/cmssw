#include <algorithm>
#include <cmath>
#include <memory>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "TracksterLinkingbyCornetto.h"

using namespace ticl;

TracksterLinkingbyCornetto::TracksterLinkingbyCornetto(const edm::ParameterSet& conf,
                                                       edm::ConsumesCollector iC,
                                                       cms::Ort::ONNXRuntime const* onnxRuntime)
    : TracksterLinkingAlgoBase(conf, iC, onnxRuntime),
      pairTest_(cornetto::makePairTest(conf)),
      minEmittedEnergy_(conf.getParameter<float>("minEmittedEnergy")),
      minEmittedPt_(conf.getParameter<float>("minEmittedPt")),
      minEmittedEnergyEM_(conf.getParameter<float>("minEmittedEnergyEM")),
      minEmittedPtEM_(conf.getParameter<float>("minEmittedPtEM")),
      minEmittedEnergyAnyPt_(conf.getParameter<float>("minEmittedEnergyAnyPt")),
      mergeFragments_(conf.getParameter<bool>("mergeFragments")),
      mergeDropped_(conf.getParameter<bool>("mergeDropped")),
      mergeHadronicOnly_(conf.getParameter<bool>("mergeHadronicOnly")),
      mergeMinParentEnergy_(conf.getParameter<float>("mergeMinParentEnergy")),
      mergeCone_{conf.getParameter<float>("mergeMaxBackwardDistance"),
                 conf.getParameter<float>("mergeMaxLongitudinalDistance"),
                 conf.getParameter<float>("mergeTransverseRadius0"),
                 conf.getParameter<float>("mergeTransverseSlope")},
      mergeEMCone_{mergeCone_.maxBackward,
                   mergeCone_.maxLongitudinal,
                   conf.getParameter<float>("mergeEMTransverseRadius0"),
                   conf.getParameter<float>("mergeEMTransverseSlope")} {
  if (!cornetto::coneValid(mergeCone_, mergeEMCone_))
    throw cms::Exception("Configuration")
        << "TracksterLinkingbyCornetto: invalid merge cone parameters, see the parameter descriptions. merge cone "
        << mergeCone_.maxBackward << " " << mergeCone_.maxLongitudinal << " " << mergeCone_.radius0 << " "
        << mergeCone_.slope << ", merge EM cone " << mergeEMCone_.radius0 << " " << mergeEMCone_.slope;
}

bool TracksterLinkingbyCornetto::emits(float energy, float pt, float em) const {
  const bool isEM = pairTest_.typeVetoProbability > 0.f && em > pairTest_.typeVetoProbability;
  if (!(energy >= (isEM ? minEmittedEnergyEM_ : minEmittedEnergy_)))
    return false;
  return pt >= (isEM ? minEmittedPtEM_ : minEmittedPt_) ||
         (minEmittedEnergyAnyPt_ > 0.f && energy >= minEmittedEnergyAnyPt_);
}

// Each trackster searches the tile of its endcap within the window of cornetto::searchWindow, which
// contains every anchor that can pass the cone.
template <typename Allowed, typename Test>
std::vector<unsigned int> TracksterLinkingbyCornetto::findParentsIn(
    const Features& f, const Tiles& tiles, const cornetto::Cone& cone, Allowed&& allowed, Test&& test) const {
  const auto n = static_cast<unsigned int>(f.energy.size());
  const auto ranksAbove = [&f](unsigned int a, unsigned int b) { return cornetto::ranksAbove(f.energy, a, b); };
  std::vector<unsigned int> parents(n);
  for (unsigned int i = 0; i < n; ++i) {
    unsigned int best = i;
    if (f.valid[i] && f.rad[i] > 0.f) {
      const float absEta = std::abs(f.eta[i]);
      const auto window = cornetto::searchWindow(cone, absEta, f.rad[i], pairTest_.etaWindow);
      const float etaMin = std::clamp(absEta - window.dEta, TileConstants::minEta, TileConstants::maxEta);
      const float etaMax = std::clamp(absEta + window.dEta, TileConstants::minEta, TileConstants::maxEta);
      auto const& tile = tiles[f.bary[i].z() > 0.f ? 1 : 0];
      const auto box = tile.searchBoxEtaPhi(etaMin, etaMax, f.phi[i] - window.dPhi, f.phi[i] + window.dPhi);
      for (int etaBin = box[0]; etaBin <= box[1]; ++etaBin) {
        for (int phiBin = box[2]; phiBin <= box[3]; ++phiBin) {
          for (unsigned int j : tile[tile.globalBin(etaBin, phiBin % TileConstants::nPhiBins)]) {
            const bool candidate = j != i && ranksAbove(j, i) && (best == i || ranksAbove(j, best)) && allowed(j, i);
            if (candidate && test(j, i))
              best = j;
          }
        }
      }
    }
    parents[i] = best;
  }
  return parents;
}

std::vector<unsigned int> TracksterLinkingbyCornetto::findParents(const Features& f, const Tiles& tiles) const {
  return findParentsIn(
      f,
      tiles,
      pairTest_.cone,
      [](unsigned int, unsigned int) { return true; },
      [&](unsigned int a, unsigned int b) { return compatible(f, a, b); });
}

TracksterLinkingbyCornetto::Features TracksterLinkingbyCornetto::groupFeatures(
    const Features& f, const std::vector<std::vector<unsigned int>>& groups) const {
  constexpr float kInvSpeedOfLight = 1.f / 29.9792458f;  // ns/cm
  const std::size_t n = groups.size();
  Features g;
  g.bary.resize(n);
  g.energy.resize(n);
  g.eta.resize(n);
  g.phi.resize(n);
  g.rad.resize(n);
  g.em.resize(n);
  g.time.resize(n);
  g.timeError.resize(n);
  g.valid.resize(n);
  for (std::size_t k = 0; k < n; ++k) {
    float energy = 0.f, emWeight = 0.f, em = 0.f;
    math::XYZVectorF weighted;
    for (unsigned int i : groups[k]) {
      if (!f.valid[i])
        continue;
      energy += f.energy[i];
      weighted += f.bary[i] * f.energy[i];
      if (f.em[i] >= 0.f) {
        emWeight += f.energy[i];
        em += f.em[i] * f.energy[i];
      }
    }
    const auto bary = energy > 0.f ? weighted / energy : weighted;
    const float rad = std::sqrt(bary.Mag2());
    // Member times moved to the group barycenter along the direction from the origin, weighted by their errors.
    float sumW = 0.f, sumWT = 0.f;
    for (unsigned int i : groups[k]) {
      if (!f.valid[i] || !(f.timeError[i] > 0.f))
        continue;
      const float w = 1.f / (f.timeError[i] * f.timeError[i]);
      sumW += w;
      sumWT += w * (f.time[i] + (rad - f.rad[i]) * kInvSpeedOfLight);
    }
    g.bary[k] = bary;
    g.energy[k] = energy;
    g.eta[k] = bary.eta();
    g.phi[k] = bary.phi();
    g.rad[k] = rad;
    g.em[k] = emWeight > 0.f ? em / emWeight : -1.f;
    g.time[k] = sumW > 0.f ? sumWT / sumW : 0.f;
    g.timeError[k] = sumW > 0.f ? 1.f / std::sqrt(sumW) : -1.f;
    g.valid[k] = energy > 0.f && std::isfinite(bary.x()) && std::isfinite(bary.y()) && std::isfinite(bary.z()) &&
                 std::isfinite(energy);
  }
  return g;
}

bool TracksterLinkingbyCornetto::mergeCompatible(const Features& g, unsigned int a, unsigned int b) const {
  if (mergeHadronicOnly_ && pairTest_.typeVetoProbability > 0.f &&
      (g.em[a] > pairTest_.typeVetoProbability || g.em[b] > pairTest_.typeVetoProbability))
    return false;
  if (!(g.energy[a] >= mergeMinParentEnergy_))
    return false;
  return compatibleIn(mergeCone_, mergeEMCone_, g, a, b);
}

std::vector<std::vector<unsigned int>> TracksterLinkingbyCornetto::mergeGroups(
    const Features& f, std::vector<std::vector<unsigned int>> groups) const {
  if (!mergeFragments_ && !mergeDropped_)
    return groups;
  const Features g = groupFeatures(f, groups);
  const auto n = static_cast<unsigned int>(groups.size());
  std::vector<bool> emitted(n);
  for (unsigned int k = 0; k < n; ++k)
    emitted[k] = g.valid[k] && emits(g.energy[k], g.energy[k] / std::cosh(g.eta[k]), g.em[k]);
  const auto tiles = std::make_unique<Tiles>();
  for (unsigned int k = 0; k < n; ++k)
    if (g.valid[k])
      (*tiles)[g.bary[k].z() > 0.f ? 1 : 0].fill(g.eta[k], g.phi[k], k);
  const auto allowed = [&](unsigned int a, unsigned int b) {
    return emitted[a] && (emitted[b] ? mergeFragments_ : mergeDropped_);
  };
  const auto test = [&](unsigned int a, unsigned int b) { return mergeCompatible(g, a, b); };
  const auto merged = cornetto::groups(cornetto::roots(findParentsIn(g, *tiles, mergeCone_, allowed, test)),
                                       [&](unsigned int root, unsigned int k) { return test(root, k); });
  std::vector<std::vector<unsigned int>> out;
  out.reserve(merged.size());
  for (auto const& m : merged) {
    std::vector<unsigned int> members;
    for (unsigned int k : m)
      members.insert(members.end(), groups[k].begin(), groups[k].end());
    std::sort(members.begin(), members.end());
    out.push_back(std::move(members));
  }
  std::sort(out.begin(), out.end(), [](auto const& x, auto const& y) { return x.front() < y.front(); });
  return out;
}

std::vector<std::vector<unsigned int>> TracksterLinkingbyCornetto::firstPassGroups(const Features& f) const {
  // About 200 kB, so on the heap rather than on the stack.
  const auto tiles = std::make_unique<Tiles>();
  for (unsigned int i = 0; i < f.energy.size(); ++i)
    if (f.valid[i])
      (*tiles)[f.bary[i].z() > 0.f ? 1 : 0].fill(f.eta[i], f.phi[i], i);

  // A trackster that is not compatible with its root leaves it.
  return cornetto::groups(cornetto::roots(findParents(f, *tiles)),
                          [&](unsigned int root, unsigned int i) { return compatible(f, root, i); });
}

void TracksterLinkingbyCornetto::linkTracksters(
    const Inputs& input,
    std::vector<Trackster>& resultTracksters,
    std::vector<std::vector<unsigned int>>& linkedResultTracksters,
    std::vector<std::vector<unsigned int>>& linkedTracksterIdToInputTracksterId) {
  const Features f = cornetto::computeFeatures(input.tracksters);
  auto groups = mergeGroups(f, firstPassGroups(f));

  resultTracksters.reserve(resultTracksters.size() + groups.size());
  linkedResultTracksters.reserve(linkedResultTracksters.size() + groups.size());
  linkedTracksterIdToInputTracksterId.reserve(linkedTracksterIdToInputTracksterId.size() + groups.size());
  for (auto& group : groups) {
    Trackster merged;
    merged.mergeTracksters(input.tracksters, group);
    // The merged trackster has no barycenter yet, so the pt uses the energy-weighted member barycenter.
    math::XYZVectorF weighted;
    float emWeight = 0.f, em = 0.f;
    for (unsigned int i : group) {
      weighted += input.tracksters[i].barycenter() * input.tracksters[i].raw_energy();
      if (f.em[i] >= 0.f) {
        emWeight += f.energy[i];
        em += f.em[i] * f.energy[i];
      }
    }
    const float pt = merged.raw_energy() / std::cosh(weighted.eta());
    if (!emits(merged.raw_energy(), pt, emWeight > 0.f ? em / emWeight : -1.f))
      continue;
    linkedResultTracksters.push_back({static_cast<unsigned int>(resultTracksters.size())});
    resultTracksters.push_back(std::move(merged));
    linkedTracksterIdToInputTracksterId.push_back(std::move(group));
  }
}
