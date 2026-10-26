#include <cmath>

#include "RecoTICL/Linking/interface/CornettoFeatures.h"

ticl::cornetto::Features ticl::cornetto::computeFeatures(const edm::MultiSpan<Trackster>& tracksters) {
  const std::size_t n = tracksters.size();
  Features f;
  f.bary.resize(n);
  f.energy.resize(n);
  f.eta.resize(n);
  f.phi.resize(n);
  f.rad.resize(n);
  f.em.resize(n);
  f.time.resize(n);
  f.timeError.resize(n);
  f.valid.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    auto const& ts = tracksters[i];
    const auto& b = ts.barycenter();
    f.bary[i] = b;
    f.energy[i] = ts.raw_energy();
    f.valid[i] = std::isfinite(b.x()) && std::isfinite(b.y()) && std::isfinite(b.z()) && std::isfinite(f.energy[i]);
    f.eta[i] = b.eta();
    f.phi[i] = b.phi();
    f.rad[i] = std::sqrt(b.Mag2());
    f.time[i] = ts.time();
    f.timeError[i] = ts.timeError();
    // -1 for a trackster with no PID: no probability, or only unknown (a merged trackster).
    float sum = 0.f;
    for (float v : ts.id_probabilities())
      sum += v;
    f.em[i] = sum > ts.id_probability(Trackster::ParticleType::unknown)
                  ? (ts.id_probability(Trackster::ParticleType::photon) +
                     ts.id_probability(Trackster::ParticleType::electron)) /
                        sum
                  : -1.f;
  }
  return f;
}
