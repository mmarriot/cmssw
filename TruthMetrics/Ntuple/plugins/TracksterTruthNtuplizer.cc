// Writes, per event, everything the trackster-building metrics need, so that the clustering can be
// re-run and scored outside CMSSW (the ticltune python package):
//
//   layer clusters   energy, position, layer, side, detector, seed DetId (CLUEstering tie-break tag),
//                    the iteration mask, and the energy that no in-time truth particle accounts for
//   truth table      sparse (layer cluster, base particle, energy): the base particle's share of the
//                    layer cluster, rechit-energy weighted (see below)
//   rechits          EVERY HGCAL rechit (energy, position, layer, detector, the layer cluster holding most of
//                    it or -1, its no-truth energy) and the layer-cluster membership with fractions (lch_*)
//   atoms            every truth-graph particle with energy in HGCAL rechits, inside the calorimeter
//                    included: its parent atom, unit and interaction, and the sparse (rechit, atom, energy)
//                    table  s(c, a) = E_rechit(c) * E_sim(a, c) / E_sim(all in-time, c)        (tra_*)
//   units            what the atoms are grouped into (tuning/V3_TRUTH_AND_METRICS.md, 4.2): the atom's nearest
//                    reconstructableFinalState ancestor-or-self, except below a pi0, where it is the pi0's decay
//                    daughter; with no such ancestor, the root (un_kind 0 / 1 / 2)
//   base particles   the caloBoundary particles that own calorimeter hits, with their interaction
//                    (signal or which in-time pileup collision) and their origin
//   CMSSW tracksters the CLUEstering assignment per layer cluster and any trackster collections, for
//                    the closure test of the emulator; with tracksterDetails also the regressed energy,
//                    barycenter, time, PID probabilities and the first linked track of every trackster
//   TICL candidates  (optional) pdg, charge, four-momentum, raw energy, time, PID probabilities, the
//                    tracksters (keys into the candidate producer's merged trackster collection) and the
//                    track (index into the tracks collection) of every TICLCandidate
//   tracks           (optional) the tracks above |eta| and pT thresholds: kinematics, charge, quality
//   gen jets         (optional) pT, eta, phi, E of every gen jet above a pT threshold
//   final state      the signal particles of the graph's reconstructableFinalState level (pdg, momentum),
//                    and for every base particle the row of its origin in that table (bp_originRow)
//
// Base particle of a truth particle p (baseLevel = "crossing", the default): its nearest ancestor-or-self that
// crossed into the calorimeter (a boundary checkpoint in the graph, not back-scattered). Every crossing counts,
// also below another one: an electron that crossed after radiating keeps its own deposits, and each brem photon or
// conversion leg that crossed on its own keeps its own. The graph's caloBoundary *level* is an antichain (only the
// outermost crossing of nested ones), which would give the electron everything; baseLevel = "caloBoundary" uses it.
// A particle with no such ancestor (rare: hits from below the boundary that no crossing particle owns) falls back to
// its nearest reconstructableFinalState ancestor-or-self, then to its root, and is flagged by `bp_kind`.
//
// Energy of base particle b in layer cluster l:
//   s(l, b) = sum over the cells c of l:  fraction_l(c) * E_rechit(c) * E_sim(b, c) / E_sim(all in-time, c)
// A cell with a rechit but no in-time sim energy (out-of-time pileup, noise) goes to the layer
// cluster's no-truth energy. Out-of-time energy in a cell that also has in-time sim energy cannot be
// told apart and is shared among the in-time particles, as in the TICL associators.

#include <algorithm>
#include <cstdlib>
#include <tuple>
#include <string>
#include <unordered_map>
#include <vector>

#include "TTree.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/transform.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/TruthInfo/interface/LogicalGraphHitIndex.h"

namespace {
  constexpr int kNone = -1;
  enum BaseKind : int { kCaloBoundary = 0, kReconstructableFinalState = 1, kRoot = 2 };
  enum UnitKind : int { kUnitFinalState = 0, kUnitPi0Daughter = 1, kUnitRoot = 2 };
  constexpr int kPi0 = 111;

  // First parent through the first production vertex that has an incoming particle; -1 for a root.
  int firstParent(truth::Graph const& graph, uint32_t p) {
    for (auto v : graph.productionVertices(p)) {
      auto in = graph.incomingParticles(v);
      if (!in.empty())
        return static_cast<int>(in[0]);
    }
    return kNone;
  }
}  // namespace

class TracksterTruthNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TracksterTruthNtuplizer(edm::ParameterSet const&);
  void beginJob() override;
  void analyze(edm::Event const&, edm::EventSetup const&) override;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  struct TracksterBranches {
    std::vector<float> energy;
    std::vector<int> nLC;
    std::vector<int> lc;
    std::vector<float> mult;
    // tracksterDetails
    std::vector<float> regE, bx, by, bz, time, timeErr, prob;
    std::vector<int> trk;
    void clear() {
      for (auto* v : {&energy, &mult, &regE, &bx, &by, &bz, &time, &timeErr, &prob})
        v->clear();
      for (auto* v : {&nLC, &lc, &trk})
        v->clear();
    }
  };

  const edm::EDGetTokenT<truth::Graph> graphToken_;
  const edm::EDGetTokenT<truth::LogicalGraphHitIndex> hitIndexToken_;
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> layerClustersToken_;
  const edm::EDGetTokenT<std::vector<float>> maskToken_;
  const edm::InputTag assignmentTag_;
  edm::EDGetTokenT<std::vector<int32_t>> assignmentToken_;
  const std::vector<edm::InputTag> trackstersTags_;
  const std::vector<edm::EDGetTokenT<std::vector<ticl::Trackster>>> trackstersTokens_;
  const std::vector<edm::EDGetTokenT<HGCRecHitCollection>> recHitTokens_;
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geomToken_;
  const bool crossing_;  // baseLevel == "crossing": every boundary crossing, not only the caloBoundary level
  const bool tracksterDetails_;
  const edm::InputTag candidatesTag_, tracksTag_, genJetsTag_;
  edm::EDGetTokenT<std::vector<TICLCandidate>> candidatesToken_;
  edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet>> genJetsToken_;
  const double trackMinAbsEta_, trackMinPt_, genJetMinPt_;
  hgcal::RecHitTools rhtools_;

  TTree* tree_ = nullptr;
  unsigned run_ = 0, lumi_ = 0;
  unsigned long long event_ = 0;
  // layer clusters
  std::vector<float> lcE_, lcX_, lcY_, lcZ_, lcMask_, lcNoTruthE_, lcRecE_;
  std::vector<int> lcLayer_, lcSide_, lcDet_, lcNHits_, lcAlgo_, lcMissingRecHits_;
  std::vector<unsigned> lcSeed_;
  // truth table
  std::vector<int> trLC_, trBP_;
  std::vector<float> trE_;
  // every HGCAL rechit, the layer-cluster membership, atoms with their (rechit, atom) table, units
  std::vector<unsigned> rhId_;
  std::vector<float> rhE_, rhX_, rhY_, rhZ_, rhNoTruthE_, lchFrac_, traE_, atDepE_, unE_, unEta_, unPhi_, unDepE_;
  std::vector<int> rhLayer_, rhDet_, rhLC_, lchLC_, lchRH_, traRH_, traAT_;
  std::vector<int> atId_, atPdg_, atParent_, atUnit_, atSignal_, atBX_, atEvt_;
  std::vector<int> unId_, unPdg_, unKind_, unOrigin_, unFsRow_, unSignal_, unBX_, unEvt_;
  // base particles
  std::vector<int> bpId_, bpPdg_, bpSignal_, bpBX_, bpEvt_, bpKind_, bpOrigin_, bpOriginPdg_, bpOriginRow_;
  std::vector<float> bpE_, bpEta_, bpPhi_, bpSimE_, bpOriginE_, bpX_, bpY_, bpZ_;
  // signal final-state particles (reconstructableFinalState level)
  std::vector<int> fsId_, fsPdg_;
  std::vector<float> fsE_, fsPt_, fsEta_, fsPhi_;
  // CMSSW reference
  std::vector<int> clueAssignment_;
  std::vector<TracksterBranches> tracksters_;
  // TICL candidates, tracks, gen jets
  std::vector<int> candPdg_, candCharge_, candNTs_, candTs_, candTrk_;
  std::vector<float> candE_, candPt_, candEta_, candPhi_, candRawE_, candTime_, candTimeErr_, candProb_;
  std::vector<int> trkIdx_, trkCharge_, trkNHits_, trkHighPurity_;
  std::vector<float> trkPt_, trkEta_, trkPhi_, trkChi2n_, trkPtErr_;
  std::vector<float> gjPt_, gjEta_, gjPhi_, gjE_;
};

TracksterTruthNtuplizer::TracksterTruthNtuplizer(edm::ParameterSet const& cfg)
    : graphToken_(consumes(cfg.getParameter<edm::InputTag>("graph"))),
      hitIndexToken_(consumes(cfg.getParameter<edm::InputTag>("hitIndex"))),
      layerClustersToken_(consumes(cfg.getParameter<edm::InputTag>("layerClusters"))),
      maskToken_(consumes(cfg.getParameter<edm::InputTag>("layerClusterMask"))),
      assignmentTag_(cfg.getParameter<edm::InputTag>("clueAssignment")),
      trackstersTags_(cfg.getParameter<std::vector<edm::InputTag>>("tracksters")),
      trackstersTokens_(edm::vector_transform(
          trackstersTags_, [this](edm::InputTag const& t) { return consumes<std::vector<ticl::Trackster>>(t); })),
      recHitTokens_(edm::vector_transform(cfg.getParameter<std::vector<edm::InputTag>>("recHits"),
                                          [this](edm::InputTag const& t) { return consumes<HGCRecHitCollection>(t); })),
      geomToken_(esConsumes()),
      crossing_(cfg.getParameter<std::string>("baseLevel") == "crossing"),
      tracksterDetails_(cfg.getParameter<bool>("tracksterDetails")),
      candidatesTag_(cfg.getParameter<edm::InputTag>("candidates")),
      tracksTag_(cfg.getParameter<edm::InputTag>("tracks")),
      genJetsTag_(cfg.getParameter<edm::InputTag>("genJets")),
      trackMinAbsEta_(cfg.getParameter<double>("trackMinAbsEta")),
      trackMinPt_(cfg.getParameter<double>("trackMinPt")),
      genJetMinPt_(cfg.getParameter<double>("genJetMinPt")) {
  if (!crossing_ && cfg.getParameter<std::string>("baseLevel") != "caloBoundary")
    throw cms::Exception("Configuration") << "baseLevel must be crossing or caloBoundary";
  usesResource(TFileService::kSharedResource);
  if (!assignmentTag_.label().empty())
    assignmentToken_ = consumes(assignmentTag_);
  if (!candidatesTag_.label().empty())
    candidatesToken_ = consumes(candidatesTag_);
  if (!tracksTag_.label().empty())
    tracksToken_ = consumes(tracksTag_);
  if (!genJetsTag_.label().empty())
    genJetsToken_ = consumes(genJetsTag_);
  tracksters_.resize(trackstersTags_.size());
}

void TracksterTruthNtuplizer::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("events", "trackster-building truth ntuple");
  tree_->Branch("run", &run_);
  tree_->Branch("lumi", &lumi_);
  tree_->Branch("event", &event_);
  tree_->Branch("lc_E", &lcE_);
  tree_->Branch("lc_x", &lcX_);
  tree_->Branch("lc_y", &lcY_);
  tree_->Branch("lc_z", &lcZ_);
  tree_->Branch("lc_layer", &lcLayer_);
  tree_->Branch("lc_side", &lcSide_);
  tree_->Branch("lc_det", &lcDet_);
  tree_->Branch("lc_nhits", &lcNHits_);
  tree_->Branch("lc_algo", &lcAlgo_);
  tree_->Branch("lc_seed", &lcSeed_);
  tree_->Branch("lc_mask", &lcMask_);
  tree_->Branch("lc_noTruthE", &lcNoTruthE_);
  tree_->Branch("lc_recE", &lcRecE_);
  tree_->Branch("lc_missingRecHits", &lcMissingRecHits_);
  tree_->Branch("tr_lc", &trLC_);
  tree_->Branch("tr_bp", &trBP_);
  tree_->Branch("tr_E", &trE_);
  tree_->Branch("rh_id", &rhId_);
  tree_->Branch("rh_E", &rhE_);
  tree_->Branch("rh_x", &rhX_);
  tree_->Branch("rh_y", &rhY_);
  tree_->Branch("rh_z", &rhZ_);
  tree_->Branch("rh_layer", &rhLayer_);
  tree_->Branch("rh_det", &rhDet_);
  tree_->Branch("rh_lc", &rhLC_);
  tree_->Branch("rh_noTruthE", &rhNoTruthE_);
  tree_->Branch("lch_lc", &lchLC_);
  tree_->Branch("lch_rh", &lchRH_);
  tree_->Branch("lch_frac", &lchFrac_);
  tree_->Branch("at_id", &atId_);
  tree_->Branch("at_pdg", &atPdg_);
  tree_->Branch("at_parent", &atParent_);
  tree_->Branch("at_unit", &atUnit_);
  tree_->Branch("at_signal", &atSignal_);
  tree_->Branch("at_bx", &atBX_);
  tree_->Branch("at_evt", &atEvt_);
  tree_->Branch("at_depE", &atDepE_);
  tree_->Branch("tra_rh", &traRH_);
  tree_->Branch("tra_at", &traAT_);
  tree_->Branch("tra_E", &traE_);
  tree_->Branch("un_id", &unId_);
  tree_->Branch("un_pdg", &unPdg_);
  tree_->Branch("un_kind", &unKind_);
  tree_->Branch("un_origin", &unOrigin_);
  tree_->Branch("un_fsRow", &unFsRow_);
  tree_->Branch("un_signal", &unSignal_);
  tree_->Branch("un_bx", &unBX_);
  tree_->Branch("un_evt", &unEvt_);
  tree_->Branch("un_E", &unE_);
  tree_->Branch("un_eta", &unEta_);
  tree_->Branch("un_phi", &unPhi_);
  tree_->Branch("un_depE", &unDepE_);
  tree_->Branch("bp_id", &bpId_);
  tree_->Branch("bp_pdg", &bpPdg_);
  tree_->Branch("bp_signal", &bpSignal_);
  tree_->Branch("bp_bx", &bpBX_);
  tree_->Branch("bp_evt", &bpEvt_);
  tree_->Branch("bp_kind", &bpKind_);
  tree_->Branch("bp_origin", &bpOrigin_);
  tree_->Branch("bp_originPdg", &bpOriginPdg_);
  tree_->Branch("bp_originE", &bpOriginE_);
  tree_->Branch("bp_E", &bpE_);
  tree_->Branch("bp_eta", &bpEta_);
  tree_->Branch("bp_phi", &bpPhi_);
  tree_->Branch("bp_simE", &bpSimE_);
  tree_->Branch("bp_x", &bpX_);
  tree_->Branch("bp_y", &bpY_);
  tree_->Branch("bp_z", &bpZ_);
  tree_->Branch("bp_originRow", &bpOriginRow_);
  tree_->Branch("fs_id", &fsId_);
  tree_->Branch("fs_pdg", &fsPdg_);
  tree_->Branch("fs_E", &fsE_);
  tree_->Branch("fs_pt", &fsPt_);
  tree_->Branch("fs_eta", &fsEta_);
  tree_->Branch("fs_phi", &fsPhi_);
  if (!assignmentTag_.label().empty())
    tree_->Branch("clue_assignment", &clueAssignment_);
  for (std::size_t s = 0; s < trackstersTags_.size(); ++s) {
    const std::string n = "ts_" + trackstersTags_[s].label();
    auto& b = tracksters_[s];
    tree_->Branch((n + "_E").c_str(), &b.energy);
    tree_->Branch((n + "_nLC").c_str(), &b.nLC);
    tree_->Branch((n + "_lc").c_str(), &b.lc);
    tree_->Branch((n + "_mult").c_str(), &b.mult);
    if (tracksterDetails_) {
      tree_->Branch((n + "_regE").c_str(), &b.regE);
      tree_->Branch((n + "_bx").c_str(), &b.bx);
      tree_->Branch((n + "_by").c_str(), &b.by);
      tree_->Branch((n + "_bz").c_str(), &b.bz);
      tree_->Branch((n + "_time").c_str(), &b.time);
      tree_->Branch((n + "_timeErr").c_str(), &b.timeErr);
      tree_->Branch((n + "_prob").c_str(), &b.prob);
      tree_->Branch((n + "_trk").c_str(), &b.trk);
    }
  }
  if (!candidatesTag_.label().empty()) {
    tree_->Branch("cand_pdg", &candPdg_);
    tree_->Branch("cand_charge", &candCharge_);
    tree_->Branch("cand_E", &candE_);
    tree_->Branch("cand_pt", &candPt_);
    tree_->Branch("cand_eta", &candEta_);
    tree_->Branch("cand_phi", &candPhi_);
    tree_->Branch("cand_rawE", &candRawE_);
    tree_->Branch("cand_time", &candTime_);
    tree_->Branch("cand_timeErr", &candTimeErr_);
    tree_->Branch("cand_prob", &candProb_);
    tree_->Branch("cand_nTs", &candNTs_);
    tree_->Branch("cand_ts", &candTs_);
    tree_->Branch("cand_trk", &candTrk_);
  }
  if (!tracksTag_.label().empty()) {
    tree_->Branch("trk_idx", &trkIdx_);
    tree_->Branch("trk_pt", &trkPt_);
    tree_->Branch("trk_eta", &trkEta_);
    tree_->Branch("trk_phi", &trkPhi_);
    tree_->Branch("trk_charge", &trkCharge_);
    tree_->Branch("trk_ptErr", &trkPtErr_);
    tree_->Branch("trk_chi2n", &trkChi2n_);
    tree_->Branch("trk_nhits", &trkNHits_);
    tree_->Branch("trk_highPurity", &trkHighPurity_);
  }
  if (!genJetsTag_.label().empty()) {
    tree_->Branch("gj_pt", &gjPt_);
    tree_->Branch("gj_eta", &gjEta_);
    tree_->Branch("gj_phi", &gjPhi_);
    tree_->Branch("gj_E", &gjE_);
  }
}

void TracksterTruthNtuplizer::analyze(edm::Event const& event, edm::EventSetup const& setup) {
  rhtools_.setGeometry(setup.getData(geomToken_));
  auto const& graph = event.get(graphToken_);
  auto const& hitIndex = event.get(hitIndexToken_);
  auto const& lcs = event.get(layerClustersToken_);
  auto const& mask = event.get(maskToken_);
  const uint32_t nP = graph.nParticles();
  run_ = event.id().run();
  lumi_ = event.id().luminosityBlock();
  event_ = event.id().event();

  // --- base particle and origin of every truth particle -------------------------------------------
  std::vector<int> parent(nP);
  for (uint32_t p = 0; p < nP; ++p)
    parent[p] = firstParent(graph, p);
  auto crossed = [&](uint32_t p) {
    return !graph.particles()[p].backscattered && truth::Particle(&graph, p).checkpoint(0).has_value();
  };
  auto nearestCrossing = [&](uint32_t p) {
    for (int u = static_cast<int>(p); u >= 0; u = parent[u])
      if (crossed(u))
        return u;
    return kNone;
  };
  auto nearestAt = [&](uint32_t p, truth::LevelFlag flag) {
    for (int u = static_cast<int>(p); u >= 0; u = parent[u])
      if (graph.particles()[u].isAtLevel(flag))
        return u;
    return kNone;
  };
  auto rootOf = [&](uint32_t p) {
    int u = static_cast<int>(p);
    while (parent[u] >= 0)
      u = parent[u];
    return u;
  };
  std::vector<int> baseOf(nP, kNone), kindOf(nP, kRoot);
  for (uint32_t p = 0; p < nP; ++p) {
    if (int b = crossing_ ? nearestCrossing(p) : nearestAt(p, truth::LevelFlag::CaloBoundary); b >= 0) {
      baseOf[p] = b;
      kindOf[p] = kCaloBoundary;
    } else if (int r = nearestAt(p, truth::LevelFlag::ReconstructableFinalState); r >= 0) {
      baseOf[p] = r;
      kindOf[p] = kReconstructableFinalState;
    } else {
      baseOf[p] = rootOf(p);
      kindOf[p] = kRoot;
    }
  }

  // --- in-time sim energy per cell, per base particle ---------------------------------------------
  std::unordered_map<int, int> bpIndex;  // graph particle id -> row of the base-particle table
  std::vector<double> bpSim;
  std::unordered_map<uint32_t, std::vector<std::pair<int, float>>> cellSim;
  std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, float>>> cellSimAtom;  // per graph particle
  std::unordered_map<uint32_t, double> cellTotal;
  for (uint32_t p = 0; p < nP; ++p) {
    auto hits = hitIndex.directHits(truth::HitChannel::Calo, p);
    if (hits.empty())
      continue;
    const int b = baseOf[p];
    auto [it, inserted] = bpIndex.emplace(b, static_cast<int>(bpIndex.size()));
    if (inserted)
      bpSim.push_back(0.);
    const int row = it->second;
    for (auto const& h : hits) {
      if (h.energy <= 0.f)
        continue;
      cellSim[h.detId].emplace_back(row, h.energy);
      cellSimAtom[h.detId].emplace_back(p, h.energy);
      cellTotal[h.detId] += h.energy;
      bpSim[row] += h.energy;
    }
  }

  std::unordered_map<uint32_t, float> recHitEnergy;
  for (auto const& token : recHitTokens_)
    for (auto const& rh : event.get(token))
      recHitEnergy[rh.detid().rawId()] = rh.energy();

  // --- layer clusters and the truth table --------------------------------------------------------
  for (auto* v : {&lcE_, &lcX_, &lcY_, &lcZ_, &lcMask_, &lcNoTruthE_, &lcRecE_, &trE_})
    v->clear();
  for (auto* v : {&lcLayer_, &lcSide_, &lcDet_, &lcNHits_, &lcAlgo_, &lcMissingRecHits_, &trLC_, &trBP_})
    v->clear();
  lcSeed_.clear();
  std::unordered_map<uint32_t, std::pair<int, float>> rhBestLC;  // rechit -> layer cluster with the largest fraction
  std::vector<std::tuple<int, uint32_t, float>> lchRaw;            // (layer cluster, rechit DetId, fraction)
  std::vector<double> accum(bpIndex.size(), 0.);
  std::vector<int> touched;
  for (std::size_t l = 0; l < lcs.size(); ++l) {
    auto const& lc = lcs[l];
    const DetId seed = lc.seed();
    lcE_.push_back(lc.energy());
    lcX_.push_back(lc.x());
    lcY_.push_back(lc.y());
    lcZ_.push_back(lc.z());
    lcLayer_.push_back(static_cast<int>(rhtools_.getLayerWithOffset(seed)));
    lcSide_.push_back(lc.z() > 0 ? 1 : 0);
    lcDet_.push_back(static_cast<int>(seed.det()));
    lcNHits_.push_back(static_cast<int>(lc.hitsAndFractions().size()));
    lcAlgo_.push_back(static_cast<int>(lc.algo()));
    lcSeed_.push_back(seed.rawId());
    lcMask_.push_back(l < mask.size() ? mask[l] : 0.f);

    double noTruth = 0., recSum = 0.;
    int missing = 0;
    touched.clear();
    for (auto const& [detId, fraction] : lc.hitsAndFractions()) {
      auto rh = recHitEnergy.find(detId.rawId());
      if (rh == recHitEnergy.end()) {
        ++missing;
        continue;
      }
      const double share = static_cast<double>(fraction) * rh->second;
      recSum += share;
      lchRaw.emplace_back(static_cast<int>(l), detId.rawId(), fraction);
      if (auto [b, ins] = rhBestLC.emplace(detId.rawId(), std::make_pair(static_cast<int>(l), fraction));
          !ins && fraction > b->second.second)
        b->second = {static_cast<int>(l), fraction};
      auto tot = cellTotal.find(detId.rawId());
      if (tot == cellTotal.end() || tot->second <= 0.) {
        noTruth += share;
        continue;
      }
      for (auto const& [row, e] : cellSim[detId.rawId()]) {
        if (accum[row] == 0.)
          touched.push_back(row);
        accum[row] += share * e / tot->second;
      }
    }
    std::sort(touched.begin(), touched.end());
    for (int row : touched) {
      trLC_.push_back(static_cast<int>(l));
      trBP_.push_back(row);
      trE_.push_back(static_cast<float>(accum[row]));
      accum[row] = 0.;
    }
    lcNoTruthE_.push_back(static_cast<float>(noTruth));
    lcRecE_.push_back(static_cast<float>(recSum));
    lcMissingRecHits_.push_back(missing);
  }

  // --- every HGCAL rechit (sorted by DetId: reproducible order), the layer-cluster membership, -----
  // --- and the (rechit, atom) table; atoms = graph particles with sim energy in a rechit cell ------
  for (auto* v : {&rhE_, &rhX_, &rhY_, &rhZ_, &rhNoTruthE_, &lchFrac_, &traE_})
    v->clear();
  for (auto* v : {&rhLayer_, &rhDet_, &rhLC_, &lchLC_, &lchRH_, &traRH_, &traAT_})
    v->clear();
  rhId_.clear();
  std::vector<uint32_t> rhIds;
  rhIds.reserve(recHitEnergy.size());
  for (auto const& [id, e] : recHitEnergy)
    rhIds.push_back(id);
  std::sort(rhIds.begin(), rhIds.end());
  std::vector<uint32_t> atoms;
  for (uint32_t id : rhIds)
    if (auto c = cellSimAtom.find(id); c != cellSimAtom.end())
      for (auto const& ps : c->second)
        atoms.push_back(ps.first);
  std::sort(atoms.begin(), atoms.end());
  atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
  std::unordered_map<uint32_t, int> atRow, rhRow;
  for (std::size_t k = 0; k < atoms.size(); ++k)
    atRow.emplace(atoms[k], static_cast<int>(k));
  rhRow.reserve(rhIds.size());
  std::vector<double> atAccum(atoms.size(), 0.), atDep(atoms.size(), 0.);
  std::vector<char> atTouched(atoms.size(), 0);
  for (uint32_t id : rhIds) {
    const DetId det(id);
    const float e = recHitEnergy[id];
    const auto pos = rhtools_.getPosition(det);
    const int row = static_cast<int>(rhE_.size());
    rhRow.emplace(id, row);
    rhId_.push_back(id);
    rhE_.push_back(e);
    rhX_.push_back(pos.x());
    rhY_.push_back(pos.y());
    rhZ_.push_back(pos.z());
    rhLayer_.push_back(static_cast<int>(rhtools_.getLayerWithOffset(det)));
    rhDet_.push_back(static_cast<int>(det.det()));
    auto best = rhBestLC.find(id);
    rhLC_.push_back(best == rhBestLC.end() ? kNone : best->second.first);
    auto tot = cellTotal.find(id);
    if (tot == cellTotal.end() || tot->second <= 0.) {
      rhNoTruthE_.push_back(e);
      continue;
    }
    rhNoTruthE_.push_back(0.f);
    touched.clear();
    for (auto const& [p, s] : cellSimAtom[id]) {
      const int a = atRow[p];
      if (!atTouched[a]) {
        atTouched[a] = 1;
        touched.push_back(a);
      }
      atAccum[a] += static_cast<double>(e) * s / tot->second;
    }
    std::sort(touched.begin(), touched.end());
    for (int a : touched) {
      traRH_.push_back(row);
      traAT_.push_back(a);
      traE_.push_back(static_cast<float>(atAccum[a]));
      atDep[a] += atAccum[a];
      atAccum[a] = 0.;
      atTouched[a] = 0;
    }
  }
  for (auto const& [l, id, fraction] : lchRaw) {
    lchLC_.push_back(l);
    lchRH_.push_back(rhRow.at(id));
    lchFrac_.push_back(fraction);
  }

  // --- signal final-state table (the origin level of the base particles) -------------------------
  for (auto* v : {&fsId_, &fsPdg_})
    v->clear();
  for (auto* v : {&fsE_, &fsPt_, &fsEta_, &fsPhi_})
    v->clear();
  std::unordered_map<int, int> fsRow;  // graph particle id -> row of the final-state table
  for (uint32_t p = 0; p < nP; ++p) {
    auto const& pd = graph.particles()[p];
    if (!pd.isSignal() || !pd.isAtLevel(truth::LevelFlag::ReconstructableFinalState))
      continue;
    fsRow.emplace(static_cast<int>(p), static_cast<int>(fsId_.size()));
    fsId_.push_back(static_cast<int>(p));
    fsPdg_.push_back(pd.pdgId);
    fsE_.push_back(pd.momentum.energy());
    fsPt_.push_back(pd.momentum.pt());
    fsEta_.push_back(pd.momentum.pt() > 0 ? pd.momentum.eta() : 0.f);
    fsPhi_.push_back(pd.momentum.phi());
  }

  // --- atoms and units ---------------------------------------------------------------------------
  // unit of a particle: its nearest reconstructableFinalState ancestor-or-self; below a pi0 the pi0's decay
  // daughter it descends from (downstream reconstructs photons, not pi0s); with no such ancestor, its root
  auto unitOf = [&](uint32_t p) -> std::pair<int, int> {
    const int o = nearestAt(p, truth::LevelFlag::ReconstructableFinalState);
    if (o < 0)
      return {rootOf(p), kUnitRoot};
    if (std::abs(graph.particles()[o].pdgId) == kPi0 && static_cast<int>(p) != o) {
      int u = static_cast<int>(p);
      while (parent[u] >= 0 && parent[u] != o)
        u = parent[u];
      if (parent[u] == o)
        return {u, kUnitPi0Daughter};
    }
    return {o, kUnitFinalState};
  };
  for (auto* v : {&atId_, &atPdg_, &atParent_, &atUnit_, &atSignal_, &atBX_, &atEvt_})
    v->clear();
  for (auto* v : {&unId_, &unPdg_, &unKind_, &unOrigin_, &unFsRow_, &unSignal_, &unBX_, &unEvt_})
    v->clear();
  for (auto* v : {&atDepE_, &unE_, &unEta_, &unPhi_, &unDepE_})
    v->clear();
  std::vector<std::pair<int, int>> atUnit(atoms.size());
  std::vector<int> unitIds;
  for (std::size_t k = 0; k < atoms.size(); ++k) {
    atUnit[k] = unitOf(atoms[k]);
    unitIds.push_back(atUnit[k].first);
  }
  std::sort(unitIds.begin(), unitIds.end());
  unitIds.erase(std::unique(unitIds.begin(), unitIds.end()), unitIds.end());
  std::unordered_map<int, int> unRow, unKindOf;
  for (std::size_t u = 0; u < unitIds.size(); ++u)
    unRow.emplace(unitIds[u], static_cast<int>(u));
  std::vector<double> unDep(unitIds.size(), 0.);
  for (std::size_t k = 0; k < atoms.size(); ++k) {
    const uint32_t p = atoms[k];
    auto const& pd = graph.particles()[p];
    int q = parent[p];
    while (q >= 0 && !atRow.count(static_cast<uint32_t>(q)))
      q = parent[q];
    const int u = unRow.at(atUnit[k].first);
    unKindOf[u] = atUnit[k].second;
    atId_.push_back(static_cast<int>(p));
    atPdg_.push_back(pd.pdgId);
    atParent_.push_back(q >= 0 ? atRow.at(static_cast<uint32_t>(q)) : kNone);
    atUnit_.push_back(u);
    atSignal_.push_back(pd.isSignal() ? 1 : 0);
    atBX_.push_back(pd.bunchCrossing());
    atEvt_.push_back(pd.eventIndex());
    atDepE_.push_back(static_cast<float>(atDep[k]));
    unDep[u] += atDep[k];
  }
  for (std::size_t u = 0; u < unitIds.size(); ++u) {
    const int gid = unitIds[u];
    auto const& pd = graph.particles()[gid];
    const int origin = nearestAt(gid, truth::LevelFlag::ReconstructableFinalState);
    auto f = fsRow.find(origin);
    unId_.push_back(gid);
    unPdg_.push_back(pd.pdgId);
    unKind_.push_back(unKindOf[static_cast<int>(u)]);
    unOrigin_.push_back(origin);
    unFsRow_.push_back(f == fsRow.end() ? kNone : f->second);
    unSignal_.push_back(pd.isSignal() ? 1 : 0);
    unBX_.push_back(pd.bunchCrossing());
    unEvt_.push_back(pd.eventIndex());
    unE_.push_back(pd.momentum.energy());
    unEta_.push_back(pd.momentum.pt() > 0 ? pd.momentum.eta() : 0.f);
    unPhi_.push_back(pd.momentum.phi());
    unDepE_.push_back(static_cast<float>(unDep[u]));
  }

  // --- base-particle table --------------------------------------------------------------------------
  for (auto* v : {&bpId_, &bpPdg_, &bpSignal_, &bpBX_, &bpEvt_, &bpKind_, &bpOrigin_, &bpOriginPdg_, &bpOriginRow_})
    v->assign(bpIndex.size(), kNone);
  for (auto* v : {&bpE_, &bpEta_, &bpPhi_, &bpSimE_, &bpOriginE_, &bpX_, &bpY_, &bpZ_})
    v->assign(bpIndex.size(), 0.f);
  for (auto const& [gid, row] : bpIndex) {
    auto const& pd = graph.particles()[gid];
    bpId_[row] = gid;
    bpPdg_[row] = pd.pdgId;
    bpSignal_[row] = pd.isSignal() ? 1 : 0;
    bpBX_[row] = pd.bunchCrossing();
    bpEvt_[row] = pd.eventIndex();
    bpKind_[row] = (crossing_ ? crossed(gid) : pd.isAtLevel(truth::LevelFlag::CaloBoundary)) ? kCaloBoundary
                   : pd.isAtLevel(truth::LevelFlag::ReconstructableFinalState)               ? kReconstructableFinalState
                                                                                             : kRoot;
    // energy and direction where it entered the calorimeter, if recorded
    if (auto cp = truth::Particle(&graph, gid).checkpoint(0)) {
      bpE_[row] = cp->momentum.energy();
      bpEta_[row] = cp->position.eta();
      bpPhi_[row] = cp->position.phi();
      bpX_[row] = cp->position.x();
      bpY_[row] = cp->position.y();
      bpZ_[row] = cp->position.z();
    } else {
      bpE_[row] = pd.momentum.energy();
      bpEta_[row] = pd.momentum.pt() > 0 ? pd.momentum.eta() : 0.f;
      bpPhi_[row] = pd.momentum.phi();
    }
    bpSimE_[row] = static_cast<float>(bpSim[row]);
    if (int o = nearestAt(gid, truth::LevelFlag::ReconstructableFinalState); o >= 0) {
      bpOrigin_[row] = o;
      bpOriginPdg_[row] = graph.particles()[o].pdgId;
      bpOriginE_[row] = graph.particles()[o].momentum.energy();
      if (auto f = fsRow.find(o); f != fsRow.end())
        bpOriginRow_[row] = f->second;
    }
  }

  // --- CMSSW reference ----------------------------------------------------------------------------
  if (!assignmentTag_.label().empty()) {
    auto const& a = event.get(assignmentToken_);
    clueAssignment_.assign(a.begin(), a.end());
  }
  for (std::size_t s = 0; s < trackstersTokens_.size(); ++s) {
    auto& out = tracksters_[s];
    out.clear();
    for (auto const& t : event.get(trackstersTokens_[s])) {
      out.energy.push_back(t.raw_energy());
      out.nLC.push_back(static_cast<int>(t.vertices().size()));
      for (std::size_t k = 0; k < t.vertices().size(); ++k) {
        out.lc.push_back(static_cast<int>(t.vertices()[k]));
        out.mult.push_back(k < t.vertex_multiplicity().size() ? t.vertex_multiplicity()[k] : 1.f);
      }
      if (tracksterDetails_) {
        out.regE.push_back(t.regressed_energy());
        out.bx.push_back(t.barycenter().x());
        out.by.push_back(t.barycenter().y());
        out.bz.push_back(t.barycenter().z());
        out.time.push_back(t.time());
        out.timeErr.push_back(t.timeError());
        for (float p : t.id_probabilities())
          out.prob.push_back(p);
        out.trk.push_back(t.trackIdx());
      }
    }
  }

  // --- TICL candidates: tracksters are keys into the candidate producer's merged collection -------
  if (!candidatesTag_.label().empty()) {
    for (auto* v : {&candPdg_, &candCharge_, &candNTs_, &candTs_, &candTrk_})
      v->clear();
    for (auto* v : {&candE_, &candPt_, &candEta_, &candPhi_, &candRawE_, &candTime_, &candTimeErr_, &candProb_})
      v->clear();
    for (auto const& c : event.get(candidatesToken_)) {
      candPdg_.push_back(c.pdgId());
      candCharge_.push_back(c.charge());
      candE_.push_back(c.energy());
      candPt_.push_back(c.pt());
      candEta_.push_back(c.eta());
      candPhi_.push_back(c.phi());
      candRawE_.push_back(c.rawEnergy());
      candTime_.push_back(c.time());
      candTimeErr_.push_back(c.timeError());
      for (float p : c.idProbabilities())
        candProb_.push_back(p);
      candNTs_.push_back(static_cast<int>(c.tracksters().size()));
      for (auto const& t : c.tracksters())
        candTs_.push_back(static_cast<int>(t.key()));
      candTrk_.push_back(c.trackPtrs().empty() ? kNone : static_cast<int>(c.trackPtrs()[0].key()));
    }
  }

  // --- tracks (the candidate's track index refers to the full collection: trk_idx keeps it) -------
  if (!tracksTag_.label().empty()) {
    for (auto* v : {&trkIdx_, &trkCharge_, &trkNHits_, &trkHighPurity_})
      v->clear();
    for (auto* v : {&trkPt_, &trkEta_, &trkPhi_, &trkChi2n_, &trkPtErr_})
      v->clear();
    auto const& tracks = event.get(tracksToken_);
    for (std::size_t i = 0; i < tracks.size(); ++i) {
      auto const& t = tracks[i];
      if (std::abs(t.eta()) < trackMinAbsEta_ || t.pt() < trackMinPt_)
        continue;
      trkIdx_.push_back(static_cast<int>(i));
      trkPt_.push_back(t.pt());
      trkEta_.push_back(t.eta());
      trkPhi_.push_back(t.phi());
      trkCharge_.push_back(t.charge());
      trkPtErr_.push_back(t.ptError());
      trkChi2n_.push_back(t.normalizedChi2());
      trkNHits_.push_back(static_cast<int>(t.numberOfValidHits()));
      trkHighPurity_.push_back(t.quality(reco::TrackBase::highPurity) ? 1 : 0);
    }
  }

  // --- gen jets -----------------------------------------------------------------------------------
  if (!genJetsTag_.label().empty()) {
    for (auto* v : {&gjPt_, &gjEta_, &gjPhi_, &gjE_})
      v->clear();
    for (auto const& j : event.get(genJetsToken_)) {
      if (j.pt() < genJetMinPt_)
        continue;
      gjPt_.push_back(j.pt());
      gjEta_.push_back(j.eta());
      gjPhi_.push_back(j.phi());
      gjE_.push_back(j.energy());
    }
  }
  tree_->Fill();
}

void TracksterTruthNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("graph", edm::InputTag("truthLogicalGraphProducer"));
  desc.add<edm::InputTag>("hitIndex", edm::InputTag("truthLogicalGraphHitIndexProducer"));
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<std::string>("baseLevel", "crossing")
      ->setComment("crossing: hits go to the nearest ancestor that crossed into the calorimeter (every crossing); "
                   "caloBoundary: to the nearest member of the caloBoundary level (outermost crossings only)");
  desc.add<edm::InputTag>("layerClusterMask", edm::InputTag("filteredLayerClustersCLUE3DHigh", "CLUE3DHigh"))
      ->setComment("The mask of the trackster-building iteration: the layer clusters it may use");
  desc.add<edm::InputTag>("clueAssignment", edm::InputTag("ticlTrackstersCLUEsteringAssignment"))
      ->setComment("CLUEstering trackster index per layer cluster (-1 = outlier or masked); empty label = skip");
  desc.add<std::vector<edm::InputTag>>("tracksters", {edm::InputTag("ticlTrackstersCLUE3DHigh")});
  desc.add<bool>("tracksterDetails", true)
      ->setComment("also store regressed energy, barycenter, time, PID probabilities and track index per trackster");
  desc.add<edm::InputTag>("candidates", edm::InputTag(""))
      ->setComment("TICLCandidate collection (e.g. ticlCandidate); empty label = skip. cand_ts are keys into the "
                   "candidate producer's trackster collection: list it in `tracksters` under the same label");
  desc.add<edm::InputTag>("tracks", edm::InputTag(""))->setComment("track collection (e.g. generalTracks); empty = skip");
  desc.add<double>("trackMinAbsEta", 1.2);
  desc.add<double>("trackMinPt", 0.5);
  desc.add<edm::InputTag>("genJets", edm::InputTag(""))->setComment("gen jet collection (e.g. ak4GenJetsNoNu); empty = skip");
  desc.add<double>("genJetMinPt", 5.);
  desc.add<std::vector<edm::InputTag>>("recHits",
                                       {edm::InputTag("HGCalRecHit", "HGCEERecHits"),
                                        edm::InputTag("HGCalRecHit", "HGCHEFRecHits"),
                                        edm::InputTag("HGCalRecHit", "HGCHEBRecHits")});
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(TracksterTruthNtuplizer);
