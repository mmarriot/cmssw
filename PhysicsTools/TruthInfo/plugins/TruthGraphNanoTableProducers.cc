// Flat nanoAOD tables of the MC-truth graph itself.
//
// TruthGraphFlatTableProducer writes a selection of the particles of the logical graph (the
// signal interaction, every particle at chosen truth levels, pile-up included, every particle
// a reco object is matched to, and their ancestors), the vertices that produce or end them
// and the truth vertices of the vertex association maps, and the parent -> child links
// between the particles. It publishes the graph ids of its particle and
// vertex rows, in row order, which is ascending graph id: the key space with which
// TruthBranchAssociationFlatTableProducer and TruthBranchVertexAssociationFlatTableProducer
// (TruthBranchNanoTableProducers.cc) remap the association maps of
// SimGeneral/TruthGraphAssociatorProducers, so their links index these tables directly.
//
// TruthGraphHitFlatTableProducer writes, for the exported particles at chosen truth levels,
// the simulated energy their subgraph (the particle and its descendants) left in each hit of
// a hit table, from truth::LogicalGraphHitIndex. The hits are found through a uint32 DetId
// column of that table, so any hit table with one works.
//
// A missing input gives empty tables with all their columns (skipNonExistingSrc, as the
// table producers of TruthBranchNanoTableProducers.cc), so a NanoAOD output books the same
// branches whichever event comes first.

#include <algorithm>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/allowedValues.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/transform.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"
#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/TruthInfo/interface/LogicalGraphHitIndex.h"
#include "SimDataFormats/TruthInfo/interface/Particle.h"

namespace {

  // The particle's endpoint: of its decay vertices, the latest in time, the higher graph id
  // on a tie; -1 for none. A Geant4 track that interacts several times keeps several decay
  // vertices, in no particular order.
  int32_t endVertex(truth::Graph const& g, uint32_t particle) {
    bool found = false;
    uint32_t best = 0;
    for (const uint32_t v : g.decayVertices(particle)) {
      if (v >= g.nVertices())
        continue;
      const double t = g.vertices()[v].position.t();
      const double tBest = g.vertices()[best].position.t();
      if (!found || t > tBest || (t == tBest && v > best))
        best = v;
      found = true;
    }
    return found ? static_cast<int32_t>(best) : -1;
  }

  // The first production vertex, -1 for none. The graphs of the validation samples give a
  // particle at most one.
  int32_t productionVertex(truth::Graph const& g, uint32_t particle) {
    const auto prod = g.productionVertices(particle);
    return !prod.empty() && prod[0] < g.nVertices() ? static_cast<int32_t>(prod[0]) : -1;
  }

}  // namespace

class TruthGraphFlatTableProducer : public edm::global::EDProducer<> {
public:
  using SharedEnergyMap = ticl::TICLAssociationMap<ticl::mapWithSharedEnergyAndScore>;
  using FractionMap = ticl::TICLAssociationMap<ticl::mapWithFractionAndScore>;

  explicit TruthGraphFlatTableProducer(edm::ParameterSet const& cfg)
      : graphToken_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        recoToTruthTokens_(edm::vector_transform(cfg.getParameter<std::vector<edm::InputTag>>("keepMatchedBy"),
                                                 [this](edm::InputTag const& tag) { return mapToken(tag); })),
        truthToRecoTokens_(
            edm::vector_transform(cfg.getParameter<std::vector<edm::InputTag>>("keepMatchedByTruthToReco"),
                                  [this](edm::InputTag const& tag) { return mapToken(tag); })),
        recoToTruthVertexTokens_(
            edm::vector_transform(cfg.getParameter<std::vector<edm::InputTag>>("keepVerticesMatchedBy"),
                                  [this](edm::InputTag const& tag) { return consumes<FractionMap>(tag); })),
        truthToRecoVertexTokens_(
            edm::vector_transform(cfg.getParameter<std::vector<edm::InputTag>>("keepVerticesMatchedByTruthToReco"),
                                  [this](edm::InputTag const& tag) { return consumes<FractionMap>(tag); })),
        name_(cfg.getParameter<std::string>("name")),
        doc_(cfg.getParameter<std::string>("doc")),
        vertexName_(cfg.getParameter<std::string>("vertexName")),
        childName_(cfg.getParameter<std::string>("childName")),
        keepSignal_(cfg.getParameter<bool>("keepSignal")),
        keepAncestors_(cfg.getParameter<bool>("keepAncestors")),
        keepLevelMask_(cfg.getParameter<unsigned int>("keepLevelMask")),
        skipNonExistingSrc_(cfg.getParameter<bool>("skipNonExistingSrc")),
        precision_(cfg.getParameter<int>("precision")),
        positionPrecision_(cfg.getParameter<int>("positionPrecision")),
        particleToken_(produces<nanoaod::FlatTable>()),
        vertexToken_(produces<nanoaod::FlatTable>("vertices")),
        childToken_(produces<nanoaod::FlatTable>("children")),
        selectedToken_(produces<std::vector<unsigned int>>("selected")),
        selectedVerticesToken_(produces<std::vector<unsigned int>>("selectedVertices")) {}

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  edm::EDGetTokenT<SharedEnergyMap> mapToken(edm::InputTag const& tag) { return consumes<SharedEnergyMap>(tag); }

  const edm::EDGetTokenT<truth::Graph> graphToken_;
  const std::vector<edm::EDGetTokenT<SharedEnergyMap>> recoToTruthTokens_, truthToRecoTokens_;
  const std::vector<edm::EDGetTokenT<FractionMap>> recoToTruthVertexTokens_, truthToRecoVertexTokens_;
  const std::string name_, doc_, vertexName_, childName_;
  const bool keepSignal_, keepAncestors_;
  const uint32_t keepLevelMask_;
  const bool skipNonExistingSrc_;
  const int precision_, positionPrecision_;
  const edm::EDPutTokenT<nanoaod::FlatTable> particleToken_, vertexToken_, childToken_;
  const edm::EDPutTokenT<std::vector<unsigned int>> selectedToken_, selectedVerticesToken_;
};

void TruthGraphFlatTableProducer::produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const {
  // Without a graph the selection is empty and the tables keep their columns.
  static const truth::Graph noGraph;
  auto const graph = event.getHandle(graphToken_);
  if (!graph.isValid() && !skipNonExistingSrc_)
    throw cms::Exception("ProductNotFound")
        << "TruthGraphFlatTableProducer '" << name_ << "': the truth graph is missing";
  auto const& g = graph.isValid() ? *graph : noGraph;
  auto const& particles = g.particles();
  const uint32_t n = g.nParticles();
  const uint32_t nv = g.nVertices();

  // Selection: the signal interaction, the particles at the keepLevelMask levels, the
  // particles a reco object is matched to in either direction, then the ancestors of all.
  std::vector<uint8_t> keep(n, 0);
  for (uint32_t i = 0; i < n; ++i)
    keep[i] = (keepSignal_ && particles[i].isSignal()) || (particles[i].levelFlags & keepLevelMask_) != 0;
  // Marks the truth objects (particles, or vertices for the vertex maps) the maps of `tokens`
  // match: the entries of a reco -> truth map, the non-empty rows of a truth -> reco map.
  auto const markMatched = [&](auto const& tokens, bool truthToReco, std::vector<uint8_t>& marked) {
    const size_t nTruth = marked.size();
    for (auto const& token : tokens) {
      auto const map = event.getHandle(token);
      if (!map.isValid()) {
        if (!skipNonExistingSrc_)
          throw cms::Exception("ProductNotFound") << "TruthGraphFlatTableProducer '" << name_
                                                  << "': an association map to keep matched truth objects by is "
                                                  << "missing";
        continue;
      }
      auto const& rows = map->getMap();
      if (truthToReco) {  // a row per truth object
        const size_t nRows = std::min(rows.size(), nTruth);
        for (size_t id = 0; id < nRows; ++id)
          marked[id] |= !rows[id].empty();
      } else {  // a row per reco object
        for (auto const& row : rows)
          for (auto const& element : row)
            if (element.index() < nTruth)
              marked[element.index()] = 1;
      }
    }
  };
  markMatched(recoToTruthTokens_, false, keep);
  markMatched(truthToRecoTokens_, true, keep);
  if (keepAncestors_) {
    std::vector<uint32_t> todo;
    for (uint32_t i = 0; i < n; ++i)
      if (keep[i])
        todo.push_back(i);
    while (!todo.empty()) {
      const uint32_t i = todo.back();
      todo.pop_back();
      for (const uint32_t v : g.productionVertices(i))
        for (const uint32_t p : g.incomingParticles(v))
          if (p < n && !keep[p]) {
            keep[p] = 1;
            todo.push_back(p);
          }
    }
  }
  auto selected = std::make_unique<std::vector<unsigned int>>();  // graph id of each row, ascending
  std::vector<int32_t> rowOf(n, -1);
  for (uint32_t i = 0; i < n; ++i)
    if (keep[i]) {
      rowOf[i] = selected->size();
      selected->push_back(i);
    }

  // The vertex rows, in ascending graph id: the production vertex and the endpoint of every
  // exported particle, and the truth vertices the vertex maps match, so no vertex link is lost.
  const size_t m = selected->size();
  std::vector<int32_t> prodGraphId(m), endGraphId(m);
  std::vector<uint8_t> keepVertex(nv, 0);
  markMatched(recoToTruthVertexTokens_, false, keepVertex);
  markMatched(truthToRecoVertexTokens_, true, keepVertex);
  for (size_t r = 0; r < m; ++r) {
    prodGraphId[r] = productionVertex(g, (*selected)[r]);
    endGraphId[r] = endVertex(g, (*selected)[r]);
    for (const int32_t v : {prodGraphId[r], endGraphId[r]})
      if (v >= 0)
        keepVertex[v] = 1;
  }
  auto selectedVertices = std::make_unique<std::vector<unsigned int>>();  // graph id of each row, ascending
  std::vector<int32_t> vertexRowOf(nv, -1);
  for (uint32_t v = 0; v < nv; ++v)
    if (keepVertex[v]) {
      vertexRowOf[v] = selectedVertices->size();
      selectedVertices->push_back(v);
    }

  std::vector<uint32_t> levelFlags(m);
  std::vector<int32_t> pdgId(m), prodVertexIdx(m, -1), decayVertexIdx(m, -1), bunchCrossing(m), eventIndex(m);
  std::vector<int16_t> status(m);
  std::vector<uint16_t> statusFlags(m);
  std::vector<uint8_t> role(m);
  std::vector<bool> isSignal(m), hasGen(m), hasSim(m), hasBoundary(m);
  std::vector<float> charge(m), pt(m), eta(m), phi(m), energy(m);
  std::vector<float> boundaryX(m), boundaryY(m), boundaryZ(m), boundaryT(m);
  std::vector<float> boundaryPt(m), boundaryEta(m), boundaryPhi(m), boundaryE(m);
  std::vector<int32_t> parentIdx, childIdx;
  for (size_t r = 0; r < m; ++r) {
    const uint32_t i = (*selected)[r];
    auto const& p = particles[i];
    pdgId[r] = p.pdgId;
    status[r] = p.status;
    statusFlags[r] = p.statusFlags;
    levelFlags[r] = p.levelFlags;
    role[r] = p.role;
    isSignal[r] = p.isSignal();
    hasGen[r] = p.hasGen();
    hasSim[r] = p.hasSim();
    bunchCrossing[r] = p.bunchCrossing();
    eventIndex[r] = p.eventIndex();
    charge[r] = g.particle(i).charge();
    pt[r] = p.momentum.pt();
    eta[r] = p.momentum.eta();
    phi[r] = p.momentum.phi();
    energy[r] = p.momentum.E();
    if (prodGraphId[r] >= 0)
      prodVertexIdx[r] = vertexRowOf[prodGraphId[r]];
    if (endGraphId[r] >= 0)
      decayVertexIdx[r] = vertexRowOf[endGraphId[r]];
    // Checkpoint 0 is the tracker-calorimeter boundary crossing the SimTrack recorded
    // (TruthLogicalGraphProducer); a back-scattered track has one too. The columns follow
    // TruthParticleFlatTableProducer: zeros and hasBoundary false without one.
    if (const auto boundary = g.particle(i).checkpoint(0)) {
      auto const& pos = boundary->position;
      // The checkpoint is stored in float. ROOT's eta, log(z/rho + sqrt((z/rho)^2 + 1)),
      // cancels for negative z: in float it is off by ~1e-4 at eta = -4.6, so the angles
      // are computed in double.
      const math::XYZTLorentzVectorD mom(boundary->momentum);
      hasBoundary[r] = true;
      boundaryX[r] = pos.x();
      boundaryY[r] = pos.y();
      boundaryZ[r] = pos.z();
      boundaryT[r] = pos.t();
      boundaryE[r] = mom.energy();
      if (mom.px() != 0. || mom.py() != 0. || mom.pz() != 0.) {
        boundaryPt[r] = mom.pt();
        boundaryEta[r] = mom.eta();
        boundaryPhi[r] = mom.phi();
      }
    }
    for (const uint32_t v : g.decayVertices(i))
      if (v < nv)
        for (const uint32_t child : g.outgoingParticles(v))
          if (child < n && rowOf[child] >= 0) {
            parentIdx.push_back(r);
            childIdx.push_back(rowOf[child]);
          }
  }

  auto table = std::make_unique<nanoaod::FlatTable>(m, name_, false);
  table->setDoc(doc_);
  table->addColumn<uint32_t>("graphId", *selected, "Particle id in the truth graph (key of the association maps)");
  table->addColumn<int32_t>("pdgId", pdgId, "PDG id");
  table->addColumn<float>("charge", charge, "Electric charge [e]");
  table->addColumn<float>("pt", pt, "Transverse momentum [GeV]", precision_);
  table->addColumn<float>(
      "eta", eta, "Pseudorapidity (ROOT convention for pt = 0: +-(22756 + |pz|/GeV), 0 for p = 0)", precision_);
  table->addColumn<float>("phi", phi, "Azimuth [rad]", precision_);
  table->addColumn<float>("energy", energy, "Energy [GeV]", precision_);
  table->addColumn<int16_t>("status", status, "Generator status (0 for a SIM-only particle)");
  table->addColumn<uint16_t>("statusFlags", statusFlags, "Packed reco::GenStatusFlags bits (0 if not available)");
  table->addColumn<uint32_t>(
      "levelFlags", levelFlags, "Truth levels the particle belongs to (truth::LevelFlag bits; bit 3 CaloBoundary)");
  table->addColumn<bool>("isSignal", isSignal, "True for the signal interaction, false for pile-up");
  table->addColumn<int32_t>("bunchCrossing", bunchCrossing, "Bunch crossing of the interaction (0 = in time)");
  table->addColumn<int32_t>("eventIndex", eventIndex, "Interaction index inside the bunch crossing (0 = signal)");
  table->addColumn<bool>("hasGen", hasGen, "Has a generator particle");
  table->addColumn<bool>("hasSim", hasSim, "Has a Geant4 track (SimTrack)");
  table->addColumn<uint8_t>("role", role, "truth::ParticleRole: 0 Normal, 1 Connector, 2 SignalStandIn");
  table->addColumn<int32_t>("prodVertexIdx",
                            prodVertexIdx,
                            "Row of the production vertex in " + vertexName_ + " (its position is there), -1 if none");
  table->addColumn<int32_t>("decayVertexIdx",
                            decayVertexIdx,
                            "Row of the endpoint in " + vertexName_ +
                                ": of the decay or interaction vertices (where its saved secondaries were made, "
                                "not necessarily where it stopped), the latest in time; -1 if none. The children of "
                                "every end vertex are in " +
                                childName_);
  table->addColumn<bool>(
      "hasBoundary",
      hasBoundary,
      "Checkpoint at the tracker-calorimeter boundary recorded (truth::Checkpoint 0 of the SimTrack; the boundary "
      "columns are 0 without it)");
  table->addColumn<float>("boundaryX", boundaryX, "x at the calo boundary [cm]", positionPrecision_);
  table->addColumn<float>("boundaryY", boundaryY, "y at the calo boundary [cm]", positionPrecision_);
  table->addColumn<float>("boundaryZ", boundaryZ, "z at the calo boundary [cm]", positionPrecision_);
  table->addColumn<float>("boundaryT", boundaryT, "t at the calo boundary [ns]", positionPrecision_);
  table->addColumn<float>("boundaryPt", boundaryPt, "pT at the calo boundary [GeV]", precision_);
  table->addColumn<float>("boundaryEta",
                          boundaryEta,
                          "eta of the momentum at the calo boundary (convention of eta for pt = 0)",
                          precision_);
  table->addColumn<float>("boundaryPhi", boundaryPhi, "phi of the momentum at the calo boundary", precision_);
  table->addColumn<float>("boundaryE", boundaryE, "Energy at the calo boundary [GeV]", precision_);
  event.put(particleToken_, std::move(table));

  const size_t nRows = selectedVertices->size();
  std::vector<float> vx(nRows), vy(nRows), vz(nRows), vt(nRows);
  std::vector<uint8_t> vReason(nRows), vRole(nRows);
  std::vector<bool> vSignal(nRows);
  for (size_t r = 0; r < nRows; ++r) {
    auto const& v = g.vertices()[(*selectedVertices)[r]];
    vx[r] = v.position.x();
    vy[r] = v.position.y();
    vz[r] = v.position.z();
    vt[r] = v.position.t();
    vReason[r] = v.reason;
    vRole[r] = v.role;
    vSignal[r] = v.isSignal();
  }
  auto vertices = std::make_unique<nanoaod::FlatTable>(nRows, vertexName_, false);
  vertices->setDoc("MC-truth graph vertices: the production vertex and the endpoint of each particle of " + name_ +
                   ", and the truth vertices the vertex association maps match");
  vertices->addColumn<uint32_t>(
      "graphId", *selectedVertices, "Vertex id in the truth graph (key of the vertex association maps)");
  vertices->addColumn<float>("x", vx, "x [cm]", positionPrecision_);
  vertices->addColumn<float>("y", vy, "y [cm]", positionPrecision_);
  vertices->addColumn<float>("z", vz, "z [cm]", positionPrecision_);
  vertices->addColumn<float>("t", vt, "Time [ns]", positionPrecision_);
  vertices->addColumn<uint8_t>(
      "reason",
      vReason,
      "truth::VertexReason (1 primary, 2 decay, 3 brem, 5 conversion, 11 hadronic inelastic, 17 hard scatter, ...)");
  vertices->addColumn<uint8_t>(
      "role", vRole, "truth::VertexRole: 0 Normal, 1 InitialState, 2 UnderlyingEvent, 3 Interaction, 4 BeamSideInput");
  vertices->addColumn<bool>("isSignal", vSignal, "True for the signal interaction, false for pile-up");
  event.put(vertexToken_, std::move(vertices));

  auto children = std::make_unique<nanoaod::FlatTable>(parentIdx.size(), childName_, false);
  children->setDoc("Parent -> child links between the particles of " + name_ +
                   ", through every decay or interaction vertex of the parent, ordered by parent");
  children->addColumn<int32_t>("parentIdx", parentIdx, "Row of the parent in " + name_);
  children->addColumn<int32_t>("childIdx", childIdx, "Row of the child in " + name_);
  event.put(childToken_, std::move(children));

  event.put(selectedToken_, std::move(selected));
  event.put(selectedVerticesToken_, std::move(selectedVertices));
}

void TruthGraphFlatTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("truthLogicalGraphProducer"))->setComment("The truth::Graph");
  desc.add<std::string>("name", "TruthParticles")->setComment("Name of the particle table");
  desc.add<std::string>("doc", "MC-truth graph particles (truth::Graph)")->setComment("Doc of the particle table");
  desc.add<std::string>("vertexName", "TruthVertices")->setComment("Name of the vertex table (instance 'vertices')");
  desc.add<std::string>("childName", "TruthParticleChildren")
      ->setComment("Name of the parent -> child link table (instance 'children')");
  desc.add<bool>("keepSignal", true)->setComment("Export every particle of the signal interaction");
  desc.add<unsigned int>("keepLevelMask", 0)
      ->setComment(
          "Export every particle with any of these truth::LevelFlag bits, pile-up included (1 << 3: CaloBoundary, "
          "every particle that entered the calorimeter)");
  desc.add<std::vector<edm::InputTag>>("keepMatchedBy", {})
      ->setComment(
          "Reco -> truth ticl::TICLAssociationMap<ticl::mapWithSharedEnergyAndScore> maps: export every particle a "
          "reco object is matched to");
  desc.add<std::vector<edm::InputTag>>("keepMatchedByTruthToReco", {})
      ->setComment(
          "Truth -> reco maps of the same type (a row per graph particle): export every particle with an entry");
  desc.add<std::vector<edm::InputTag>>("keepVerticesMatchedBy", {})
      ->setComment(
          "Reco -> truth ticl::TICLAssociationMap<ticl::mapWithFractionAndScore> vertex maps: export every truth "
          "vertex a reco vertex is matched to");
  desc.add<std::vector<edm::InputTag>>("keepVerticesMatchedByTruthToReco", {})
      ->setComment("Truth -> reco vertex maps (a row per graph vertex): export every vertex with an entry");
  desc.add<bool>("keepAncestors", true)->setComment("Export the ancestors of every exported particle");
  desc.add<bool>("skipNonExistingSrc", true)
      ->setComment("Write empty tables when the graph is missing, and ignore a missing map (false: both are errors)");
  // A float has 23 mantissa bits. 14 bits are a relative precision of 6e-5, 16 bits of 1.5e-5
  // (15 um at 1 m); the zeros of absent values stay exact.
  desc.add<int>("precision", 14)
      ->setComment("Mantissa bits kept for pt, eta, phi and energy, at production and at the boundary (-1: all)");
  desc.add<int>("positionPrecision", 16)->setComment("Mantissa bits kept for positions and times (-1: all)");
  descriptions.add("truthGraphFlatTable", desc);
}

class TruthGraphHitFlatTableProducer : public edm::global::EDProducer<> {
public:
  explicit TruthGraphHitFlatTableProducer(edm::ParameterSet const& cfg)
      : graphToken_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        hitIndexToken_(consumes(cfg.getParameter<edm::InputTag>("hitIndex"))),
        particlesToken_(consumes(cfg.getParameter<edm::InputTag>("particles"))),
        hitTableToken_(consumes(cfg.getParameter<edm::InputTag>("hitTable"))),
        name_(cfg.getParameter<std::string>("name")),
        doc_(cfg.getParameter<std::string>("doc")),
        hitIdColumn_(cfg.getParameter<std::string>("hitIdColumn")),
        channel_(cfg.getParameter<std::string>("channel") == "Muon" ? truth::HitChannel::Muon
                                                                    : truth::HitChannel::Calo),
        levelMask_(cfg.getParameter<unsigned int>("levelMask")),
        skipNonExistingSrc_(cfg.getParameter<bool>("skipNonExistingSrc")),
        precision_(cfg.getParameter<int>("precision")),
        token_(produces<nanoaod::FlatTable>()) {}

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  const edm::EDGetTokenT<truth::Graph> graphToken_;
  const edm::EDGetTokenT<truth::LogicalGraphHitIndex> hitIndexToken_;
  const edm::EDGetTokenT<std::vector<unsigned int>> particlesToken_;
  const edm::EDGetTokenT<nanoaod::FlatTable> hitTableToken_;
  const std::string name_, doc_, hitIdColumn_;
  const truth::HitChannel channel_;
  const uint32_t levelMask_;
  const bool skipNonExistingSrc_;
  const int precision_;
  const edm::EDPutTokenT<nanoaod::FlatTable> token_;
};

void TruthGraphHitFlatTableProducer::produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const {
  std::vector<int32_t> particleIdx, hitIdx;
  std::vector<float> energy;
  auto const graph = event.getHandle(graphToken_);
  auto const hitIndex = event.getHandle(hitIndexToken_);
  auto const particles = event.getHandle(particlesToken_);
  auto const hits = event.getHandle(hitTableToken_);
  const bool present = graph.isValid() && hitIndex.isValid() && particles.isValid() && hits.isValid();
  if (!present && !skipNonExistingSrc_)
    throw cms::Exception("ProductNotFound") << "TruthGraphHitFlatTableProducer '" << name_
                                            << "': the truth graph, the hit index, the particle ids or the hit table "
                                               "is missing";
  if (present) {
    const int column = hits->columnIndex(hitIdColumn_);
    if (column < 0 || hits->columnType(column) != nanoaod::FlatTable::ColumnType::UInt32)
      throw cms::Exception("Configuration") << "TruthGraphHitFlatTableProducer '" << name_ << "': hit table '"
                                            << hits->name() << "' has no uint32 column '" << hitIdColumn_ << "'";
    // (DetId, row) of every hit, sorted once; each lookup is a binary search. Several rows
    // can share a DetId: on the Muon channel the simulation's unit (a CSC layer, an RPC roll,
    // a GEM partition) holds several rechits, and each gets the unit's simulated energy.
    std::vector<std::pair<uint32_t, int32_t>> hitRows;
    hitRows.reserve(hits->nRows());
    for (int32_t row = 0; const uint32_t detId : hits->columnData<uint32_t>(column))
      hitRows.emplace_back(detId, row++);
    std::sort(hitRows.begin(), hitRows.end());
    auto const byDetId = [](std::pair<uint32_t, int32_t> const& a, std::pair<uint32_t, int32_t> const& b) {
      return a.first < b.first;
    };

    std::vector<truth::LogicalGraphHitIndex::Hit> subgraph;
    std::vector<std::pair<int32_t, float>> found;  // (hit row, sim energy), summed per row below
    for (uint32_t r = 0; r < particles->size(); ++r) {
      const uint32_t id = (*particles)[r];
      if (id >= graph->nParticles() || !(graph->particles()[id].levelFlags & levelMask_))
        continue;
      subgraph.clear();
      hitIndex->appendSubgraphHits(channel_, id, subgraph);
      found.clear();
      for (auto const& hit : subgraph) {
        const auto [first, last] =
            std::equal_range(hitRows.begin(), hitRows.end(), std::pair(hit.detId, int32_t{0}), byDetId);
        for (auto it = first; it != last; ++it)
          found.emplace_back(it->second, hit.energy);
      }
      std::sort(found.begin(), found.end());
      for (size_t k = 0; k < found.size(); ++k) {
        if (k > 0 && found[k].first == found[k - 1].first) {
          energy.back() += found[k].second;
          continue;
        }
        particleIdx.push_back(r);
        hitIdx.push_back(found[k].first);
        energy.push_back(found[k].second);
      }
    }
  }
  auto table = std::make_unique<nanoaod::FlatTable>(particleIdx.size(), name_, false);
  table->setDoc(doc_);
  table->addColumn<int32_t>("particleIdx", particleIdx, "Row of the truth particle in its table");
  table->addColumn<int32_t>("hitIdx", hitIdx, "Row of the hit in its table");
  table->addColumn<float>("simEnergy",
                          energy,
                          "Simulated energy the particle and its descendants left in the hit's cell [GeV]; on the "
                          "Muon channel the energy of the hit's unit, given to every hit of that unit (do not sum "
                          "over hits sharing a DetId). Energy in cells without a hit-table row is not written",
                          precision_);
  event.put(token_, std::move(table));
}

void TruthGraphHitFlatTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("truthLogicalGraphProducer"))->setComment("The truth::Graph");
  desc.add<edm::InputTag>("hitIndex", edm::InputTag("truthLogicalGraphHitIndexProducer"))
      ->setComment("The truth::LogicalGraphHitIndex of the graph");
  desc.add<edm::InputTag>("particles", edm::InputTag("truthGraphTable", "selected"))
      ->setComment(
          "Graph ids of the rows of the particle table, in row order (the 'selected' product of "
          "TruthGraphFlatTableProducer)");
  desc.add<edm::InputTag>("hitTable")->setComment("A hit table with a uint32 DetId column");
  desc.add<std::string>("hitIdColumn", "detId")->setComment("The DetId column of the hit table");
  desc.ifValue(edm::ParameterDescription<std::string>("channel", "Calo", true),
               edm::allowedValues<std::string>("Calo", "Muon"))
      ->setComment(
          "Channel of the hit index: Calo (HGCAL, ECAL, HCAL) or Muon (keyed by DT wire, CSC layer, RPC roll, "
          "GEM eta partition)");
  desc.add<unsigned int>("levelMask", 1u << 3)
      ->setComment(
          "Exported particles with any of these truth::LevelFlag bits get their subgraph hits (default "
          "CaloBoundary: an antichain, so no hit is counted twice)");
  desc.add<std::string>("name", "TruthHits")->setComment("Name of the table");
  desc.add<std::string>("doc", "Simulated energy of truth-graph particles in hits")->setComment("Doc of the table");
  desc.add<bool>("skipNonExistingSrc", true)->setComment("Write an empty table when an input is missing");
  desc.add<int>("precision", 12)
      ->setComment("Mantissa bits kept for the simulated energy (-1: all; 12 bits: a relative precision of 2.4e-4)");
  descriptions.add("truthGraphHitFlatTable", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TruthGraphFlatTableProducer);
DEFINE_FWK_MODULE(TruthGraphHitFlatTableProducer);
