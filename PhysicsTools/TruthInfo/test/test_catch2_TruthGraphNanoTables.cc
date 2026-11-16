// Unit tests of the truth-graph NanoAOD tables (plugins/TruthGraphNanoTableProducers.cc) and of
// their use as the key space of the branch association tables (plugins/TruthBranchNanoTableProducers.cc).

#include "catch2/catch_all.hpp"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/TestProcessor/interface/TestProcessor.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "PhysicsTools/TruthInfo/interface/LogicalGraphHitIndexBuilder.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/TruthInfo/interface/LogicalGraphHitIndex.h"

namespace {

  using EnergyMap = ticl::TICLAssociationMap<ticl::mapWithSharedEnergyAndScore>;
  using FractionMap = ticl::TICLAssociationMap<ticl::mapWithFractionAndScore>;
  constexpr uint32_t kCaloBoundary = 1u << 3;
  constexpr uint32_t kCHadrons = 1u << 9;
  const uint64_t kPileup = EncodedEventId(0, 3).rawId();
  constexpr uint32_t kVertices = 9;

  // The CSR graph builder of the other tests in this package (LevelFlags_t.cpp).
  struct GraphBuilder {
    GraphBuilder(uint32_t nParticles, uint32_t nVertices) {
      graph.particles().resize(nParticles);
      graph.vertices().resize(nVertices);
    }
    void addDecay(uint32_t particleId, uint32_t vertexId) {
      d2v.emplace_back(particleId, vertexId);
      v2i.emplace_back(vertexId, particleId);
    }
    void addProduction(uint32_t vertexId, uint32_t particleId) {
      v2o.emplace_back(vertexId, particleId);
      p2v.emplace_back(particleId, vertexId);
    }
    static void csr(uint32_t n,
                    std::vector<std::pair<uint32_t, uint32_t>>& pairs,
                    std::vector<uint32_t>& off,
                    std::vector<uint32_t>& flat) {
      std::sort(pairs.begin(), pairs.end());
      off.assign(n + 1, 0);
      for (auto const& pr : pairs)
        ++off[pr.first + 1];
      for (uint32_t i = 1; i <= n; ++i)
        off[i] += off[i - 1];
      flat.assign(pairs.size(), 0);
      auto cur = off;
      for (auto const& pr : pairs)
        flat[cur[pr.first]++] = pr.second;
    }
    truth::Graph finish() {
      csr(graph.nParticles(), d2v, graph.particleToDecayVertexOffsets(), graph.particleToDecayVertices());
      csr(graph.nParticles(), p2v, graph.particleToProductionVertexOffsets(), graph.particleToProductionVertices());
      csr(graph.nVertices(), v2o, graph.vertexToOutgoingParticleOffsets(), graph.vertexToOutgoingParticles());
      csr(graph.nVertices(), v2i, graph.vertexToIncomingParticleOffsets(), graph.vertexToIncomingParticles());
      REQUIRE(graph.isConsistent());
      return graph;
    }
    truth::Graph graph;
    std::vector<std::pair<uint32_t, uint32_t>> d2v, p2v, v2o, v2i;
  };

  // Signal: beam proton p9 -> v0 -> top p0 -> v1 -> {pi+ p1 (CaloBoundary), nu p2 (pt 0)};
  //         p1 -> v2 -> e p3 (SIM only); p1 also ends at v6 and v7 (same time, later than v2),
  //         so its endpoint is v7; p3 has a second production vertex v8 (unused: the first wins).
  // Pile-up: v3 -> {p4, p6, p7}; p4 (CHadrons bit) -> v4 -> pi- p5 (CaloBoundary); v5 -> p8.
  // Vertex v has position (v, 0, v + 1) and time v, except v7 (time 6, as v6).
  // Kept with keepLevelMask CaloBoundary and a map matching p7: p0-p5, p7, p9 (p4 as the
  // ancestor of p5); p6 and p8 are dropped. Vertices kept: v0-v4 and v7.
  truth::Graph buildGraph() {
    GraphBuilder b(10, kVertices);
    auto particle = [&b](uint32_t id, int pdgId, uint64_t eventId, bool gen, bool sim, math::XYZTLorentzVectorD p) {
      auto& d = b.graph.particles()[id];
      d.pdgId = pdgId;
      d.eventId = eventId;
      d.genNode = gen ? 100 + id : -1;
      d.simNode = sim ? 200 + id : -1;
      d.momentum = p;
    };
    particle(0, 6, 0, true, false, {10., 0., 50., 180.});
    particle(1, 211, 0, true, true, {20., 0., 30., 40.});
    particle(2, 12, 0, true, false, {0., 0., 5., 5.});
    particle(3, 11, 0, false, true, {1., 1., 1., 2.});
    particle(4, 2212, kPileup, true, false, {1., 0., 1., 2.});
    particle(5, -211, kPileup, false, true, {1., 0., -1., 2.});
    particle(6, 22, kPileup, true, false, {1., 0., 0., 1.});
    particle(7, 22, kPileup, true, false, {0., 1., 0., 1.});
    particle(8, 22, kPileup, true, false, {0., 0., 1., 1.});
    particle(9, 2212, 0, true, false, {0., 0., 7000., 7000.});
    b.graph.particles()[0].statusFlags = 0x101;
    b.graph.particles()[1].levelFlags = kCaloBoundary;
    b.graph.particles()[1].checkpoints.push_back(truth::Checkpoint{0, {120., 0., 200., 5.}, {20., 0., 30., 40.}});
    b.graph.particles()[4].levelFlags = kCHadrons;
    b.graph.particles()[5].levelFlags = kCaloBoundary;
    b.graph.particles()[5].checkpoints.push_back(truth::Checkpoint{0, {120., 10., -200., 4.}, {1., 0., -1., 2.}});
    for (uint32_t v = 0; v < kVertices; ++v) {
      b.graph.vertices()[v].eventId = (v >= 3 && v <= 5) ? kPileup : 0;
      b.graph.vertices()[v].position = math::XYZTLorentzVectorD(v, 0., v + 1., v == 7 ? 6. : v);
    }
    b.graph.vertices()[0].role = 1;    // truth::VertexRole::InitialState
    b.graph.vertices()[1].reason = 2;  // truth::VertexReason::Decay
    b.addDecay(9, 0);
    b.addProduction(0, 0);
    b.addDecay(0, 1);
    b.addProduction(1, 1);
    b.addProduction(1, 2);
    b.addDecay(1, 2);
    b.addDecay(1, 6);
    b.addDecay(1, 7);
    b.addProduction(2, 3);
    b.addProduction(8, 3);
    b.addProduction(3, 4);
    b.addProduction(3, 6);
    b.addProduction(3, 7);
    b.addDecay(4, 4);
    b.addProduction(4, 5);
    b.addProduction(5, 8);
    return b.finish();
  }

  // Reco object 0 matched to p7, reco object 1 to p1.
  std::unique_ptr<EnergyMap> recoToTruthMap() {
    auto map = std::make_unique<EnergyMap>(2);
    map->insert(0, 7, 3.f, 0.2f);
    map->insert(1, 1, 5.f, 0.1f);
    return map;
  }

  // A truth -> reco map, a row per graph particle: p6 shares with reco object 0.
  std::unique_ptr<EnergyMap> truthToRecoMap() {
    auto map = std::make_unique<EnergyMap>(10);
    map->insert(6, 0, 0.5f, 0.5f);
    return map;
  }

  std::string graphTable(std::string const& extra = "") {
    return R"_(
process.truthGraphTable = cms.EDProducer("TruthGraphFlatTableProducer",
    src = cms.InputTag("graph"), keepLevelMask = cms.uint32(1 << 3), keepMatchedBy = cms.VInputTag("matches")
    )_" + extra +
           ")\n";
  }

  const std::string kHeader = "from FWCore.TestProcessor.TestProcess import *\nprocess = TestProcess()\n";

  template <typename T>
  std::vector<T> column(nanoaod::FlatTable const& table, std::string const& name) {
    const int c = table.columnIndex(name);
    REQUIRE(c >= 0);
    auto const data = table.columnData<T>(c);
    return std::vector<T>(data.begin(), data.end());
  }

  std::string columnDoc(nanoaod::FlatTable const& table, std::string const& name) {
    const int c = table.columnIndex(name);
    REQUIRE(c >= 0);
    return table.columnDoc(c);
  }

  std::vector<std::pair<std::string, nanoaod::FlatTable::ColumnType>> layout(nanoaod::FlatTable const& table) {
    std::vector<std::pair<std::string, nanoaod::FlatTable::ColumnType>> out;
    for (unsigned int c = 0; c < table.nColumns(); ++c)
      out.emplace_back(table.columnName(c), table.columnType(c));
    return out;
  }

  // The graph ids of the rows of TruthGraphFlatTableProducer configured with `extra`.
  std::vector<unsigned int> selection(std::string const& extra, bool withTruthToReco = false) {
    edm::test::TestProcessor::Config config{kHeader + graphTable(extra) +
                                            "process.moduleToTest(process.truthGraphTable)\n"};
    auto graphToken = config.produces<truth::Graph>("graph");
    auto mapToken = config.produces<EnergyMap>("matches");
    auto truthToRecoToken = config.produces<EnergyMap>("truthToReco");
    edm::test::TestProcessor tester(config);
    auto event = withTruthToReco ? tester.test(std::make_pair(graphToken, std::make_unique<truth::Graph>(buildGraph())),
                                               std::make_pair(mapToken, recoToTruthMap()),
                                               std::make_pair(truthToRecoToken, truthToRecoMap()))
                                 : tester.test(std::make_pair(graphToken, std::make_unique<truth::Graph>(buildGraph())),
                                               std::make_pair(mapToken, recoToTruthMap()));
    return *event.get<std::vector<unsigned int>>("selected");
  }

}  // namespace

TEST_CASE("TruthGraphFlatTableProducer selects, links and keeps its columns", "[TruthGraphFlatTableProducer]") {
  edm::test::TestProcessor::Config config{kHeader + graphTable() + "process.moduleToTest(process.truthGraphTable)\n"};
  auto graphToken = config.produces<truth::Graph>("graph");
  auto mapToken = config.produces<EnergyMap>("matches");
  edm::test::TestProcessor tester(config);

  auto event = tester.test(std::make_pair(graphToken, std::make_unique<truth::Graph>(buildGraph())),
                           std::make_pair(mapToken, recoToTruthMap()));
  auto const& particles = *event.get<nanoaod::FlatTable>();
  auto const& vertices = *event.get<nanoaod::FlatTable>("vertices");
  auto const& children = *event.get<nanoaod::FlatTable>("children");

  SECTION("selection: signal, CaloBoundary, matched, ancestors; rows in ascending graph id") {
    const std::vector<unsigned int> expected{0, 1, 2, 3, 4, 5, 7, 9};
    REQUIRE(*event.get<std::vector<unsigned int>>("selected") == expected);
    REQUIRE(column<uint32_t>(particles, "graphId") == std::vector<uint32_t>(expected.begin(), expected.end()));
  }
  SECTION("vertices: the production vertex and the endpoint of each exported particle, ascending") {
    REQUIRE(*event.get<std::vector<unsigned int>>("selectedVertices") == std::vector<unsigned int>{0, 1, 2, 3, 4, 7});
    REQUIRE(column<uint8_t>(vertices, "reason")[1] == 2);  // v1: decay
    REQUIRE(column<uint8_t>(vertices, "role")[0] == 1);    // v0: initial state
    REQUIRE(column<bool>(vertices, "isSignal") == std::vector<bool>{true, true, true, false, false, true});
    REQUIRE(column<uint32_t>(vertices, "graphId") == std::vector<uint32_t>{0, 1, 2, 3, 4, 7});
    // rows: p0 p1 p2 p3 p4 p5 p7 p9. p3: first production vertex (v2, not v8).
    REQUIRE(column<int32_t>(particles, "prodVertexIdx") == std::vector<int32_t>{0, 1, 1, 2, 3, 4, 3, -1});
    // p1 ends at v2 (t 2), v6 and v7 (t 6): the latest, the higher id on the tie, is v7 (row 5).
    REQUIRE(column<int32_t>(particles, "decayVertexIdx") == std::vector<int32_t>{1, 5, -1, -1, 4, -1, -1, 0});
  }
  SECTION("parent -> child links through every end vertex") {
    REQUIRE(column<int32_t>(children, "parentIdx") == std::vector<int32_t>{0, 0, 1, 4, 7});
    REQUIRE(column<int32_t>(children, "childIdx") == std::vector<int32_t>{1, 2, 3, 5, 0});
  }
  SECTION("column layout") {
    using enum nanoaod::FlatTable::ColumnType;
    auto names = [](nanoaod::FlatTable const& table) {
      std::vector<std::string> out;
      for (auto const& [name, type] : layout(table))
        out.push_back(name);
      return out;
    };
    REQUIRE(names(particles) ==
            std::vector<std::string>{"graphId",     "pdgId",         "charge",         "pt",          "eta",
                                     "phi",         "energy",        "status",         "statusFlags", "levelFlags",
                                     "isSignal",    "bunchCrossing", "eventIndex",     "hasGen",      "hasSim",
                                     "role",        "prodVertexIdx", "decayVertexIdx", "hasBoundary", "boundaryX",
                                     "boundaryY",   "boundaryZ",     "boundaryT",      "boundaryPt",  "boundaryEta",
                                     "boundaryPhi", "boundaryE"});
    REQUIRE(layout(vertices) == decltype(layout(vertices)){{"graphId", UInt32},
                                                           {"x", Float},
                                                           {"y", Float},
                                                           {"z", Float},
                                                           {"t", Float},
                                                           {"reason", UInt8},
                                                           {"role", UInt8},
                                                           {"isSignal", Bool}});
    REQUIRE(layout(children) == decltype(layout(children)){{"parentIdx", Int32}, {"childIdx", Int32}});
  }
  SECTION("boundary from checkpoint 0, zeros and hasBoundary false elsewhere; production position") {
    auto const has = column<bool>(particles, "hasBoundary");
    auto const bx = column<float>(particles, "boundaryX");
    for (size_t r = 0; r < has.size(); ++r) {
      REQUIRE(has[r] == (r == 1 || r == 5));
      if (!has[r])
        REQUIRE(bx[r] == 0.f);
    }
    REQUIRE(bx[1] == 120.f);
    REQUIRE(column<float>(particles, "boundaryZ")[5] == -200.f);
    REQUIRE(column<float>(particles, "boundaryT")[1] == 5.f);
    REQUIRE(column<float>(particles, "boundaryE")[5] == 2.f);
    REQUIRE(column<float>(particles, "boundaryPt")[1] == 20.f);
    REQUIRE(column<float>(particles, "boundaryPhi")[1] == 0.f);
    REQUIRE(column<float>(particles, "boundaryEta")[5] == Catch::Approx(-0.8814).margin(1e-3));  // asinh(-1)
    // p3 is produced at v2 (row 2), at z = 3; the production position is a TruthVertices row.
    REQUIRE(column<float>(vertices, "z")[column<int32_t>(particles, "prodVertexIdx")[3]] == 3.f);
  }
  SECTION("kinematics and flags") {
    REQUIRE(column<float>(particles, "eta")[2] == Catch::Approx(22761.f));     // pt = 0: ROOT's +(22756 + pz)
    REQUIRE(column<float>(particles, "charge")[0] == Catch::Approx(2. / 3.));  // top
    REQUIRE(column<float>(particles, "charge")[5] == -1.f);                    // pi-
    REQUIRE(column<float>(particles, "pt")[1] == 20.f);
    REQUIRE(column<float>(particles, "energy")[1] == 40.f);
    REQUIRE(column<float>(particles, "phi")[7] == 0.f);
    REQUIRE(column<bool>(particles, "isSignal") ==
            std::vector<bool>{true, true, true, true, false, false, false, true});
    REQUIRE(column<int32_t>(particles, "eventIndex")[4] == 3);
    REQUIRE(column<uint16_t>(particles, "statusFlags")[0] == 0x101);
    REQUIRE(column<uint32_t>(particles, "levelFlags")[5] == kCaloBoundary);
  }

  SECTION("a missing graph, and a present but empty one, give empty tables with the same columns") {
    // Taken before the next event, which reuses the products of this one.
    const std::vector<decltype(layout(particles))> filled{layout(particles), layout(vertices), layout(children)};
    for (bool present : {false, true}) {
      auto empty = present ? tester.test(std::make_pair(graphToken, std::make_unique<truth::Graph>())) : tester.test();
      for (size_t k = 0; std::string const instance : {"", "vertices", "children"}) {
        auto const& none = *empty.get<nanoaod::FlatTable>(instance);
        REQUIRE(none.nRows() == 0);
        REQUIRE(layout(none) == filled[k++]);
      }
      REQUIRE(empty.get<std::vector<unsigned int>>("selected")->empty());
    }
  }
}

TEST_CASE("TruthGraphFlatTableProducer selection switches", "[TruthGraphFlatTableProducer]") {
  SECTION("keepSignal false: CaloBoundary, matched and their ancestors") {
    REQUIRE(selection(", keepSignal = cms.bool(False)") == std::vector<unsigned int>{0, 1, 4, 5, 7, 9});
  }
  SECTION("keepAncestors false: no p4 above p5") {
    REQUIRE(selection(", keepAncestors = cms.bool(False)") == std::vector<unsigned int>{0, 1, 2, 3, 5, 7, 9});
  }
  SECTION("both false: CaloBoundary and matched only") {
    REQUIRE(selection(", keepSignal = cms.bool(False), keepAncestors = cms.bool(False)") ==
            std::vector<unsigned int>{1, 5, 7});
  }
  SECTION("keepMatchedByTruthToReco: the truth rows with an entry are kept") {
    REQUIRE(selection(", keepMatchedByTruthToReco = cms.VInputTag('truthToReco')", true) ==
            std::vector<unsigned int>{0, 1, 2, 3, 4, 5, 6, 7, 9});
  }
  SECTION("skipNonExistingSrc false: a missing graph is an error") {
    edm::test::TestProcessor::Config config{kHeader + graphTable(", skipNonExistingSrc = cms.bool(False)") +
                                            "process.moduleToTest(process.truthGraphTable)\n"};
    auto mapToken = config.produces<EnergyMap>("matches");
    config.produces<truth::Graph>("graph");
    edm::test::TestProcessor tester(config);
    REQUIRE_THROWS_AS(tester.test(std::make_pair(mapToken, recoToTruthMap())), cms::Exception);
  }
}

TEST_CASE("Branch association tables index the truth-graph tables", "[TruthBranchAssociationFlatTableProducer]") {
  const std::string particleMap = kHeader + graphTable() + R"_(
process.toTest = cms.EDProducer("TruthBranchAssociationFlatTableProducer", src = cms.InputTag("matches"),
    targets = cms.InputTag("truthGraphTable", "selected"), name = cms.string("RecoToTruth")EXTRA)
process.moduleToTest(process.toTest, cms.Task(process.truthGraphTable))
)_";
  auto runParticleMap = [&particleMap](std::string const& extra) -> nanoaod::FlatTable {
    std::string cfg = particleMap;
    cfg.replace(cfg.find("EXTRA"), 5, extra);
    edm::test::TestProcessor::Config config{cfg};
    auto graphToken = config.produces<truth::Graph>("graph");
    auto mapToken = config.produces<EnergyMap>("matches");
    edm::test::TestProcessor tester(config);
    return *tester
                .test(std::make_pair(graphToken, std::make_unique<truth::Graph>(buildGraph())),
                      std::make_pair(mapToken, recoToTruthMap()))
                .get<nanoaod::FlatTable>("RecoToTruthLinksTable");
  };
  SECTION("particle map, reco -> truth: links are rows of TruthParticles; default payload columns") {
    auto const links = runParticleMap("");
    REQUIRE(column<uint32_t>(links, "index") == std::vector<uint32_t>{6, 1});  // p7 is row 6, p1 row 1
    REQUIRE(column<uint32_t>(links, "graphId") == std::vector<uint32_t>{7, 1});
    REQUIRE(column<float>(links, "sharedEnergy") == std::vector<float>{3.f, 5.f});
    REQUIRE(columnDoc(links, "sharedEnergy") == "Energy the reco object shares with the truth branch [GeV]");
    REQUIRE(columnDoc(links, "index") == "Row in the truth particle table");
  }
  SECTION("track payload: a shared hit count, with its own name and docs") {
    auto const links = runParticleMap(
        ", recoToTruthPayload = cms.PSet(name = cms.string('nSharedHits'), doc = cms.string('Number of hits'), "
        "scoreDoc = cms.string('1 - hit fraction'))");
    REQUIRE(links.columnIndex("sharedEnergy") < 0);
    REQUIRE(column<float>(links, "nSharedHits") == std::vector<float>{3.f, 5.f});
    REQUIRE(columnDoc(links, "nSharedHits") == "Number of hits");
    REQUIRE(columnDoc(links, "score") == "1 - hit fraction");
  }
  SECTION("vertex map, both directions: links and rows are rows of TruthVertices; matched vertices kept") {
    // Without keepVerticesMatchedBy a link to a vertex that is not exported is left out; with
    // it, the matched vertex v5 gets a row and the link is kept.
    for (bool keepMatched : {false, true}) {
      for (bool recoToTruth : {true, false}) {
        const std::string keep =
            !keepMatched ? ""
                         : (recoToTruth ? ", keepVerticesMatchedBy = cms.VInputTag('vertexMap')"
                                        : ", keepVerticesMatchedByTruthToReco = cms.VInputTag('vertexMap')");
        const std::string cfg = kHeader + graphTable(keep) + R"_(
process.toTest = cms.EDProducer("TruthBranchVertexAssociationFlatTableProducer", src = cms.InputTag("vertexMap"),
    targets = cms.InputTag("truthGraphTable", "selectedVertices"), name = cms.string("VertexAssociation"),
    direction = cms.string(")_" +
                                std::string(recoToTruth ? "recoToTruth" : "truthToReco") + R"_("))
process.moduleToTest(process.toTest, cms.Task(process.truthGraphTable))
)_";
        edm::test::TestProcessor::Config config{cfg};
        auto graphToken = config.produces<truth::Graph>("graph");
        auto mapToken = config.produces<EnergyMap>("matches");
        auto vertexToken = config.produces<FractionMap>("vertexMap");
        edm::test::TestProcessor tester(config);
        // recoToTruth: one reco vertex sharing with v3 (exported, row 3) and v5 (not a production
        // vertex or endpoint of an exported particle). truthToReco: a row per graph vertex;
        // v3 -> reco 0, v5 -> reco 0.
        auto map = std::make_unique<FractionMap>(recoToTruth ? 1 : kVertices);
        if (recoToTruth) {
          map->insert(0, 3, 0.75f, 0.25f);
          map->insert(0, 5, 0.25f, 0.75f);
        } else {
          map->insert(3, 0, 0.5f, 0.5f);
          map->insert(5, 0, 1.f, 0.f);
        }
        auto event = tester.test(std::make_pair(graphToken, std::make_unique<truth::Graph>(buildGraph())),
                                 std::make_pair(mapToken, recoToTruthMap()),
                                 std::make_pair(vertexToken, std::move(map)));
        auto const& main = *event.get<nanoaod::FlatTable>();
        auto const& links = *event.get<nanoaod::FlatTable>("VertexAssociationLinksTable");
        // with the keep, v5 is row 5 (rows: v0 v1 v2 v3 v4 v5 v7)
        REQUIRE(column<uint32_t>(links, "graphId") ==
                (keepMatched ? std::vector<uint32_t>{3, 5} : std::vector<uint32_t>{3}));
        if (recoToTruth) {
          REQUIRE(main.nRows() == 1);
          REQUIRE(column<float>(links, "fraction") ==
                  (keepMatched ? std::vector<float>{0.75f, 0.25f} : std::vector<float>{0.75f}));
          REQUIRE(column<uint32_t>(links, "index") ==
                  (keepMatched ? std::vector<uint32_t>{3, 5} : std::vector<uint32_t>{3}));
          REQUIRE(columnDoc(links, "index") == "Row in the truth vertex table");
        } else {
          REQUIRE(column<int32_t>(main, "nVertexAssociationLinks") ==
                  (keepMatched ? std::vector<int32_t>{0, 0, 0, 1, 0, 1, 0} : std::vector<int32_t>{0, 0, 0, 1, 0, 0}));
          REQUIRE(column<float>(links, "fraction") ==
                  (keepMatched ? std::vector<float>{0.5f, 1.f} : std::vector<float>{0.5f}));
          REQUIRE(column<uint32_t>(links, "index") ==
                  (keepMatched ? std::vector<uint32_t>{0, 0} : std::vector<uint32_t>{0}));
        }
      }
    }
  }
}

TEST_CASE("TruthGraphHitFlatTableProducer sums subgraph sim energy per hit row", "[TruthGraphHitFlatTableProducer]") {
  const std::string python = kHeader + graphTable() + R"_(
process.toTest = cms.EDProducer("TruthGraphHitFlatTableProducer", src = cms.InputTag("graph"),
    hitIndex = cms.InputTag("hitIndex"), hitTable = cms.InputTag("hits"), hitIdColumn = cms.string("COLUMN"),
    levelMask = cms.uint32(MASK))
process.moduleToTest(process.toTest, cms.Task(process.truthGraphTable))
)_";
  // p1 (CaloBoundary) and its child p3 share cell 1000; p3 also hits 2000; p5 (CaloBoundary)
  // hits 3000; p7 (exported, not CaloBoundary) hits 1000; p6 (not exported) hits 2000. The
  // GEN-only p4 is the hit-index parent of p5 and p7, which the SIM root p6 separates in the
  // DFS order, so its subgraph is two ranges.
  auto buildIndex = []() {
    truth::LogicalGraphHitIndexBuilder builder(10);
    for (uint32_t id : {1u, 3u})
      builder.setSimTrackForParticle(id, 0, 100 + id);
    for (uint32_t id : {5u, 6u, 7u})
      builder.setSimTrackForParticle(id, kPileup, 100 + id);
    builder.addParticleChild(1, 3);
    builder.addParticleChild(4, 5);
    builder.addParticleChild(4, 7);
    builder.addHit(truth::HitChannel::Calo, 0, 101, 1000, 1.f, 0);
    builder.addHit(truth::HitChannel::Calo, 0, 103, 1000, 0.5f, 0);
    builder.addHit(truth::HitChannel::Calo, 0, 103, 2000, 2.f, 1);
    builder.addHit(truth::HitChannel::Calo, kPileup, 105, 3000, 4.f, 2);
    builder.addHit(truth::HitChannel::Calo, kPileup, 106, 2000, 7.f, 1);
    builder.addHit(truth::HitChannel::Calo, kPileup, 107, 1000, 9.f, 0);
    return builder.finish();
  };
  // The table, copied out: an edm::test::Event does not outlive its TestProcessor.
  auto run = [&python, &buildIndex](std::string const& hitIdColumn, uint32_t mask, bool withIndex) {
    std::string cfg = python;
    cfg.replace(cfg.find("COLUMN"), 6, hitIdColumn);
    cfg.replace(cfg.find("MASK"), 4, std::to_string(mask));
    edm::test::TestProcessor::Config config{cfg};
    auto graphToken = config.produces<truth::Graph>("graph");
    auto mapToken = config.produces<EnergyMap>("matches");
    auto indexToken = config.produces<truth::LogicalGraphHitIndex>("hitIndex");
    auto hitsToken = config.produces<nanoaod::FlatTable>("hits");
    edm::test::TestProcessor tester(config);
    // Two rows hold DetId 1000: both get its energy.
    auto hits = std::make_unique<nanoaod::FlatTable>(5, "Hits", false);
    hits->addColumn<uint32_t>("detId", std::vector<uint32_t>{3000, 1000, 2000, 4000, 1000}, "");
    hits->addColumn<float>("energy", std::vector<float>(5, 1.f), "");
    if (withIndex)
      return nanoaod::FlatTable(
          *tester
               .test(std::make_pair(graphToken, std::make_unique<truth::Graph>(buildGraph())),
                     std::make_pair(mapToken, recoToTruthMap()),
                     std::make_pair(indexToken, std::make_unique<truth::LogicalGraphHitIndex>(buildIndex())),
                     std::make_pair(hitsToken, std::move(hits)))
               .get<nanoaod::FlatTable>());
    return nanoaod::FlatTable(*tester
                                   .test(std::make_pair(graphToken, std::make_unique<truth::Graph>(buildGraph())),
                                         std::make_pair(mapToken, recoToTruthMap()),
                                         std::make_pair(hitsToken, std::move(hits)))
                                   .get<nanoaod::FlatTable>());
  };

  SECTION("subgraph energy per (particle row, hit row), CaloBoundary particles") {
    auto const table = run("detId", kCaloBoundary, true);
    REQUIRE(column<int32_t>(table, "particleIdx") == std::vector<int32_t>{1, 1, 1, 5});
    REQUIRE(column<int32_t>(table, "hitIdx") == std::vector<int32_t>{1, 2, 4, 0});
    REQUIRE(column<float>(table, "simEnergy") == std::vector<float>{1.5f, 2.f, 1.5f, 4.f});
  }
  SECTION("a GEN-only particle whose subgraph spans several ranges") {
    REQUIRE(buildIndex().subgraphRanges(4).size() == 2);
    auto const table = run("detId", kCaloBoundary | kCHadrons, true);
    REQUIRE(column<int32_t>(table, "particleIdx") == std::vector<int32_t>{1, 1, 1, 4, 4, 4, 5});
    REQUIRE(column<int32_t>(table, "hitIdx") == std::vector<int32_t>{1, 2, 4, 0, 1, 4, 0});
    REQUIRE(column<float>(table, "simEnergy") == std::vector<float>{1.5f, 2.f, 1.5f, 4.f, 9.f, 9.f, 4.f});
  }
  SECTION("a missing hit index gives an empty table with the same columns") {
    auto const filled = run("detId", kCaloBoundary, true);
    auto const empty = run("detId", kCaloBoundary, false);
    REQUIRE(empty.nRows() == 0);
    REQUIRE(layout(empty) == layout(filled));
  }
  SECTION("a hit table without the configured uint32 column is a configuration error") {
    REQUIRE_THROWS_AS(run("energy", kCaloBoundary, true), cms::Exception);
  }
}
