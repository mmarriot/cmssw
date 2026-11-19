#include "catch2/catch_all.hpp"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/TestProcessor/interface/TestProcessor.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace {
  std::string config(std::string const& extra = "") {
    return R"(from FWCore.TestProcessor.TestProcess import *
from PhysicsTools.NanoAOD.vertexTrackTableProducer_cfi import vertexTrackTableProducer
process = TestProcess()
process.toTest = vertexTrackTableProducer.clone(src = "vertices", tracks = "tracks", name = "PVTrack")
)" + extra +
           "process.moduleToTest(process.toTest)\n";
  }

  reco::Vertex vertex(double z) { return reco::Vertex(reco::Vertex::Point(0, 0, z), reco::Vertex::Error(), 1., 1., 0); }
}  // namespace

TEST_CASE("VertexTrackTableProducer", "[VertexTrackTableProducer]") {
  SECTION("no vertices") {
    edm::test::TestProcessor::Config cfg{config()};
    auto vertices = cfg.produces<reco::VertexCollection>("vertices");
    auto tracks = cfg.produces<reco::TrackCollection>("tracks");
    edm::test::TestProcessor tester(cfg);
    auto event = tester.test(std::make_pair(vertices, std::make_unique<reco::VertexCollection>()),
                             std::make_pair(tracks, std::make_unique<reco::TrackCollection>()));
    auto const& table = *event.get<nanoaod::FlatTable>();
    REQUIRE(table.name() == "PVTrack");
    REQUIRE(table.nRows() == 0);
    REQUIRE(table.columnIndex("vertexIdx") >= 0);
    REQUIRE(table.columnIndex("trackIdx") >= 0);
    REQUIRE(table.columnIndex("weight") >= 0);
  }

  SECTION("vertices without tracks, with a cut and maxLen") {
    edm::test::TestProcessor::Config cfg{config("process.toTest.cut = 'z > 0'\nprocess.toTest.maxLen = 2\n")};
    auto vertices = cfg.produces<reco::VertexCollection>("vertices");
    auto tracks = cfg.produces<reco::TrackCollection>("tracks");
    edm::test::TestProcessor tester(cfg);
    auto input = std::make_unique<reco::VertexCollection>();
    for (double z : {-1., 1., 2., 3.})
      input->push_back(vertex(z));
    auto event = tester.test(std::make_pair(vertices, std::move(input)),
                             std::make_pair(tracks, std::make_unique<reco::TrackCollection>()));
    REQUIRE(event.get<nanoaod::FlatTable>()->nRows() == 0);
  }

  SECTION("a track that is not in 'tracks'") {
    edm::test::TestProcessor::Config cfg{config()};
    auto vertices = cfg.produces<reco::VertexCollection>("vertices");
    auto tracks = cfg.produces<reco::TrackCollection>("tracks");
    edm::test::TestProcessor tester(cfg);
    reco::TrackCollection other(1);
    auto input = std::make_unique<reco::VertexCollection>(1, vertex(0.));
    input->front().add(reco::TrackBaseRef(reco::TrackRef(&other, 0)), 0.5f);
    REQUIRE_THROWS_AS(tester.test(std::make_pair(vertices, std::move(input)),
                                  std::make_pair(tracks, std::make_unique<reco::TrackCollection>(1))),
                      cms::Exception);
  }

  SECTION("missing vertices") {
    edm::test::TestProcessor::Config cfg{config()};
    edm::test::TestProcessor tester(cfg);
    REQUIRE_THROWS_AS(tester.test(), cms::Exception);
  }

  SECTION("missing vertices, skipNonExistingSrc") {
    edm::test::TestProcessor::Config cfg{config("process.toTest.skipNonExistingSrc = True\n")};
    edm::test::TestProcessor tester(cfg);
    REQUIRE(tester.test().get<nanoaod::FlatTable>()->nRows() == 0);
  }
}

// Vertices with tracks, made by VertexTrackTestInputProducer (this directory).
namespace {
  nanoaod::FlatTable links(std::string const& parameters = "", std::string const& vertices = "'input'") {
    edm::test::TestProcessor tester(edm::test::TestProcessor::Config{
        R"(from FWCore.TestProcessor.TestProcess import *
from PhysicsTools.NanoAOD.vertexTrackTableProducer_cfi import vertexTrackTableProducer
process = TestProcess()
process.input = cms.EDProducer('VertexTrackTestInputProducer')
process.toTest = vertexTrackTableProducer.clone(src = cms.InputTag()" +
        vertices + R"(), tracks = cms.InputTag('input', 'tracks'), name = 'PVTrack')" + parameters + R"()
process.moduleToTest(process.toTest, cms.Task(process.input))
)"});
    return *tester.test().get<nanoaod::FlatTable>();
  }

  template <typename T>
  std::vector<T> column(nanoaod::FlatTable const& table, std::string const& name) {
    REQUIRE(table.columnIndex(name) >= 0);
    auto const data = table.columnData<T>(table.columnIndex(name));
    return {data.begin(), data.end()};
  }
}  // namespace

TEST_CASE("VertexTrackTableProducer: rows of the links", "[VertexTrackTableProducer]") {
  SECTION("all vertices") {
    auto const table = links();
    REQUIRE(column<int32_t>(table, "vertexIdx") == std::vector<int32_t>{0, 0, 1, 2, 2});
    REQUIRE(column<int32_t>(table, "trackIdx") == std::vector<int32_t>{0, 2, 1, 2, 0});
    auto const weights = column<float>(table, "weight");
    const std::vector<float> expected{1.f, 0.5f, 0.25f, 1.f, 0.75f};
    REQUIRE(weights.size() == expected.size());
    for (std::size_t i = 0; i < weights.size(); ++i)
      REQUIRE(weights[i] == Catch::Approx(expected[i]).margin(1. / 255));  // stored in steps of 1/255
  }

  SECTION("vertexIdx is the row of the vertex after the cut") {
    auto const table = links(", cut = 'z > 0'");
    REQUIRE(column<int32_t>(table, "vertexIdx") == std::vector<int32_t>{0, 1, 1});
    REQUIRE(column<int32_t>(table, "trackIdx") == std::vector<int32_t>{1, 2, 0});
  }

  SECTION("maxLen") {
    auto const table = links(", cut = 'z > 0', maxLen = cms.uint32(1)");
    REQUIRE(column<int32_t>(table, "vertexIdx") == std::vector<int32_t>{0});
    REQUIRE(column<int32_t>(table, "trackIdx") == std::vector<int32_t>{1});
  }

  SECTION("a track of another collection") { REQUIRE_THROWS_AS(links("", "'input', 'otherVertices'"), cms::Exception); }
}
