#include "catch2/catch_all.hpp"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/TestProcessor/interface/TestProcessor.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace {
  std::string config(std::string const& extra = "") {
    return R"(from FWCore.TestProcessor.TestProcess import *
from DPGAnalysis.MuonTools.muDTSegmentMatchTableProducer_cfi import muDTSegmentMatchTableProducer
process = TestProcess()
process.toTest = muDTSegmentMatchTableProducer.clone(src = "muons", segmentSrc = "segments")
)" + extra +
           "process.moduleToTest(process.toTest)\n";
  }

  reco::Muon muon(double pt) { return reco::Muon(1, reco::Muon::LorentzVector(pt, 0, 0, pt)); }
}  // namespace

TEST_CASE("MuDTSegmentMatchTableProducer", "[MuSegmentMatchTableProducer]") {
  SECTION("muons without matches") {
    edm::test::TestProcessor::Config cfg{config()};
    auto muons = cfg.produces<reco::MuonCollection>("muons");
    auto segments = cfg.produces<DTRecSegment4DCollection>("segments");
    edm::test::TestProcessor tester(cfg);
    auto event = tester.test(std::make_pair(muons, std::make_unique<reco::MuonCollection>(3, muon(10.))),
                             std::make_pair(segments, std::make_unique<DTRecSegment4DCollection>()));
    auto const& table = *event.get<nanoaod::FlatTable>();
    REQUIRE(table.name() == "muonDTSegmentMatch");
    REQUIRE(table.nRows() == 0);
    for (auto const* column : {"muonIdx", "segmentIdx", "mask"})
      REQUIRE(table.columnIndex(column) >= 0);
  }

  SECTION("a segment that is not in 'segmentSrc'") {
    edm::test::TestProcessor::Config cfg{config()};
    auto muons = cfg.produces<reco::MuonCollection>("muons");
    auto segments = cfg.produces<DTRecSegment4DCollection>("segments");
    edm::test::TestProcessor tester(cfg);
    DTRecSegment4DCollection other;
    reco::MuonSegmentMatch match;
    match.dtSegmentRef = DTRecSegment4DRef(&other, 0);
    reco::MuonChamberMatch chamber;
    chamber.id = DTChamberId(0, 1, 1);
    chamber.segmentMatches.push_back(match);
    auto input = std::make_unique<reco::MuonCollection>(1, muon(10.));
    input->front().setMatches({chamber});
    REQUIRE_THROWS_AS(tester.test(std::make_pair(muons, std::move(input)),
                                  std::make_pair(segments, std::make_unique<DTRecSegment4DCollection>())),
                      cms::Exception);
  }

  SECTION("a segment match of a muon that the cut removes") {
    edm::test::TestProcessor::Config cfg{config("process.toTest.cut = 'pt > 20'\n")};
    auto muons = cfg.produces<reco::MuonCollection>("muons");
    auto segments = cfg.produces<DTRecSegment4DCollection>("segments");
    edm::test::TestProcessor tester(cfg);
    DTRecSegment4DCollection other;
    reco::MuonSegmentMatch match;
    match.dtSegmentRef = DTRecSegment4DRef(&other, 0);
    reco::MuonChamberMatch chamber;
    chamber.id = DTChamberId(0, 1, 1);
    chamber.segmentMatches.push_back(match);
    auto input = std::make_unique<reco::MuonCollection>(1, muon(10.));
    input->front().setMatches({chamber});
    auto event = tester.test(std::make_pair(muons, std::move(input)),
                             std::make_pair(segments, std::make_unique<DTRecSegment4DCollection>()));
    REQUIRE(event.get<nanoaod::FlatTable>()->nRows() == 0);
  }

  SECTION("missing segments") {
    edm::test::TestProcessor::Config cfg{config()};
    auto muons = cfg.produces<reco::MuonCollection>("muons");
    edm::test::TestProcessor tester(cfg);
    REQUIRE_THROWS_AS(tester.test(std::make_pair(muons, std::make_unique<reco::MuonCollection>(1, muon(10.)))),
                      cms::Exception);
  }

  SECTION("missing segments, skipNonExistingSrc") {
    edm::test::TestProcessor::Config cfg{config("process.toTest.skipNonExistingSrc = True\n")};
    auto muons = cfg.produces<reco::MuonCollection>("muons");
    edm::test::TestProcessor tester(cfg);
    auto event = tester.test(std::make_pair(muons, std::make_unique<reco::MuonCollection>(1, muon(10.))));
    REQUIRE(event.get<nanoaod::FlatTable>()->nRows() == 0);
  }
}

// Segments and muons with matches to them, made by MuSegmentMatchTestInputProducer (this directory).
namespace {
  std::string matchConfig(std::string const& plugin,
                          std::string const& parameters,
                          std::string const& muons = "'input'") {
    return R"(from FWCore.TestProcessor.TestProcess import *
process = TestProcess()
process.input = cms.EDProducer('MuSegmentMatchTestInputProducer')
process.toTest = cms.EDProducer(')" +
           plugin + "', src = cms.InputTag(" + muons + "), " + parameters + R"()
process.moduleToTest(process.toTest, cms.Task(process.input))
)";
  }

  template <typename T>
  std::vector<T> column(nanoaod::FlatTable const& table, std::string const& name) {
    REQUIRE(table.columnIndex(name) >= 0);
    auto const data = table.columnData<T>(table.columnIndex(name));
    return {data.begin(), data.end()};
  }

  nanoaod::FlatTable matches(std::string const& plugin, std::string const& parameters) {
    edm::test::TestProcessor tester(edm::test::TestProcessor::Config{matchConfig(plugin, parameters)});
    return *tester.test().get<nanoaod::FlatTable>();
  }
}  // namespace

TEST_CASE("MuSegmentMatchTableProducer: rows of the matches", "[MuSegmentMatchTableProducer]") {
  SECTION("DT") {
    auto const table = matches("MuDTSegmentMatchTableProducer", "segmentSrc = cms.InputTag('input', 'dtSegments')");
    REQUIRE(table.name() == "muonDTSegmentMatch");
    REQUIRE(column<int32_t>(table, "muonIdx") == std::vector<int32_t>{0, 0, 2});
    REQUIRE(column<int32_t>(table, "segmentIdx") == std::vector<int32_t>{1, 0, 2});
    REQUIRE(column<uint32_t>(table, "mask") == std::vector<uint32_t>{1 << 10, 0, 1 << 18});
  }

  SECTION("CSC: a segment match without a segment is left out") {
    auto const table = matches("MuCSCSegmentMatchTableProducer", "segmentSrc = cms.InputTag('input', 'cscSegments')");
    REQUIRE(column<int32_t>(table, "muonIdx") == std::vector<int32_t>{0, 2});
    REQUIRE(column<int32_t>(table, "segmentIdx") == std::vector<int32_t>{1, 0});
    REQUIRE(column<uint32_t>(table, "mask") == std::vector<uint32_t>{2, 1});
  }

  SECTION("GEM") {
    auto const table = matches("MuGEMSegmentMatchTableProducer", "segmentSrc = cms.InputTag('input', 'gemSegments')");
    REQUIRE(column<int32_t>(table, "muonIdx") == std::vector<int32_t>{0});
    REQUIRE(column<int32_t>(table, "segmentIdx") == std::vector<int32_t>{0});
    REQUIRE(column<uint32_t>(table, "mask") == std::vector<uint32_t>{4});
  }

  SECTION("muonIdx is the row of the muon after the cut") {
    auto const table = matches("MuDTSegmentMatchTableProducer",
                               "segmentSrc = cms.InputTag('input', 'dtSegments'), cut = cms.string('pt > 15')");
    REQUIRE(column<int32_t>(table, "muonIdx") == std::vector<int32_t>{1});
    REQUIRE(column<int32_t>(table, "segmentIdx") == std::vector<int32_t>{2});
  }

  SECTION("maxLen") {
    auto const table = matches(
        "MuDTSegmentMatchTableProducer",
        "segmentSrc = cms.InputTag('input', 'dtSegments'), cut = cms.string('pt > 15'), maxLen = cms.uint32(1)");
    REQUIRE(table.nRows() == 0);
  }

  SECTION("a segment of another collection") {
    edm::test::TestProcessor tester(edm::test::TestProcessor::Config{matchConfig(
        "MuDTSegmentMatchTableProducer", "segmentSrc = cms.InputTag('input', 'dtSegments')", "'input', 'otherMuons'")});
    REQUIRE_THROWS_AS(tester.test(), cms::Exception);
  }
}
