#include "catch2/catch_all.hpp"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/TestProcessor/interface/TestProcessor.h"
#include "FWCore/Utilities/interface/Exception.h"

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
