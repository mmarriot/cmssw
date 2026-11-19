#include "catch2/catch_all.hpp"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/TestProcessor/interface/TestProcessor.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <cstdint>
#include <string>
#include <vector>

// The input is made by PFCandidateLinkTestInputProducer (this directory): six candidates of pt 10 to 60, with
// tracks 0, -, 2, -, -, 1; GSF tracks -, 1, -, -, -, -; muons -, -, 1, -, -, -; the last two copies of the TICL PF
// candidates; block elements of the clusters ecal 2 and hcal 0 (candidate 0), ecal 0 (1), other 0 (2).
namespace {
  std::string configuration(std::string const& parameters = "") {
    return R"cfg(from FWCore.TestProcessor.TestProcess import *
process = TestProcess()
process.input = cms.EDProducer('PFCandidateLinkTestInputProducer')
process.toTest = cms.EDProducer('PFCandidateLinkTableProducer', src = cms.InputTag('input'), name = cms.string('PFCand'),
    tracks = cms.InputTag('input', 'tracks'), gsfTracks = cms.InputTag('input', 'gsfTracks'),
    muons = cms.InputTag('input', 'muons'), ticlCandidates = cms.InputTag('input', 'ticl'),
    blocks = cms.InputTag('input', 'blocks'),
    clusters = cms.VPSet(
        cms.PSet(src = cms.InputTag('input', 'ecal'), name = cms.string('PFCandECAL'), doc = cms.string('ECAL')),
        cms.PSet(src = cms.InputTag('input', 'hcal'), name = cms.string('PFCandHCAL'))))
)cfg" + parameters +
           R"cfg(
process.moduleToTest(process.toTest, cms.Task(process.input))
)cfg";
  }

  std::vector<int32_t> column(nanoaod::FlatTable const& table, std::string const& name) {
    REQUIRE(table.columnIndex(name) >= 0);
    auto const data = table.columnData<int32_t>(table.columnIndex(name));
    return {data.begin(), data.end()};
  }
}  // namespace

TEST_CASE("PFCandidateLinkTableProducer: rows of the objects of each candidate", "[PFCandidateLinkTableProducer]") {
  edm::test::TestProcessor tester(edm::test::TestProcessor::Config{configuration()});
  auto const event = tester.test();

  auto const& table = *event.get<nanoaod::FlatTable>();
  REQUIRE(table.name() == "PFCand");
  REQUIRE(table.extension());
  REQUIRE(table.nRows() == 6);
  REQUIRE(table.nColumns() == 4);
  REQUIRE(column(table, "trackIdx") == std::vector<int32_t>{0, -1, 2, -1, -1, 1});
  REQUIRE(column(table, "gsfTrackIdx") == std::vector<int32_t>{-1, 1, -1, -1, -1, -1});
  REQUIRE(column(table, "muonIdx") == std::vector<int32_t>{-1, -1, 1, -1, -1, -1});
  REQUIRE(column(table, "ticlCandidateIdx") == std::vector<int32_t>{-1, -1, -1, -1, 0, 1});

  // the element of the "other" clusters, not configured, and the null cluster are left out
  auto const& ecal = *event.get<nanoaod::FlatTable>("PFCandECAL");
  REQUIRE(ecal.name() == "PFCandECAL");
  REQUIRE(ecal.doc() == "ECAL");
  REQUIRE(column(ecal, "pfCandIdx") == std::vector<int32_t>{0, 1});
  REQUIRE(column(ecal, "clusterIdx") == std::vector<int32_t>{2, 0});
  auto const& hcal = *event.get<nanoaod::FlatTable>("PFCandHCAL");
  REQUIRE(column(hcal, "pfCandIdx") == std::vector<int32_t>{0});
  REQUIRE(column(hcal, "clusterIdx") == std::vector<int32_t>{0});
}

TEST_CASE("PFCandidateLinkTableProducer: the selection of the candidate table", "[PFCandidateLinkTableProducer]") {
  auto rows = [](std::string const& parameters) {
    edm::test::TestProcessor tester(edm::test::TestProcessor::Config{configuration(parameters)});
    auto const event = tester.test();
    return std::make_tuple(column(*event.get<nanoaod::FlatTable>(), "trackIdx"),
                           column(*event.get<nanoaod::FlatTable>(), "ticlCandidateIdx"),
                           column(*event.get<nanoaod::FlatTable>("PFCandECAL"), "pfCandIdx"));
  };
  SECTION("cut: the candidate index of the cluster pairs is the row among the selected candidates") {
    auto const [tracks, ticl, ecal] = rows("process.toTest.cut = cms.string('pt > 15')");
    REQUIRE(tracks == std::vector<int32_t>{-1, 2, -1, -1, 1});
    REQUIRE(ticl == std::vector<int32_t>{-1, -1, -1, 0, 1});
    REQUIRE(ecal == std::vector<int32_t>{0});
  }
  SECTION("maxLen") {
    auto const [tracks, ticl, ecal] = rows("process.toTest.maxLen = cms.uint32(2)");
    REQUIRE(tracks == std::vector<int32_t>{0, -1});
    REQUIRE(ticl == std::vector<int32_t>{-1, -1});
    REQUIRE(ecal == std::vector<int32_t>{0, 1});
  }
}

TEST_CASE("PFCandidateLinkTableProducer: optional inputs", "[PFCandidateLinkTableProducer]") {
  edm::test::TestProcessor tester(edm::test::TestProcessor::Config{
      configuration("process.toTest.tracks = ''\nprocess.toTest.gsfTracks = ''\nprocess.toTest.muons = ''\n"
                    "process.toTest.ticlCandidates = ''\nprocess.toTest.clusters = []")});
  auto const event = tester.test();
  auto const& table = *event.get<nanoaod::FlatTable>();
  REQUIRE(table.nRows() == 6);
  REQUIRE(table.nColumns() == 0);
}

TEST_CASE("PFCandidateLinkTableProducer: inputs that do not match", "[PFCandidateLinkTableProducer]") {
  SECTION("references into another collection than the configured one") {
    edm::test::TestProcessor tester(
        edm::test::TestProcessor::Config{configuration("process.toTest.muons = cms.InputTag('input', 'otherMuons')")});
    REQUIRE_THROWS_AS(tester.test(), cms::Exception);
  }
  SECTION("candidates that do not end with the TICL candidates") {
    edm::test::TestProcessor tester(edm::test::TestProcessor::Config{
        configuration("process.toTest.ticlCandidates = cms.InputTag('input', 'reversedTicl')")});
    REQUIRE_THROWS_AS(tester.test(), cms::Exception);
  }
  SECTION("blocks that are not those of the candidates") {
    edm::test::TestProcessor tester(edm::test::TestProcessor::Config{
        configuration("process.toTest.blocks = cms.InputTag('input', 'otherBlocks')")});
    REQUIRE_THROWS_AS(tester.test(), cms::Exception);
  }
  SECTION("a cluster table named as the candidate table") {
    REQUIRE_THROWS_AS(edm::test::TestProcessor(edm::test::TestProcessor::Config{
                          configuration("process.toTest.clusters[0].name = 'PFCand'")}),
                      cms::Exception);
  }
}
