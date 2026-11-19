#include "catch2/catch_all.hpp"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/TestProcessor/interface/TestProcessor.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <cstdint>
#include <string>
#include <vector>

// The input is made by MuonTrackLinkTestInputProducer (this directory): three muons with inner tracks general 1, 0, -;
// outer tracks standAloneAtVertex 1, standAlone 0, -; global tracks global 0, -, -.
namespace {
  std::string configuration(std::string const& parameters = "") {
    return R"cfg(from FWCore.TestProcessor.TestProcess import *
process = TestProcess()
process.input = cms.EDProducer('MuonTrackLinkTestInputProducer')
process.toTest = cms.EDProducer('MuonTrackLinkTableProducer', src = cms.InputTag('input', 'muons'),
    name = cms.string('Muon'),
    tracks = cms.VPSet(
        cms.PSet(name = cms.string('innerTrackIdx'), type = cms.string('inner'), src = cms.InputTag('input', 'general'),
                 doc = cms.string('inner track')),
        cms.PSet(name = cms.string('outerTrackIdx'), type = cms.string('outer'),
                 src = cms.InputTag('input', 'standAloneAtVertex')),
        cms.PSet(name = cms.string('outerTrackNotAtVertexIdx'), type = cms.string('outer'),
                 src = cms.InputTag('input', 'standAlone')),
        cms.PSet(name = cms.string('globalTrackIdx'), type = cms.string('global'), src = cms.InputTag('input', 'global'))))
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

TEST_CASE("MuonTrackLinkTableProducer: rows of the tracks of each muon", "[MuonTrackLinkTableProducer]") {
  edm::test::TestProcessor tester(edm::test::TestProcessor::Config{configuration()});
  auto const event = tester.test();
  auto const& table = *event.get<nanoaod::FlatTable>();
  REQUIRE(table.name() == "Muon");
  REQUIRE(table.extension());
  REQUIRE(table.nColumns() == 4);
  REQUIRE(table.columnDoc(table.columnIndex("innerTrackIdx")) == "inner track");
  REQUIRE(column(table, "innerTrackIdx") == std::vector<int32_t>{1, 0, -1});
  REQUIRE(column(table, "outerTrackIdx") == std::vector<int32_t>{1, -1, -1});
  REQUIRE(column(table, "outerTrackNotAtVertexIdx") == std::vector<int32_t>{-1, 0, -1});
  REQUIRE(column(table, "globalTrackIdx") == std::vector<int32_t>{0, -1, -1});
}

TEST_CASE("MuonTrackLinkTableProducer: the selection of the muon table", "[MuonTrackLinkTableProducer]") {
  edm::test::TestProcessor tester(
      edm::test::TestProcessor::Config{configuration("process.toTest.maxLen = cms.uint32(2)")});
  auto const event = tester.test();
  REQUIRE(column(*event.get<nanoaod::FlatTable>(), "innerTrackIdx") == std::vector<int32_t>{1, 0});
}

TEST_CASE("MuonTrackLinkTableProducer: tracks in collections that are not configured", "[MuonTrackLinkTableProducer]") {
  // by default -1 with a warning; with strict, an error
  const std::string strict = "process.toTest.strict = cms.bool(True)\n";
  SECTION("an outer track in neither outer-track collection") {
    const std::string other = "process.toTest.src = cms.InputTag('input', 'otherMuons')\n";
    edm::test::TestProcessor tester(edm::test::TestProcessor::Config{configuration(other)});
    auto const event = tester.test();
    REQUIRE(column(*event.get<nanoaod::FlatTable>(), "outerTrackIdx") == std::vector<int32_t>{-1});
    REQUIRE(column(*event.get<nanoaod::FlatTable>(), "outerTrackNotAtVertexIdx") == std::vector<int32_t>{-1});
    edm::test::TestProcessor strictTester(edm::test::TestProcessor::Config{configuration(other + strict)});
    REQUIRE_THROWS_AS(strictTester.test(), cms::Exception);
  }
  SECTION("only one of the two outer-track collections") {
    const std::string one = "del process.toTest.tracks[2]\n";
    edm::test::TestProcessor tester(edm::test::TestProcessor::Config{configuration(one)});
    auto const event = tester.test();
    REQUIRE(column(*event.get<nanoaod::FlatTable>(), "outerTrackIdx") == std::vector<int32_t>{1, -1, -1});
    edm::test::TestProcessor strictTester(edm::test::TestProcessor::Config{configuration(one + strict)});
    REQUIRE_THROWS_AS(strictTester.test(), cms::Exception);
  }
}

TEST_CASE("MuonTrackLinkTableProducer: missing muons", "[MuonTrackLinkTableProducer]") {
  const std::string missing = "process.toTest.src = cms.InputTag('input', 'missing')\n";
  SECTION("fail by default") {
    edm::test::TestProcessor tester(edm::test::TestProcessor::Config{configuration(missing)});
    REQUIRE_THROWS_AS(tester.test(), cms::Exception);
  }
  SECTION("an empty table with skipNonExistingSrc") {
    edm::test::TestProcessor tester(edm::test::TestProcessor::Config{
        configuration(missing + "process.toTest.skipNonExistingSrc = cms.bool(True)")});
    auto const event = tester.test();
    REQUIRE(event.get<nanoaod::FlatTable>()->nRows() == 0);
    REQUIRE(event.get<nanoaod::FlatTable>()->nColumns() == 4);
  }
}
